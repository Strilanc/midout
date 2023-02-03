import collections
from typing import Dict, Optional, List, Tuple, DefaultDict, TYPE_CHECKING, Iterable

if TYPE_CHECKING:
    from midout.gen import Builder

# These rotations form a group.
# This is the multiplication table of the group.
# We use the stim operation names of the elements.
XZ_BASE_TRANSITIONS = [
    ['I', 'H_YZ', 'H_XY', 'C_ZYX', 'C_XYZ', 'H'],
    ['H_YZ', 'I', 'C_ZYX', 'H_XY', 'H', 'C_XYZ'],
    ['H_XY', 'C_XYZ', 'I', 'H', 'H_YZ', 'C_ZYX'],
    ['C_XYZ', 'H_XY', 'H', 'I', 'C_ZYX', 'H_YZ'],
    ['C_ZYX', 'H', 'H_YZ', 'C_XYZ', 'I', 'H_XY'],
    ['H', 'C_ZYX', 'C_XYZ', 'H_YZ', 'H_XY', 'I'],
]


SINGLE_QUBIT_MAP = {
    (XZ_BASE_TRANSITIONS[u][0], XZ_BASE_TRANSITIONS[0][v]): XZ_BASE_TRANSITIONS[u][v]
    for u in range(6)
    for v in range(6)
}
SINGLE_QUBIT_EQUIV_MAP = {
    'S': 'H_XY',
    'S_DAG': 'H_XY',
    'SQRT_Y': 'H',
    'SQRT_Y_DAG': 'H',
    'SQRT_X': 'H_YZ',
    'SQRT_X_DAG': 'H_YZ',
}


def _compute_xy_base_transition_map() -> Dict[str, str]:
    # There are 6 single qubit stabilizer rotations (up to Paulis).
    # Here we name them based on where they map the X and Z axes.
    # The first character says where the X axis is mapped, and the
    # second character says where the Z axis is mapped.
    # For example, in "YX" the first character is 'Y' so an X axis
    # input is mapped to a Y axis output and the second character
    # is 'X' so a Z axis input is mapped to an X axis output.
    XZ_BASES = ['XY', 'XZ', 'YX', 'YZ', 'ZX', 'ZY']

    # For more convenient lookup we store the multiplication table into a dictionary.
    result = {
        f'{XZ_BASES[inp]} -> {XZ_BASES[out]}': XZ_BASE_TRANSITIONS[out][inp]
        for inp in range(len(XZ_BASES))
        for out in range(len(XZ_BASES))
    }

    return result


XZ_BASE_TRANSITION_MAP: Dict[str, str] = _compute_xy_base_transition_map()
DESIRED_Z_TO_ORIENTATION: Dict[str, str] = {
    'X': 'ZX',
    'Y': 'ZY',
    'Z': 'XZ',
}


class SingleQubitGatesPlanner():
    def __init__(self, builder: 'Builder'):
        self.builder = builder
        self.gates = {}

    def gate(self, gate: str, qubits: Iterable[complex]):
        gate = SINGLE_QUBIT_EQUIV_MAP.get(gate, gate)
        for q in qubits:
            prev = self.gates.get(q)
            if prev is None:
                self.gates[q] = gate
            else:
                combined = SINGLE_QUBIT_MAP[(gate, prev)]
                if combined == 'I':
                    del self.gates[q]
                else:
                    self.gates[q] = combined

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        gate_groups = collections.defaultdict(list)
        for q, g in self.gates.items():
            gate_groups[g].append(q)
        for g in sorted(gate_groups.keys()):
            targets = gate_groups[g]
            if g == "H_XY":
                g = "S"
            self.builder.gate(g, targets)
        if gate_groups:
            self.builder.tick()


class InteractionPlanner:
    """This class is responsible for compiling a series of Pauli interactions.

    All PCP interactions are supported (CX, CY, CZ, XCX, etc.), and converted
    into a circuit that only performs CZ gates separates by single qubit basis
    changes.

    Normally this class is used via a Builder:

        with builder.plan_interactions(layer_count=..., ...) as planner:
            ...

    The desired interactions are added to the planner during the `with` block,
    and when the `with` block ends the requested interactions are compiled into
    operations appended to the builder's circuit.

    The compilation works by tracking the required orientation of each qubit
    over time. The gates to switch from one orientation to the other are hardcoded
    above in XZ_BASE_TRANSITION_MAP. The orientations are changed as needed in order
    for the output CZ operation to perform the desired input PCP interaction in
    each layer.

    The resulting circuit will have CZ layers separated by (optional) single
    qubit basis change layers. Sometimes a single qubit basis layer can be
    removed (e.g. a qubit needs to change orientations, but it was unused in the
    previous CZ layer so the basis change can happen in an earlier basis change
    layer, and moving it out empties the current basis change layer). The code
    has a few heuristics to find these situations and remove the unnecessary
    basis change layers, though it is by no means guaranteed to find the minimum
    solution.
    """

    def __init__(self,
                 *,
                 layer_count: int,
                 builder: 'Builder',
                 start_orientations: Optional[Dict[complex, str]],
                 end_orientations: Optional[Dict[complex, str]]):
        self.builder = builder
        self.start_orientations = start_orientations
        self.end_orientations = end_orientations
        self.layer_interactions: List[List[Tuple[complex, complex]]] = []
        self.layer_orientations: List[Dict[complex, str]] = []
        for _ in range(layer_count):
            self.layer_interactions.append([])
            self.layer_orientations.append({})

    def _do_switch_basis(self,
                         *,
                         old_orientations: Dict[complex, str],
                         partial_new_orientations: Dict[complex, str],
                         pre_tick_if_any: bool,
                         post_tick_if_any: bool) -> Dict[complex, str]:
        groups: DefaultDict[str, List[complex]] = collections.defaultdict(list)
        for q, old_b in partial_new_orientations.items():
            new_b = old_orientations.get(q, 'XZ')
            transition = XZ_BASE_TRANSITION_MAP[f'{old_b} -> {new_b}']
            groups[transition].append(q)
        if 'I' in groups:
            del groups['I']
        if groups:
            if pre_tick_if_any:
                self.builder.tick()
            for gate in sorted(groups.keys()):
                self.builder.gate(gate, groups[gate])
            if post_tick_if_any:
                self.builder.tick()
        return {**old_orientations, **partial_new_orientations}

    def pcp(self, b1: str, b2: str, q1: complex, q2: complex, *, layer: int) -> None:
        layer_orientation = self.layer_orientations[layer]
        if q1 in layer_orientation:
            raise ValueError(f"Collision at {q1}")
        if q2 in layer_orientation:
            raise ValueError(f"Collision at {q2}")
        layer_orientation[q1] = DESIRED_Z_TO_ORIENTATION[b1]
        layer_orientation[q2] = DESIRED_Z_TO_ORIENTATION[b2]
        self.layer_interactions[layer].append((q1, q2))

    def cx(self, q1: complex, q2: complex, *, layer: int) -> None:
        self.pcp('Z', 'X', q1=q1, q2=q2, layer=layer)

    def cz(self, q1: complex, q2: complex, *, layer: int) -> None:
        self.pcp('Z', 'Z', q1=q1, q2=q2, layer=layer)

    @staticmethod
    def simplify_layer_base_transitions(basis_layers: List[Dict[complex, str]]):
        basis_layers = [dict(e) for e in basis_layers]
        for k in range(len(basis_layers)):
            cur = basis_layers[k]
            prev = basis_layers[k - 1] if k > 0 else None
            if prev and all(prev.get(q, b) == b for q, b in cur.items()):
                for q, v in cur.items():
                    prev[q] = v
                cur.clear()
        return basis_layers

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        start_orientations = dict(self.start_orientations or {})
        end_orientations = dict(self.end_orientations or {})
        for q in self.builder.q2i.keys():
            if q not in start_orientations:
                start_orientations[q] = 'XZ'
            if q not in end_orientations:
                end_orientations[q] = 'XZ'

        num_layers = len(self.layer_orientations)
        layer_orientations = InteractionPlanner.simplify_layer_base_transitions([
            *self.layer_orientations,
            end_orientations,
        ])

        full_cur_orientations = dict(start_orientations)
        for k in range(num_layers):
            full_cur_orientations = self._do_switch_basis(
                old_orientations=full_cur_orientations,
                partial_new_orientations=layer_orientations[k],
                post_tick_if_any=True,
                pre_tick_if_any=False,
            )
            self.builder.cz(self.layer_interactions[k])
            if k < num_layers - 1:
                self.builder.tick()
        self._do_switch_basis(
            old_orientations=full_cur_orientations,
            partial_new_orientations=layer_orientations[-1],
            post_tick_if_any=False,
            pre_tick_if_any=True,
        )
