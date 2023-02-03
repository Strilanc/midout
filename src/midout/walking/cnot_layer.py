import dataclasses
from typing import cast, Dict, Iterable, Set, Tuple

import matplotlib as mpl
import matplotlib.pyplot as plt

from midout.walking.diamonds import Diamond, DiamondState
from midout.walking.util import Basis, Direction, Qubit, qubit_coords

CNOT_CIRCLE_RADIUS = 0.09


def _make_cnot_layer(
    qubits: Set[Qubit], center: Qubit, direction: Direction, aligned=True, flip_cnots=False
) -> Set[Tuple[Qubit, Qubit]]:
    """returns all appropriate CNOTS over a given set of qubits.

    CNOTs are generated in the following way:
        consider an integer square grid, with even (e) and odd (o) subgrids
            c o e o
            o e o e
            e o e o
        with the even subgrid defined by including the 'center' qubit (c)
        for each qubit on this grid relative to center, try to make a cnot in direction
            the cnot is included if the qubit the cnot is to is also in qubits


    first, there is a CNOT between the origin qubit and the qubit origin+direction
        The origin qubit will be the control. This means gates will be in the right direction
        for the surface code if the origin is the measure qubit of a Z plaquette

    Qubits on the same integer grid as the origin will also have CNOTs
        the 'parity' of those qubits on the grid determine if they're the target or control
        ie qubits immediately neighbouring the origin on the same integer grid will be targets

    Args:
        qubits: all qubits involved in the layer, all returned cnots will involve only these qubits
        center: the qubit to center the layer on
        direction: the direction of the cnot involving the center qubit
            should not point to a qubit on the integer grid relative to the center qubit,
        aligned: whether to align or anti-align the cnots
            if true, all cnots on the even subgrid will point in the same direction
            if false, cnots on the odd subgrid not including the origin will reverse their direction
        flip_cnots: whether to reverse the control and target of every cnot in the layer
    """
    center = cast(complex, center)
    direction = cast(complex, direction)
    if direction.real.is_integer() and direction.imag.is_integer():
        raise ValueError
    cnot_pairs = []
    for q in qubits:
        q = cast(complex, q)
        direction_from_center = q - center
        if direction_from_center.real.is_integer() and direction_from_center.imag.is_integer():
            # q is on the integer subgrid, so it should produce a cnot
            odd_subgrid = (direction_from_center.real % 2) != (direction_from_center.imag % 2)
            reverse_cnot = odd_subgrid != flip_cnots  # XOR is equiv to != for bools
            this_direction = -1 * direction if odd_subgrid and not aligned else direction
            other_cnot_qubit = q + this_direction
            if other_cnot_qubit not in qubits:
                continue
            gate_qubits = (q, q + this_direction)
            if reverse_cnot:
                gate_qubits = gate_qubits[::-1]

            cnot_pairs.append(gate_qubits)

    return cast(Set[Tuple[Qubit, Qubit]], set(cnot_pairs))


@dataclasses.dataclass
class CnotLayer:
    targets: Dict[Qubit, Qubit]
    controls: Dict[Qubit, Qubit]

    @property
    def all_qubits(self) -> Iterable[Qubit]:
        return list(self.targets.keys()) + list(self.controls.keys())

    def __contains__(self, item):
        return item in self.targets or item in self.controls

    def add_pair(self, control, target):
        if control in self or target in self:
            raise ValueError("A qubit can't be involved in 2 CNOTs in the same layer. ")
        self.controls[control] = target
        self.targets[target] = control
        return self

    def remove_gates_touching_qubit(self, qubit: Qubit, quiet=False):
        if qubit in self.targets.keys():
            other = self.targets.pop(qubit)
            self.controls.pop(other)
        elif qubit in self.controls.keys():
            other = self.controls.pop(qubit)
            self.targets.pop(other)
        elif not quiet:
            raise ValueError(f"qubit {qubit} not involved in any gates")

    @staticmethod
    def from_tuples(tuples: Iterable[Tuple[Qubit, Qubit]]):
        """construct the class from tuples of (control, target)"""
        layer = CnotLayer(targets={}, controls={})
        for control, target in tuples:
            layer.add_pair(control=control, target=target)
        return layer

    @staticmethod
    def from_layer_spec(*args, **kwargs):
        return CnotLayer.from_tuples(_make_cnot_layer(*args, **kwargs))

    @staticmethod
    def for_standard_surface_code_cycle(
        qubits, center=1 + 1j, arrangement='standard'
    ) -> Tuple['CnotLayer', 'CnotLayer', 'CnotLayer', 'CnotLayer']:
        """for the standard surface code round."""

        if arrangement == 'standard':
            # this one returns the diamonds that craig originally posted
            return (
                CnotLayer.from_layer_spec(qubits, center, direction=0.5 + 0.5j, aligned=True),
                CnotLayer.from_layer_spec(qubits, center, direction=-0.5 + 0.5j, aligned=False),
                CnotLayer.from_layer_spec(qubits, center, direction=0.5 - 0.5j, aligned=False),
                CnotLayer.from_layer_spec(qubits, center, direction=-0.5 - 0.5j, aligned=True),
            )
        if arrangement == 'hook_errors':
            # this ones flips the middle two directions
            # which gives you bad hook errors in both observables
            return (
                CnotLayer.from_layer_spec(qubits, center, direction=0.5 + 0.5j, aligned=True),
                CnotLayer.from_layer_spec(qubits, center, direction=0.5 - 0.5j, aligned=False),
                CnotLayer.from_layer_spec(qubits, center, direction=-0.5 + 0.5j, aligned=False),
                CnotLayer.from_layer_spec(qubits, center, direction=-0.5 - 0.5j, aligned=True),
            )
        elif arrangement == 'all_together':
            # this one has all 4 layers aligned, so basically, SS instead of SN
            # this gets hook errors on one direction wrong, but the other direction right
            return (
                CnotLayer.from_layer_spec(qubits, center, direction=0.5 + 0.5j, aligned=True),
                CnotLayer.from_layer_spec(qubits, center, direction=0.5 - 0.5j, aligned=True),
                CnotLayer.from_layer_spec(qubits, center, direction=-0.5 + 0.5j, aligned=True),
                CnotLayer.from_layer_spec(qubits, center, direction=-0.5 - 0.5j, aligned=True),
            )
        elif arrangement == 'stupid':
            # this ones does the checks in a circle, interlacing the stabilizers, CC rather than SN
            return (
                CnotLayer.from_layer_spec(qubits, center, direction=0.5 + 0.5j, aligned=True),
                CnotLayer.from_layer_spec(qubits, center, direction=0.5 - 0.5j, aligned=True),
                CnotLayer.from_layer_spec(qubits, center, direction=-0.5 - 0.5j, aligned=True),
                CnotLayer.from_layer_spec(qubits, center, direction=-0.5 + 0.5j, aligned=True),
            )
        elif arrangement == 'step':
            # this ones moves the center halfway through the cycle
            # I think it's important that the first two layers are the 'gets hook errors wrong' ones
            # so that the last two layers can 'get hook errors right'
            return (
                CnotLayer.from_layer_spec(qubits, center, direction=0.5 + 0.5j, aligned=True),
                CnotLayer.from_layer_spec(qubits, center, direction=0.5 - 0.5j, aligned=False),
                CnotLayer.from_layer_spec(
                    qubits, center + 0.5 + 0.5j, direction=-0.5 + 0.5j, aligned=False
                ),
                CnotLayer.from_layer_spec(
                    qubits, center + 0.5 + 0.5j, direction=0.5 + 0.5j, aligned=True
                ),
            )
        else:
            raise ValueError

    def _check_consistency(self):
        for q in self.controls.keys():
            assert q == self.targets[self.controls[q]]
        for q in self.targets.keys():
            assert q == self.controls[self.targets[q]]

    def source_and_sink(self, basis: Basis):
        """return sources and sinks mappings.

        for a given basis,
        a stabilizer on a source qubit (q) will spread to the sink (sources[q])
            unless it's already in the set, in which case the sink gets dropped
        a stabilizer on a sink qubit will go straight through,
            unless the source is in the set, in which case it gets dropped

        for the Z basis, it's like we're tracking X errors, so source is control, sink is target
        """
        if basis == Basis.z:
            return self.controls, self.targets
        else:
            return self.targets, self.controls

    def track_stabilizers_dict(
        self, qubits: Iterable[Qubit], basis: Basis
    ) -> Dict[Qubit, Iterable[Qubit]]:
        """tracks a stabilizer through the cnot layer, like propagating errors."""
        sources, sinks = self.source_and_sink(basis)
        output_qubits = {}
        for q in qubits:
            if q in sources and sources[q] not in qubits:
                output_qubits[q] = [q, sources[q]]
            elif q in sinks and sinks[q] in qubits:
                output_qubits[q] = [sinks[q]]
            else:
                output_qubits[q] = [q]
        return output_qubits

    def track_stabilizers(self, *args, **kwargs) -> Iterable[Qubit]:
        """tracks a stabilizer through the cnot layer, like propagating errors."""
        output_map = self.track_stabilizers_dict(*args, **kwargs)
        return set(q for qubits in output_map.values() for q in qubits)

    def track_diamond(self, diamond: Diamond):
        qubits = self.track_stabilizers(qubits=diamond.qubits, basis=diamond.basis)
        return Diamond(
            qubits=frozenset(qubits),
            basis=diamond.basis,
            marker_type=diamond.marker_type,
            override_color=diamond.override_color,
            measurements=diamond.measurements,
            qubit_to_reinclude=diamond.qubit_to_reinclude,
            duid=diamond.duid,
        )

    def track_diamond_state(self, diamond_state: DiamondState):
        """make a new DiamondState that results from this one being put through a cnot layer."""
        diamonds = [self.track_diamond(d) for d in diamond_state.diamonds]
        observables = {
            b: set(self.track_stabilizers(qubits=obvs, basis=b))
            for b, obvs in diamond_state.observables.items()
        }
        return DiamondState(diamonds=frozenset(diamonds), observables=observables)

    def plot(self, ax=None):
        if ax is None:
            _, ax = plt.subplots()

        for c, t in self.controls.items():
            c_coords = qubit_coords(c)
            t_coords = qubit_coords(t)
            # draw a circle
            circle = mpl.patches.Circle(
                c_coords,
                radius=CNOT_CIRCLE_RADIUS,
                facecolor="k",
                edgecolor="grey",
                linewidth=2,
                zorder=51,
            )
            ax.add_patch(circle)
            circle = mpl.patches.Circle(
                t_coords,
                radius=CNOT_CIRCLE_RADIUS,
                facecolor="w",
                edgecolor="grey",
                linewidth=2,
                zorder=51,
            )
            ax.add_patch(circle)
            ax.plot(
                [c_coords[0], t_coords[0]],
                [c_coords[1], t_coords[1]],
                color='grey',
                linewidth=5,
                zorder=50,
            )

        if not ax.yaxis_inverted():
            ax.invert_yaxis()
        ax.set_ylabel("Y = Imag = 2nd coord")
        ax.xaxis.tick_top()
        ax.xaxis.set_label_position('top')
        ax.set_xlabel("X = Real = 1st coord")

        ax.autoscale_view()
        ax.grid(True)
        ax.set_aspect('equal')

        return ax
