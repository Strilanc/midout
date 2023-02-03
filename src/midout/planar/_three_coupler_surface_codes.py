import functools
from typing import Any, FrozenSet, Tuple, Union, Literal

import stim

from midout import gen
from midout._circuit_case import CircuitCase
from midout.gen._layer_translate import to_z_basis_interaction_circuit


def make_hex_planar_surface_code(
        *,
        basis: str,
        distance: int,
        gate: Union[str, Literal[
            'CX',
            'CXSWAP',
            'CZ',
            'ISWAP',
        ]],
        wiggle: bool,
        rounds: int,
) -> CircuitCase:
    """
    Args:
        basis: The type of memory experiment. The basis of the time boundaries.
        distance: The desired code distance of the circuit.
        gate: The desired kind of circuit.
        rounds: The desired number of rounds between logical initialization and
            logical measurement.
        wiggle: Alternate between measurement sets so all qubits get reset.
    """

    kak_gate = gate
    if gate == 'CZ':
        kak_gate = 'CX'
    if gate == 'ISWAP':
        kak_gate = 'CXSWAP'
    if gate == 'CZ_MZZ':
        kak_gate = 'CX_MXX_MZZ'
    layout = ThreeCouplerLayoutHelper(
        basis=basis,
        distance=distance,
        gate=kak_gate,
        rounds=rounds,
        wiggle=wiggle,
    )
    circuit = layout.make_ideal_circuit()
    patch = layout.patch.after_coordinate_transform(lambda e: e / 2)
    circuit = gen.stim_circuit_with_transformed_coords(circuit, lambda e: e / 2)
    if gate == 'CZ' or gate == 'ISWAP' or gate == 'CZ_MZZ':
        circuit = to_z_basis_interaction_circuit(circuit)
    expected_interactions = frozenset([gate])
    if gate == 'CXSWAP':
        expected_interactions = frozenset(['CXSWAP'])
    if gate == 'CZ_MZZ':
        expected_interactions = frozenset(['CZ', 'MZZ'])
    if gate == 'CX_MXX_MZZ':
        expected_interactions = frozenset(['CX', 'MZZ', 'MXX', 'MX'])
    return CircuitCase(
        circuit=circuit,
        patches=[patch],
        expected_interactions=expected_interactions,
    )


class ThreeCouplerLayoutHelper:
    """Helper class with useful properties for using three couplers."""
    def __init__(
            self,
            *,
            basis: str,
            distance: int,
            gate: str,
            rounds: int,
            wiggle: bool,
    ):
        assert rounds >= 2
        self.basis = basis
        self.distance = distance
        self.gate = gate
        self.rounds = rounds
        self.wiggle_a = False
        self.wiggle_b = wiggle

    @property
    def use_more_spikes(self) -> bool:
        return self.gate == 'CXSWAP'

    @functools.cached_property
    def boundary_spike_qubits_x(self) -> FrozenSet[complex]:
        result = []
        for x in range(self.distance - 1):
            result.append(self.distance * 2j + x * 2 - 1j + 1)
            if self.use_more_spikes:
                result.append(x * 2 - 1j + 1)
        return frozenset(result)

    @functools.cached_property
    def boundary_spike_qubits_z(self) -> FrozenSet[complex]:
        return frozenset(q.imag + q.real * 1j for q in self.boundary_spike_qubits_x)

    @functools.cached_property
    def boundary_spike_point_tiles(self) -> gen.Patch:
        return gen.Patch([
            *[
                gen.Tile(
                    bases="X",
                    measurement_qubit=c,
                    ordered_data_qubits=[c],
                ) for c in self.boundary_spike_qubits_x
            ],
            *[
                gen.Tile(
                    bases="Z",
                    measurement_qubit=c,
                    ordered_data_qubits=[c],
                ) for c in self.boundary_spike_qubits_z
            ],
        ])

    def is_tile_center(self, q: complex) -> bool:
        if q.real % 2 == q.imag % 2:
            return False
        if min(q.real, q.imag) < 0:
            return False
        return True

    def is_data_qubit(self, q: complex) -> bool:
        return min(q.real, q.imag) >= 0 or q in self.boundary_spike_qubits_x | self.boundary_spike_qubits_z

    @functools.cached_property
    def patch(self) -> gen.Patch:
        adjacencies_x = [1j**d for d in range(4)][::-1]
        adjacencies_z = [1j**-d for d in range(3, 7)][::-1]

        tiles = [*self.boundary_spike_point_tiles.tiles]

        for x in range(-3, self.distance * 2 - 1):
            for y in range(-3, self.distance * 2 - 1):
                center = x + y * 1j
                if not self.is_tile_center(center):
                    continue
                basis = "X" if center.real % 2 == 1 else 'Z'
                if (center.real - center.imag) % 4 == 1:
                    wiggle = self.wiggle_a
                else:
                    wiggle = self.wiggle_b
                if basis == 'X':
                    adjacencies = adjacencies_x
                else:
                    adjacencies = adjacencies_z
                if wiggle:
                    adjacencies = adjacencies[::-1]
                data_qubits = [
                    center + d
                    for d in adjacencies
                    if self.is_data_qubit(center + d)
                ]
                if len(data_qubits) < 2:
                    continue
                tiles.append(gen.Tile(
                    bases=basis,
                    measurement_qubit=center,
                    ordered_data_qubits=data_qubits,
                ))

        return gen.Patch(tiles)

    @functools.cached_property
    def even_sums(self) -> FrozenSet[complex]:
        return frozenset(
            q
            for q in self.patch.data_set
            if (q.real + q.imag) % 4 == 0
        )

    @functools.cached_property
    def odd_sums(self) -> FrozenSet[complex]:
        return frozenset(
            q
            for q in self.patch.data_set
            if (q.real + q.imag) % 4 == 2
        )

    @functools.cached_property
    def even_difs(self) -> FrozenSet[complex]:
        return frozenset(
            q
            for q in self.patch.data_set
            if (q.real - q.imag) % 4 == 0
        )

    @functools.cached_property
    def odd_difs(self) -> FrozenSet[complex]:
        return frozenset(
            q
            for q in self.patch.data_set
            if (q.real - q.imag) % 4 == 2
        )

    @functools.cached_property
    def pairsA1(self) -> FrozenSet[Tuple[complex, complex]]:
        return frozenset(
            (q, p)
            for q in self.patch.data_set
            if q in self.even_difs
            if (p := q + 1 - 1j) in self.patch.data_set
        )

    @functools.cached_property
    def pairsA2(self) -> FrozenSet[Tuple[complex, complex]]:
        return frozenset(
            (q, p)[::-1 if (q in self.even_difs) ^ self.wiggle_a else +1]
            for q in self.patch.data_set
            if q.real % 2 != 1
            if (p := q + 1 + 1j) in self.patch.data_set
            if (p in self.boundary_spike_qubits_z) <= (q in self.odd_sums)
            if (p in self.boundary_spike_qubits_x) <= (q in self.even_sums)
        )

    @property
    def start_backwards(self) -> bool:
        return self.rounds % 2 == 1

    @functools.cached_property
    def measure_qubits_A(self) -> FrozenSet[complex]:
        result = frozenset(
            q
            for q in self.patch.data_set
            if q.real % 2 == 1
            and (q.imag != -1 or q.real % 4 == 3)
            and (q.real != -1 or q.imag % 4 == 1)
        )
        if self.wiggle_a:
            m = {}
            for a, b in self.pairsA2:
                m[a] = b
                m[b] = a
            result = frozenset(m.get(q, q) for q in result)

        return result

    @functools.cached_property
    def measure_qubits_B(self) -> FrozenSet[complex]:
        result = frozenset(
            q
            for q in self.patch.data_set
            if q.real % 2 == 1
            and (q.imag != -1 or q.real % 4 == 1)
            and (q.real != -1 or q.imag % 4 == 3)
        )
        if self.wiggle_b:
            m = {}
            for a, b in self.pairsB2:
                m[a] = b
                m[b] = a
            result = frozenset(m.get(q, q) for q in result)
        return result

    @functools.cached_property
    def hadamard_qubits_A(self) -> FrozenSet[complex]:
        return frozenset(
            q
            for q in self.patch.data_set
            if q in self.boundary_spike_qubits_x or (
                    q in self.measure_qubits_A
                    and q not in self.odd_difs
                    and (q not in self.even_sums) ^ self.wiggle_a
                    and q not in self.boundary_spike_qubits_z
            )
        )

    @functools.cached_property
    def pairsB1(self) -> FrozenSet[Tuple[complex, complex]]:
        return frozenset(
            (q, p)
            for q in self.patch.data_set
            if q in self.odd_difs
            if (p := q + 1 - 1j) in self.patch.data_set
        )

    @functools.cached_property
    def pairsB2(self) -> FrozenSet[Tuple[complex, complex]]:
        return frozenset(
            (q, p)[::-1 if (q in self.odd_difs) ^ self.wiggle_b else +1]
            for q in self.patch.data_set
            if q.real % 2 != 1
            if (p := q + 1 + 1j) in self.patch.data_set
            if (p in self.boundary_spike_qubits_z) <= (q in self.even_sums)
            if (p in self.boundary_spike_qubits_x) <= (q in self.odd_sums)
        )

    @functools.cached_property
    def hadamard_qubits_B(self) -> FrozenSet[complex]:
        return frozenset(
            q
            for q in self.patch.data_set
            if q in self.boundary_spike_qubits_x or (
                    q in self.measure_qubits_B
                    and q in self.odd_difs
                    and (q in self.even_sums) ^ self.wiggle_b
                    and q not in self.boundary_spike_qubits_z
            )
        )

    def make_circuit_layer_type_A_measure(self, *, out: gen.Builder, cmp_layer: Any, save_layer: Any) -> None:
        # Demolition measure the temporarily-single-qubit stabilizers.
        if self.gate == 'CX_MXX_MZZ':
            partner = {}
            for a, b in self.pairsA2:
                partner[a] = b
                partner[b] = a
            for q in self.measure_qubits_A:
                basis = 'X' if q in self.hadamard_qubits_A else 'Z'
                other = partner.get(q)
                targets = [q] if other is None else [q, other]
                out.measure_pauli_product(b2qs={basis: targets}, key=gen.AtLayer(q, save_layer))
        elif self.gate in ['CX', 'CXSWAP']:
            out.measure(self.measure_qubits_A & self.hadamard_qubits_A, save_layer=save_layer, basis='X')
            out.measure(self.measure_qubits_A - self.hadamard_qubits_A, save_layer=save_layer, basis='Z')
        else:
            raise NotImplementedError(f'{self.gate=}')

        if cmp_layer is None:
            for m in gen.sorted_complex(self.measure_qubits_A):
                if m in self.boundary_spike_qubits_z:
                    keep = self.basis == 'Z'
                elif m in self.boundary_spike_qubits_x:
                    keep = self.basis == 'X'
                elif (m in self.odd_sums) ^ self.wiggle_a:
                    keep = self.basis == 'X'
                else:
                    keep = self.basis == 'Z'
                if keep:
                    out.detector([gen.AtLayer(m, save_layer)], pos=m)
            out.shift_coords(dt=1)
        else:
            for m in gen.sorted_complex(self.measure_qubits_A):
                out.detector([gen.AtLayer(m, cmp_layer), gen.AtLayer(m, save_layer)], pos=m)
            out.shift_coords(dt=1)
        out.tick()
        if self.gate in ['CX', 'CXSWAP']:
            out.gate("R", self.measure_qubits_A - self.hadamard_qubits_A)
            out.gate("RX", self.measure_qubits_A & self.hadamard_qubits_A)
            for m in gen.sorted_complex(self.measure_qubits_A):
                out.classical_paulis(control_keys=[gen.AtLayer(m, save_layer)], targets=[m], basis='Z' if m in self.hadamard_qubits_A else 'X')
            out.tick()

    def make_circuit_layer_type_AtoB(self, *, out: gen.Builder) -> None:
        for k, p in self.transition_a_to_b:
            if k in ['CX', 'CXSWAP', 'SWAPCX']:
                out.gate2(k, p)
            else:
                out.gate(k, p)
            out.tick()

    def make_circuit_layer_type_BtoA(self, *, out: gen.Builder) -> None:
        for k, p in self.transition_a_to_b[::-1]:
            if k in ['CX', 'CXSWAP', 'SWAPCX']:
                if k == 'CXSWAP':
                    k = 'SWAPCX'
                elif k == 'SWAPCX':
                    k = 'CXSWAP'
                out.gate2(k, p)
            else:
                out.gate(k, p)
            out.tick()

    @functools.cached_property
    def transition_a_to_b(self) -> Tuple[Tuple[str, Any], ...]:
        if self.gate == 'CX':
            return (
                ("CX", self.pairsA2),
                ("CX", self.pairsA1),
                ("CX", self.pairsB1),
                ("CX", self.pairsB2),
            )
        elif self.gate == 'CXSWAP':
            return (
                ("CXSWAP", self.pairsA2),
                ("CXSWAP", self.pairsA1),
                ("SWAPCX", self.pairsB1),
                ("SWAPCX", self.pairsB2),
            )
        elif self.gate == 'CX_MXX_MZZ':
            return (
                ("CX", self.pairsA1),
                ("CX", self.pairsB1),
            )
        else:
            raise NotImplementedError(f"{self.gate=}")

    def make_circuit_layer_type_B_end(self, *, out: gen.Builder) -> None:
        # Undo interactions to restore the single qubit observables back to 4-body stabilizers.
        if self.gate == 'CX':
            out.gate2("CX", self.pairsB2)
            out.tick()
            out.gate2("CX", self.pairsB1)
            out.tick()
        elif self.gate == 'CXSWAP':
            out.gate2("CXSWAP", self.pairsB2)
            out.tick()
            out.gate2("CXSWAP", self.pairsB1)
            out.tick()
        elif self.gate == 'CX_MXX_MZZ':
            out.gate2("CX", self.pairsB1)
            out.tick()
        else:
            raise NotImplementedError(f"{self.gate=}")

    def make_circuit_layer_type_B_measure(self, *, out: gen.Builder, cmp_layer: Any, save_layer: Any) -> None:
        if self.gate in ['CX', 'CXSWAP']:
            out.measure(self.measure_qubits_B & self.hadamard_qubits_B, save_layer=save_layer, basis='X')
            out.measure(self.measure_qubits_B - self.hadamard_qubits_B, save_layer=save_layer, basis='Z')
        elif self.gate == 'CX_MXX_MZZ':
            partner = {}
            for a, b in self.pairsB2:
                partner[a] = b
                partner[b] = a
            for q in self.measure_qubits_B:
                basis = 'X' if q in self.hadamard_qubits_B else 'Z'
                other = partner.get(q)
                targets = [q] if other is None else [q, other]
                out.measure_pauli_product(b2qs={basis: targets}, key=gen.AtLayer(q, save_layer))
        else:
            raise NotImplementedError(f'{self.gate=}')
        if cmp_layer is None:
            for m in gen.sorted_complex(self.measure_qubits_B):
                if m in self.boundary_spike_qubits_z:
                    keep = self.basis == 'Z'
                elif m in self.boundary_spike_qubits_x:
                    keep = self.basis == 'X'
                elif (m in self.odd_sums) ^ self.wiggle_b:
                    keep = self.basis == 'Z'
                else:
                    keep = self.basis == 'X'
                if keep:
                    out.detector([gen.AtLayer(m, save_layer)], pos=m)
            out.shift_coords(dt=1)
        else:
            for m in gen.sorted_complex(self.measure_qubits_B):
                out.detector([gen.AtLayer(m, cmp_layer), gen.AtLayer(m, save_layer)], pos=m)
            out.shift_coords(dt=1)
        out.tick()
        if self.gate in ['CX', 'CXSWAP']:
            out.gate("R", self.measure_qubits_B - self.hadamard_qubits_B)
            out.gate("RX", self.measure_qubits_B & self.hadamard_qubits_B)
            for m in gen.sorted_complex(self.measure_qubits_B):
                out.classical_paulis(control_keys=[gen.AtLayer(m, save_layer)], targets=[m], basis='Z' if m in self.hadamard_qubits_B else 'X')
            out.tick()

    def make_ideal_circuit(self) -> stim.Circuit:
        builder = gen.Builder.for_qubits(self.patch.data_set)

        # Prepare individual qubits into the basis needed for the very first interactions.
        x_init_qubits = set()
        if self.basis == 'X':
            x_init_qubits ^= self.patch.data_set
        builder.gate("RX", self.patch.data_set & x_init_qubits)
        builder.gate("R", self.patch.data_set - x_init_qubits)
        builder.tick()

        # Perform initial rounds to get the system into a steady state for the loop.
        if self.start_backwards:
            cmp_b_first = "pre_initB"
            self.make_circuit_layer_type_B_measure(out=builder, cmp_layer=None, save_layer="pre_initB")
            self.make_circuit_layer_type_BtoA(out=builder)
        else:
            cmp_b_first = None
        self.make_circuit_layer_type_A_measure(out=builder, cmp_layer=None, save_layer="initA")
        self.make_circuit_layer_type_AtoB(out=builder)

        self.make_circuit_layer_type_B_measure(out=builder, cmp_layer=cmp_b_first, save_layer="initB")

        # Loop until the requested number of rounds have been performed.
        loop = builder.fork()
        self.make_circuit_layer_type_BtoA(out=loop)
        self.make_circuit_layer_type_A_measure(out=loop, cmp_layer="initA", save_layer="loopA")
        self.make_circuit_layer_type_AtoB(out=loop)
        self.make_circuit_layer_type_B_measure(out=loop, cmp_layer="initB", save_layer="loopB")
        builder.circuit += loop.circuit * (self.rounds // 2 - 1)

        # Rotate data qubits from the basis used for the last interaction to their measurement basis and measure them.
        self.make_circuit_layer_type_B_end(out=builder)
        builder.measure(self.patch.data_set, save_layer="end", basis=self.basis)

        # Compare data qubit measurements to known stabilizer values from previous measurements.
        for tile in self.patch.tiles:
            if tile.basis != self.basis:
                continue
            a_not_b = (tile.measurement_qubit.real - tile.measurement_qubit.imag) % 4 == 1
            offset = 1 if self.basis == 'Z' else 1j
            if (a_not_b and self.wiggle_a) or (not a_not_b and self.wiggle_b):
                offset = -1j if self.basis == 'Z' else -1
            if len(tile.data_set) == 1:
                d, = tile.data_set
                if d in self.measure_qubits_A - self.measure_qubits_B:
                    a_not_b = True
                    offset = 0
                elif d in self.measure_qubits_B - self.measure_qubits_A:
                    a_not_b = False
                    offset = 0
                elif d in self.boundary_spike_qubits_z:
                    a_not_b = d in self.odd_sums
                    offset = 0
                elif d in self.boundary_spike_qubits_x:
                    a_not_b = d in self.even_sums
                    offset = 0
            reference = gen.AtLayer(tile.measurement_qubit + offset, "loopA" if a_not_b else "loopB")
            builder.detector([
                                 gen.AtLayer(d, "end") for d in tile.ordered_data_qubits
                             ] + [reference], pos=tile.measurement_qubit)

        # Extract observables from data measurements.
        obs_xs = {q for q in self.patch.data_set if q.real == 0}
        obs_zs = {q for q in self.patch.data_set if q.imag == 0}
        obs_qs = obs_xs if self.basis == 'X' else obs_zs
        builder.obs_include([gen.AtLayer(q, "end") for q in obs_qs], obs_index=0)

        return builder.circuit.with_inlined_feedback()
