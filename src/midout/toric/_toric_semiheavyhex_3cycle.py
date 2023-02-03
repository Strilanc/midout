import functools
from typing import Union, Literal, Any, Optional, Tuple, List, FrozenSet

from midout import gen
from midout._circuit_case import CircuitCase


class BrickToricLayout:
    def __init__(self, *, distance: int):
        self.distance = distance
        self.m_couplers_A, self.m_couplers_B = self._init_m_couplers()
        self.u_couplers_A, self.u_couplers_B = self._init_u_couplers()

    def _init_m_couplers(self) -> Tuple[FrozenSet[Tuple[complex, complex]], FrozenSet[Tuple[complex, complex]]]:
        m_couplers_A = set()
        m_couplers_B = set()
        for x in range(0, 2*self.distance, 2):
            for y in range(0, 2*self.distance, 2):
                a = x + 1j*y
                b = self.wrap(a + 1 + 1j)
                if (x + y) % 4 == 0:
                    m_couplers_A.add((a, b))
                else:
                    m_couplers_B.add((a, b))
        return frozenset(m_couplers_A), frozenset(m_couplers_B)

    def _init_u_couplers(self) -> Tuple[FrozenSet[Tuple[complex, complex]], FrozenSet[Tuple[complex, complex]]]:
        u_couplers_A = set()
        u_couplers_B = set()
        for x in range(0, 2*self.distance):
            for y in range(0, 2*self.distance):
                a = x + 1j*y
                b = self.wrap(a + 1 - 1j)
                if (x - y) % 4 == 0:
                    u_couplers_A.add((a, b))
                elif (x - y) % 4 == 2:
                    u_couplers_B.add((a, b))
        return frozenset(u_couplers_A), frozenset(u_couplers_B)

    @functools.cached_property
    def data_set(self) -> FrozenSet[complex]:
        return frozenset(
            q
            for pair in self.u_couplers_A | self.u_couplers_B | self.m_couplers_A | self.m_couplers_B
            for q in pair
        )

    @functools.cached_property
    def used_set(self) -> FrozenSet[complex]:
        return self.data_set | self.measure_set

    @functools.cached_property
    def measure_set(self) -> FrozenSet[complex]:
        return frozenset(
            tile.measurement_qubit
            for tile in self.patch_A.tiles
            if len(tile.ordered_data_qubits) == 2
        )

    def _make_pair_tile(self, coupler: Tuple[complex, complex], tile_basis: str) -> gen.Tile:
        a, b = coupler
        return gen.Tile(
            bases=tile_basis,
            measurement_qubit=a + 0.5 + 0.5j,
            ordered_data_qubits=[a, b],
        )

    def _make_hex_tile(self, coupler: Tuple[complex, complex], tile_basis: str) -> List[gen.Tile]:
        d = (1 - 1j) * (+1 if tile_basis == 'X' else -1)
        a, b = coupler
        qs = [
            a,
            a + d,
            a + d * 2,
            b + d * 2,
            b + d,
            b,
        ]
        qs = [self.wrap(q) for q in qs]
        if len(qs) != 6:
            return []
        return [gen.Tile(
            bases=tile_basis,
            measurement_qubit=(a + b) / 2,
            ordered_data_qubits=qs,
        )]

    @functools.cached_property
    def coupler_patch(self) -> gen.Patch:
        return gen.Patch([
            *[gen.Tile(bases='X', measurement_qubit=(a + b) / 2, ordered_data_qubits=[a, b]) for a, b in self.m_couplers_A],
            *[gen.Tile(bases='Z', measurement_qubit=(a + b) / 2, ordered_data_qubits=[a, b]) for a, b in self.m_couplers_B],
            *[gen.Tile(bases='Y', measurement_qubit=(a + b) / 2, ordered_data_qubits=[a, b]) for a, b in self.u_couplers_A],
            *[gen.Tile(bases='XY', measurement_qubit=(a + b) / 2, ordered_data_qubits=[a, b]) for a, b in self.u_couplers_B],
        ])

    @functools.cached_property
    def cnot_layer_1A(self) -> Tuple[Tuple[complex, complex], ...]:
        result = []
        for tile in self.patch_A.tiles:
            if len(tile.ordered_data_qubits) != 2:
                continue
            q = tile.measurement_qubit
            dia = q.real - q.imag
            sign = +1 if dia % 4 == 0 else -1
            order = -1 if tile.basis == 'X' else +1
            result.append((q, self.wrap(q + (0.5 + 0.5j) * sign))[::order])
        for q in self.used_set:
            anti = q.real + q.imag
            dia = q.real - q.imag
            if dia % 4 == 0 and anti % 4 == 0:
                result.append((q, self.wrap(q + 1 - 1j)))
        return tuple(result)

    @functools.cached_property
    def cnot_layer_2A(self) -> Tuple[Tuple[complex, complex], ...]:
        result = []
        for q in self.used_set:
            anti = q.real + q.imag
            dia = q.real - q.imag
            if dia % 4 == 0 and anti % 4 == 0:
                result.append((self.wrap(q - 1 + 1j), q))
            if dia % 4 == 0 and anti % 4 == 2:
                result.append((self.wrap(q - 1 + 1j), q))
        return tuple(result)

    @functools.cached_property
    def cnot_layer_3A(self) -> Tuple[Tuple[complex, complex], ...]:
        result = []
        for q in self.used_set:
            anti = q.real + q.imag
            dia = q.real - q.imag
            if dia % 4 == 0 and anti % 4 == 2:
                result.append((q, self.wrap(q + 1 - 1j)))
        for tile in self.patch_A.tiles:
            if len(tile.ordered_data_qubits) != 2:
                continue
            q = tile.measurement_qubit
            dia = q.real - q.imag
            sign = -1 if dia % 4 == 0 else +1
            order = -1 if tile.basis == 'X' else +1
            result.append((q, self.wrap(q + (0.5 + 0.5j) * sign))[::order])
        return tuple(result)

    @functools.cached_property
    def cnot_layer_1B(self) -> Tuple[Tuple[complex, complex], ...]:
        result = []
        for tile in self.patch_B.tiles:
            if len(tile.ordered_data_qubits) != 2:
                continue
            q = tile.measurement_qubit
            dia = q.real - q.imag
            sign = +1 if dia % 4 == 0 else -1
            order = -1 if tile.basis == 'X' else +1
            result.append((q, self.wrap(q + (0.5 + 0.5j) * sign))[::order])
        for q in self.used_set:
            anti = q.real + q.imag
            dia = q.real - q.imag
            if dia % 4 == 0 and anti % 4 == 0:
                result.append((self.wrap(q - 1 + 1j), q))
        return tuple(result)

    @functools.cached_property
    def cnot_layer_2B(self) -> Tuple[Tuple[complex, complex], ...]:
        result = []
        for q in self.used_set:
            anti = q.real + q.imag
            dia = q.real - q.imag
            if dia % 4 == 0 and anti % 4 == 0:
                result.append((q, self.wrap(q + 1 - 1j)))
            if dia % 4 == 0 and anti % 4 == 2:
                result.append((q, self.wrap(q + 1 - 1j)))
        return tuple(result)

    @functools.cached_property
    def cnot_layer_3B(self) -> Tuple[Tuple[complex, complex], ...]:
        result = []
        for q in self.used_set:
            anti = q.real + q.imag
            dia = q.real - q.imag
            if dia % 4 == 0 and anti % 4 == 2:
                result.append((self.wrap(q - 1 + 1j), q))
        for tile in self.patch_B.tiles:
            if len(tile.ordered_data_qubits) != 2:
                continue
            q = tile.measurement_qubit
            dia = q.real - q.imag
            sign = -1 if dia % 4 == 0 else +1
            order = -1 if tile.basis == 'X' else +1
            result.append((q, self.wrap(q + (0.5 + 0.5j) * sign))[::order])
        return tuple(result)

    @functools.cached_property
    def patch_A(self) -> gen.Patch:
        tiles_A = []
        for c in self.m_couplers_A:
            tiles_A.append(self._make_pair_tile(c, 'X'))
            tiles_A.extend(self._make_hex_tile(c, 'Z'))
        for c in self.m_couplers_B:
            tiles_A.append(self._make_pair_tile(c, 'Z'))
            tiles_A.extend(self._make_hex_tile(c, 'X'))

        return gen.Patch(tiles_A)

    @functools.cached_property
    def patch_B(self) -> gen.Patch:
        tiles_B = []
        for c in self.m_couplers_A:
            tiles_B.append(self._make_pair_tile(c, 'Z'))
            tiles_B.extend(self._make_hex_tile(c, 'X'))
        for c in self.m_couplers_B:
            tiles_B.append(self._make_pair_tile(c, 'X'))
            tiles_B.extend(self._make_hex_tile(c, 'Z'))
        return gen.Patch(tiles_B)

    def wrap(self, q: complex) -> complex:
        r = q.real % (2 * self.distance)
        i = q.imag % (2 * self.distance)
        return r + 1j*i

    def patches(self) -> List[gen.Patch]:
        return [self.patch_A, self.patch_B, self.coupler_patch]


def make_semiheavyhex_toric_cx_3cycle_code(
        *,
        distance: int,
        basis: Union[Literal['X', 'Z'], str],
        rounds: int) -> CircuitCase:
    assert rounds >= 2
    assert rounds % 2 == 0

    layout = BrickToricLayout(distance=distance)

    def do_round(*, out: gen.Builder, save_layer: Any, cmp_layer: Optional[Any], order: bool):
        if not order:
            patch = layout.patch_B
            cx1 = layout.cnot_layer_1A
            cx2 = layout.cnot_layer_2A
            cx3 = layout.cnot_layer_3A
        else:
            patch = layout.patch_A
            cx1 = layout.cnot_layer_1B
            cx2 = layout.cnot_layer_2B
            cx3 = layout.cnot_layer_3B
        tiles2 = [tile for tile in patch.tiles if len(tile.ordered_data_qubits) == 2]

        out.gate("R", patch.measure_set)
        out.tick()
        out.gate('H', [tile.measurement_qubit for tile in tiles2 if tile.basis == 'X'])
        out.tick()

        out.gate2("CX", cx1)
        out.tick()
        out.gate2("CX", cx2)
        out.tick()
        out.gate2("CX", cx3)
        out.tick()

        out.gate('H', [tile.measurement_qubit for tile in tiles2 if tile.basis == 'X'])
        out.tick()
        out.measure(patch.measure_set, save_layer=save_layer)
        for tile in patch.tiles:
            if len(tile.ordered_data_qubits) == 2:
                m = tile.measurement_qubit
                if cmp_layer is not None:
                    out.detector([gen.AtLayer(m, save_layer), gen.AtLayer(m, cmp_layer)], pos=m)
                elif tile.basis == basis:
                    out.detector([gen.AtLayer(m, save_layer)], pos=m)
        out.shift_coords(dt=1)
        out.tick()

    builder = gen.Builder.for_qubits(layout.used_set)
    builder.gate("R", layout.used_set)
    builder.tick()
    if basis == 'X':
        builder.gate('H', layout.data_set)
        builder.tick()
    do_round(out=builder, save_layer="A_init", cmp_layer=None, order=False)
    do_round(out=builder, save_layer="B_init", cmp_layer=None, order=True)
    loop = builder.fork()
    do_round(out=loop, save_layer="A_loop", cmp_layer="A_init", order=False)
    do_round(out=loop, save_layer="B_loop", cmp_layer="B_init", order=True)
    builder.circuit += loop.circuit * (rounds // 2 - 1)

    if basis == 'X':
        builder.gate('H', layout.data_set)
        builder.tick()
    builder.measure(layout.data_set, save_layer='end')
    for tile in layout.patch_A.tiles:
        if tile.basis == basis:
            offsets = [-0.5 - 0.5j, -1.5 + 0.5j, +0.5 + 0.5j, -2.5 + 1.5j]
            sign = +1 if basis == 'Z' else -1
            from_data = [
                gen.AtLayer(layout.wrap(tile.measurement_qubit + d * sign), 'end')
                for d in offsets
            ]
            from_measure = [gen.AtLayer(
                tile.measurement_qubit,
                'B_loop' if len(tile.ordered_data_qubits) == 2 else 'A_loop')
            ]
            builder.detector(from_data + from_measure, pos=tile.measurement_qubit)

    x_qs = {q for q in layout.data_set if q.imag in [1]}
    z_qs = {q for q in layout.data_set if q.real in [0, 1]}
    assert len(x_qs & z_qs) % 2 == 1
    obs_qs = x_qs if basis == 'X' else z_qs
    builder.obs_include([gen.AtLayer(q, 'end') for q in obs_qs], obs_index=0)

    return CircuitCase(
        circuit=builder.circuit,
        patches=layout.patches(),
        expected_interactions=frozenset(['CX']),
    )
