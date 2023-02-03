from typing import Optional, Literal, Union, Any

from midout import gen
from midout._circuit_case import CircuitCase


def make_square_planar_iswap_code(
        *,
        distance: int,
        basis: Union[Literal['X', 'Z'], str],
        rounds: int,
) -> CircuitCase:
    """
    Args:
        distance: Desired code distance. Determines size of patch to make.
        basis: "X" or "Z" determining the memory experiment basis.
        rounds: Number of times to measure the measure qubits. Must be even.
    """
    assert rounds >= 2 and rounds % 2 == 0
    assert basis in ['X', 'Z']

    layout = ISwapBasicConstrictSurfaceLayoutHelper(distance=distance, basis=basis)
    builder = gen.Builder.for_qubits(layout.used_set)

    # Initial two rounds.
    layout._make_double_cycle(
        builder=builder,
        prev_layer=None,
        save_layer='init')

    # Repeat for desired number of rounds.
    loop_body = builder.fork()
    layout._make_double_cycle(
        builder=loop_body,
        prev_layer='init',
        save_layer='loop')
    builder.circuit += loop_body.circuit * (rounds // 2 - 1)

    # End with data measurement layer.
    if basis == 'X':
        builder.gate("H", {layout.iswap_pure_rev.get(e, e) for e in layout.patch_start.data_set})
        builder.tick()
    builder.measure({layout.iswap_pure_rev.get(e, e) for e in layout.patch_start.data_set},
                    save_layer='loop',
                    tracker_key=lambda e: layout.iswap_pure_ford.get(e, e))
    obs_qs = layout.obs_x_start if basis == 'X' else layout.obs_z_start
    builder.obs_include([gen.AtLayer(q, 'loop') for q in obs_qs], obs_index=0)
    for tile in layout.patch_start.tiles:
        if tile.basis == basis:
            builder.detector([gen.AtLayer(q, 'loop') for q in tile.used_set], pos=tile.measurement_qubit)
    builder.shift_coords(dt=1)
    builder.tick()

    return CircuitCase(
        circuit=builder.circuit,
        patches=[layout.patch_start_draw_nice, layout.plaqs_mid, layout.patch_end_draw_nice],
        expected_interactions=frozenset(['ISWAP']),
    )


class ISwapBasicConstrictSurfaceLayoutHelper:
    def __init__(self,
                 *,
                 distance: int,
                 basis: str):
        """
        Args:
            distance: Desired code distance. Determines size of patch to make.
            basis: Whether this is an X or Z basis memory experiment.
        """
        self.distance = distance
        self.basis = basis

        extended_measure_qubits = {
            x + y * 1j + 0.5 + 0.5j for x in range(-distance * 2, distance * 8) for y in
            range(-distance * 2, distance * 8)
        }

        xUR = +0.5 - 0.5j
        xUL = -0.5 - 0.5j
        xDR = +0.5 + 0.5j
        xDL = -0.5 + 0.5j

        self.plaqs_mid = gen.surface_code_patch(
            distance=distance,
        ).after_coordinate_transform(lambda e: e + 0.5 - 0.5j)
        self.coord_clockwise = lambda e: e + xDL * (+1 if gen.checkerboard_basis(e) == 'X' else -1)
        self.coord_mathwise = lambda e: e + xUR * (+1 if gen.checkerboard_basis(e) == 'X' else -1)
        self.obs_x_start = [self.coord_mathwise(q) for q in self.plaqs_mid.data_set if q.imag == -0.5]
        self.obs_z_start = [self.coord_mathwise(q) for q in self.plaqs_mid.data_set if q.real == 0.5]

        self.patch_start = self.plaqs_mid.after_coordinate_transform(lambda e: self.coord_clockwise(e) if e.real % 1 == 0 else self.coord_mathwise(e))
        self.patch_end = self.plaqs_mid.after_coordinate_transform(lambda e: self.coord_mathwise(e) if e.real % 1 == 0 else self.coord_clockwise(e))

        self.patch_start_draw_nice = self.plaqs_mid.after_coordinate_transform(self.coord_mathwise)
        self.patch_end_draw_nice = self.plaqs_mid.after_coordinate_transform(self.coord_clockwise)
        assert self.patch_start.data_set == self.patch_start_draw_nice.data_set
        assert self.patch_end.data_set == self.patch_end_draw_nice.data_set

        self.lay0_normal_pairs = set()
        for m in self.patch_start.measure_set:
            d = m + xUL
            if d in self.patch_start.data_set:
                self.lay0_normal_pairs.add((d, m))

        self.lay3_normal_pairs = set()
        for m in self.patch_end.measure_set:
            d = m + xDR
            if d in self.patch_end.data_set:
                self.lay3_normal_pairs.add((d, m))

        self.lay1_dropped_boundary_remap = {}
        self.lay1_normal_pairs = set()
        self.lay1_half_active_pairs = set()
        self.lay1_partners_to_dropped_boundary_qubits = set()
        self.lay1_dropped_boundary_qubits = set()
        for m in extended_measure_qubits:
            travel_direction = xUR if gen.checkerboard_basis(m) == 'X' else xDL
            d = m + travel_direction
            m_active = m in self.patch_start.measure_set
            d_active = d in self.patch_start.data_set

            if m_active and d_active:
                self.lay1_normal_pairs.add((d, m))
            elif m_active:
                self.lay1_partners_to_dropped_boundary_qubits.add(d)
                self.lay1_dropped_boundary_qubits.add(m)
                self.lay1_dropped_boundary_remap[d] = m
            elif d_active:
                self.lay1_half_active_pairs.add((m, d))

        self.lay2_dropped_boundary_remap = {}
        self.lay2_normal_pairs = set()
        self.lay2_half_active_pairs = set()
        self.lay2_dropped_boundary_qubits = set()
        self.lay2_partners_to_dropped_boundary_qubits = set()
        for m in extended_measure_qubits:
            travel_direction = xUR if gen.checkerboard_basis(m) == 'X' else xDL
            m += travel_direction
            d = m + travel_direction
            m_active = m in self.plaqs_mid.measure_set
            d_active = d in self.plaqs_mid.data_set

            if m_active and d_active:
                self.lay2_normal_pairs.add((d, m))
            elif m_active:
                self.lay2_dropped_boundary_qubits.add(d)
                self.lay2_partners_to_dropped_boundary_qubits.add(m)
                self.lay2_dropped_boundary_remap[m] = d
            elif d_active:
                self.lay2_half_active_pairs.add((m, d))

        self.dropped_boundary_qubits = self.lay1_dropped_boundary_qubits | self.lay2_dropped_boundary_qubits

        self.used_set = (self.plaqs_mid.used_set | self.patch_end.used_set | self.patch_start.used_set) - self.dropped_boundary_qubits
        self.data_set_start_or_end = (self.patch_start.data_set | self.patch_end.data_set) - self.dropped_boundary_qubits
        self.measure_set_start_or_end = (self.patch_start.measure_set | self.patch_end.measure_set) - self.dropped_boundary_qubits

        self.measure_set_end = self.used_set - self.patch_end.data_set
        self.measure_set_start = self.used_set - self.patch_start.data_set

        self.iswap_pure_rev = dict(self.lay1_dropped_boundary_remap.items())
        self.iswap_pure_ford = {}
        for k, v in self.lay0_normal_pairs:
            self.iswap_pure_ford[v] = k
            self.iswap_pure_rev[k] = v

        self.iswap_pure_rev_2 = dict(self.lay2_dropped_boundary_remap.items())
        self.iswap_pure_ford_2 = {}
        for k, v in self.lay3_normal_pairs:
            self.iswap_pure_ford_2[v] = k
            self.iswap_pure_rev_2[k] = v

    def _make_double_cycle(
            self,
            *,
            builder: gen.Builder,
            prev_layer: Optional[Any],
            save_layer: Any) -> None:

        mid_layer = ('mid', save_layer)
        self._make_forward_cycle(builder=builder, prev_layer=prev_layer, save_layer=mid_layer)
        self._make_reverse_cycle(builder=builder, prev_layer=mid_layer, save_layer=save_layer)

    def _make_forward_cycle(
            self,
            *,
            builder: gen.Builder,
            prev_layer: Optional[Any],
            save_layer: Any) -> None:
        if prev_layer is None:
            builder.gate("R", self.used_set)
        else:
            builder.gate("R", self.used_set - {self.iswap_pure_rev.get(e, e) for e in self.patch_start.data_set})
        builder.tick()

        with builder.plan_rotations() as planner:
            if prev_layer is None and self.basis == 'X':
                planner.gate("H", [self.iswap_pure_rev.get(d, d) for d in self.patch_start.data_set])
            planner.gate("H", [self.iswap_pure_rev.get(d, d) for d in self.patch_start.data_set if gen.checkerboard_basis(d) == 'Z'])
            planner.gate("H_YZ", [e[0] for e in self.lay0_normal_pairs])

        builder.gate2("ISWAP", self.lay0_normal_pairs)
        builder.tick()

        with builder.plan_rotations() as planner:
            planner.gate("S", [e[0] for e in self.lay0_normal_pairs])
            planner.gate("H", [e[0] for e in self.lay1_normal_pairs])
            planner.gate("H_YZ", self.lay1_partners_to_dropped_boundary_qubits)
            planner.gate("H", [e[1] for e in self.lay1_half_active_pairs])

        builder.gate2("ISWAP", self.lay1_normal_pairs | self.lay1_half_active_pairs)
        builder.tick()
        builder.gate2("ISWAP", self.lay2_normal_pairs | self.lay2_half_active_pairs)
        builder.tick()

        with builder.plan_rotations() as planner:
            planner.gate("H", [e[0] for e in self.lay2_normal_pairs])
            planner.gate("H", [e[1] for e in self.lay2_normal_pairs])
            planner.gate("H", [e[0] for e in self.lay2_half_active_pairs])
            planner.gate("S", [e[1] for e in self.lay2_half_active_pairs])
            planner.gate("S", [e[0] for e in self.lay3_normal_pairs])
            planner.gate("C_ZYX", [e[1] for e in self.lay3_normal_pairs])
        builder.gate2("ISWAP", self.lay3_normal_pairs)
        builder.tick()
        with builder.plan_rotations() as planner:
            planner.gate("H", [e[0] for e in self.lay3_normal_pairs])
            planner.gate("H_YZ", self.lay2_partners_to_dropped_boundary_qubits)
            planner.gate("H", [self.iswap_pure_rev_2.get(d, d) for d in self.patch_end.data_set if gen.checkerboard_basis(d) == 'Z'])

        builder.measure(
            [self.iswap_pure_ford_2.get(e, e) for e in self.measure_set_end],
            save_layer=save_layer,
            tracker_key=lambda e: self.iswap_pure_rev_2.get(e, e),
        )
        if prev_layer is None:
            for m in gen.sorted_complex(self.patch_end.measure_set):
                if gen.checkerboard_basis(m) == self.basis:
                    builder.detector([
                        gen.AtLayer(m, save_layer),
                    ], pos=m)
        else:
            for m in gen.sorted_complex(self.patch_end.measure_set):
                if gen.checkerboard_basis(m) == 'Z':
                    offset = 1 - 1j
                else:
                    offset = -1 + 1j
                builder.detector([
                    gen.AtLayer(m, save_layer),
                    gen.AtLayer(m + offset, prev_layer),
                ], pos=m)
        for m in gen.sorted_complex(self.measure_set_end - self.patch_end.measure_set - self.lay2_dropped_boundary_remap.keys()):
            builder.detector([gen.AtLayer(m, save_layer)], pos=m)
        builder.shift_coords(dt=1)
        builder.tick()

    def _make_reverse_cycle(
            self,
            *,
            builder: gen.Builder,
            prev_layer: Optional[Any],
            save_layer: Any) -> None:
        builder.gate("R", [self.iswap_pure_ford_2.get(q, q) for q in self.used_set - self.patch_end.data_set])
        builder.tick()

        with builder.plan_rotations() as planner:
            planner.gate("H", [self.iswap_pure_rev_2.get(d, d) for d in self.patch_end.data_set if gen.checkerboard_basis(d) == 'Z'])
            planner.gate("H", [e[0] for e in self.lay3_normal_pairs])

        builder.gate2("ISWAP", self.lay3_normal_pairs)
        builder.tick()

        with builder.plan_rotations() as planner:
            planner.gate("S", [e[0] for e in self.lay3_normal_pairs])
            planner.gate("C_XYZ", [e[1] for e in self.lay3_normal_pairs])
            planner.gate("H", self.lay2_partners_to_dropped_boundary_qubits)
            planner.gate("C_ZYX", [e[0] for e in self.lay2_half_active_pairs])
            planner.gate("S", [e[1] for e in self.lay2_half_active_pairs])
            planner.gate("C_ZYX", [e[0] for e in self.lay2_normal_pairs])
            planner.gate("C_ZYX", [e[1] for e in self.lay2_normal_pairs])

        builder.gate2("ISWAP", self.lay2_normal_pairs | self.lay2_half_active_pairs)
        builder.tick()
        builder.gate2("ISWAP", self.lay1_normal_pairs | self.lay1_half_active_pairs)
        builder.tick()

        with builder.plan_rotations() as planner:
            planner.gate("C_XYZ", [e[1] for e in self.lay1_half_active_pairs])
            planner.gate("C_XYZ", [e[0] for e in self.lay1_normal_pairs])
            planner.gate("C_XYZ", [e[1] for e in self.lay1_normal_pairs])
            planner.gate("H", [e[1] for e in self.lay0_normal_pairs])

        builder.gate2("ISWAP", self.lay0_normal_pairs)
        builder.tick()

        with builder.plan_rotations() as planner:
            planner.gate("C_XYZ", [e[0] for e in self.lay0_normal_pairs])
            planner.gate("S", [e[1] for e in self.lay0_normal_pairs])
            planner.gate("H", self.lay1_partners_to_dropped_boundary_qubits)
            planner.gate("H", [self.iswap_pure_rev.get(d, d) for d in self.patch_start.data_set if gen.checkerboard_basis(d) == 'Z'])

        builder.measure([self.iswap_pure_ford.get(e, e) for e in self.measure_set_start], save_layer=save_layer, tracker_key=lambda e: self.iswap_pure_rev.get(e, e))
        for m in gen.sorted_complex(self.patch_start.measure_set):
            if gen.checkerboard_basis(m) == 'X':
                offset = 1 - 1j
            else:
                offset = -1 + 1j
            builder.detector([
                gen.AtLayer(m, save_layer),
                gen.AtLayer(m + offset, prev_layer),
            ], pos=m)
        for m in gen.sorted_complex(self.measure_set_start - self.patch_start.measure_set - self.lay1_dropped_boundary_remap.keys()):
            builder.detector([gen.AtLayer(m, save_layer)], pos=m)
        builder.tick()
