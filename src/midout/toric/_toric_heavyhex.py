from typing import Union, Literal, Tuple, Any, Optional, Set

from midout import gen
from midout._circuit_case import CircuitCase


def make_heavyhex_toric_cx_code(
        *,
        distance: int,
        basis: Union[Literal['X', 'Z'], str],
        rounds: int) -> CircuitCase:
    assert rounds >= 2
    assert rounds % 2 == 0

    def wrap(q: complex) -> complex:
        r = q.real % (distance * 2)
        i = q.imag % (distance * 2)
        return r + 1j*i

    data_set = {
        x + y*1j
        for x in range(2 * distance)
        for y in range(2 * distance)
        if (x + y) % 2 == 0
    }
    route_measure_set = {wrap(q + 0.5 - 0.5j) for q in data_set}
    route_measure_set_A = {q for q in route_measure_set if (q.real - q.imag) % 4 == 1}
    route_measure_set_B = route_measure_set - route_measure_set_A
    cross_measure_set = {
        wrap(q + 0.5 + 0.5j) for q in data_set
        if q.real % 2 == 0
    }
    route_measure_set_flipped = {q for q in route_measure_set if (q.real + q.imag) % 4 == 0}
    cross_measure_set_half_1 = {q for q in cross_measure_set if (q.real - q.imag) % 4 == 0}
    measure_set = route_measure_set | cross_measure_set
    used_set = measure_set | data_set

    def cnot_measure_step(sign: int, inv: int) -> Set[Tuple[complex, complex]]:
        return {(wrap(q + (0.5 + 0.5j) * sign), q)[::inv if q in cross_measure_set_half_1 else -inv] for q in cross_measure_set}
    def cnot_across_step(across_set: Set[complex], sign: int) -> Set[Tuple[complex, complex]]:
        return {(wrap(q + (0.5 - 0.5j) * sign), q)[::sign] for q in across_set}
    def cnot_across_step_2(across_set: Set[complex], sign: int) -> Set[Tuple[complex, complex]]:
        return cnot_across_step(across_set - route_measure_set_flipped, sign) | cnot_across_step(across_set & route_measure_set_flipped, -sign)

    pairsA1 = cnot_across_step_2(route_measure_set_A, +1)
    pairsA2 = cnot_across_step_2(route_measure_set_A, -1)
    pairsA3 = cnot_measure_step(+1, +1)
    pairsA4 = cnot_measure_step(-1, +1)
    pairsB1 = cnot_across_step_2(route_measure_set_B, +1)
    pairsB2 = cnot_across_step_2(route_measure_set_B, -1)
    pairsB3 = cnot_measure_step(-1, -1)
    pairsB4 = cnot_measure_step(+1, -1)
    hsA = (cross_measure_set - cross_measure_set_half_1) ^ route_measure_set_flipped
    hsB = (cross_measure_set & cross_measure_set_half_1) ^ route_measure_set_flipped
    if basis == 'X':
        obs_qs = {q for q in used_set if q.imag == 0 and q.real % 2 == 0}
    else:
        obs_qs = {q for q in used_set if q.real == 0 and q.imag % 2 == 0}

    tiles = []
    for q in gen.sorted_complex(data_set):
        tile_basis = "X" if q.real % 2 == 0 else "Z"
        tiles.append(gen.Tile(
            bases=tile_basis,
            measurement_qubit=wrap(q + (1.5j - 0.5 if tile_basis == 'Z' else 0.5 + 0.5j)),
            ordered_data_qubits=[q, wrap(q + 1 + 1j), wrap(q + 2j), wrap(q - 1 + 1j)],
        ))

    builder = gen.Builder.for_qubits(used_set)

    def do_round_a(*, out: gen.Builder, save_layer: Any, cmp_layer: Optional[Any]):
        if cmp_layer is not None:
            out.gate2("CX", pairsB1)
        out.gate2("CX", pairsA1)
        out.tick()
        out.gate2("CX", pairsA2)
        out.tick()
        out.gate2("CX", pairsA1 | pairsA3)
        out.tick()
        out.gate2("CX", pairsA4)
        out.tick()

        out.gate("H", hsA)
        out.tick()
        out.measure(measure_set, save_layer=save_layer)
        for q in gen.sorted_complex(measure_set):
            q_basis = 'Z' if q in cross_measure_set_half_1 else 'X'
            if q in route_measure_set:
                # Flag qubit.
                out.detector([gen.AtLayer(q, save_layer)], pos=q)
            elif cmp_layer is not None:
                # Compare to previous measurement.
                out.detector([gen.AtLayer(q, save_layer), gen.AtLayer(q, cmp_layer)], pos=q)
            elif q_basis == basis:
                # Compare to initialization.
                out.detector([gen.AtLayer(q, save_layer)], pos=q)
        out.shift_coords(dt=1)
        out.tick()

        out.gate("R", measure_set)
        out.tick()
        out.gate("H", hsB)
        out.tick()

        out.gate2("CX", pairsA1)
        out.tick()
        out.gate2("CX", pairsA2)
        out.tick()
        out.gate2("CX", pairsA1)

    def do_round_b(*, out: gen.Builder, save_layer: Any, cmp_layer: Any):
        out.gate2("CX", pairsB1)
        out.tick()
        out.gate2("CX", pairsB2)
        out.tick()
        out.gate2("CX", pairsB1 | pairsB3)
        out.tick()
        out.gate2("CX", pairsB4)
        out.tick()

        out.gate("H", hsB)
        out.tick()
        out.measure(measure_set, save_layer=save_layer)
        for q in gen.sorted_complex(measure_set):
            q_basis = 'X' if q in cross_measure_set_half_1 else 'Z'
            if q in route_measure_set:
                # Flag qubit.
                out.detector([gen.AtLayer(q, save_layer)], pos=q)
            elif cmp_layer is not None:
                # Compare to previous measurement.
                out.detector([gen.AtLayer(q, save_layer), gen.AtLayer(q, cmp_layer)], pos=q)
            elif q_basis == basis:
                # Compare to initialization.
                out.detector([gen.AtLayer(q, save_layer)], pos=q)
        out.shift_coords(dt=1)
        out.tick()

    builder.gate("R", used_set)
    builder.tick()
    builder.gate('H', hsA)
    if basis == 'X':
        builder.gate('H', data_set)
    builder.tick()
    do_round_a(out=builder, save_layer="A_init", cmp_layer=None)
    do_round_b(out=builder, save_layer="B_init", cmp_layer=None)

    loop = builder.fork()
    loop.gate("R", measure_set)
    loop.tick()
    loop.gate("H", hsA)
    loop.tick()
    loop.gate2("CX", pairsB1)
    loop.tick()
    loop.gate2("CX", pairsB2)
    loop.tick()
    do_round_a(out=loop, save_layer="A_loop", cmp_layer="A_init")
    do_round_b(out=loop, save_layer="B_loop", cmp_layer="B_init")
    builder.circuit += loop.circuit * (rounds // 2 - 1)
    builder.gate("R", measure_set)
    builder.tick()
    builder.gate("H", hsA)
    builder.tick()
    builder.gate2("CX", pairsB1)
    builder.tick()
    builder.gate2("CX", pairsB2)
    builder.tick()

    builder.gate2("CX", pairsB1)
    builder.tick()
    if basis == 'X':
        builder.gate('H', data_set)
    builder.gate('H', hsA)
    builder.tick()
    builder.measure(used_set, save_layer='end')
    for q in gen.sorted_complex(measure_set):
        builder.detector([gen.AtLayer(q, 'end')], pos=q)
    for tile in tiles:
        if tile.basis == basis:
            layer = 'B_loop' if (tile.measurement_qubit in cross_measure_set_half_1) ^ (basis == 'Z') else 'A_loop'
            builder.detector([gen.AtLayer(q, 'end') for q in tile.ordered_data_qubits] + [gen.AtLayer(tile.measurement_qubit, layer)], pos=tile.measurement_qubit)

    builder.obs_include([gen.AtLayer(q, 'end') for q in obs_qs], obs_index=0)

    patch = gen.Patch([tile for tile in tiles if all(abs(tile.ordered_data_qubits[0] - q) < 5 for q in tile.ordered_data_qubits)])
    return CircuitCase(
        circuit=builder.circuit,
        patches=[patch],
        expected_interactions=frozenset(['CX']),
    )
