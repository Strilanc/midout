from typing import Union, Literal

from midout import gen
from midout._circuit_case import CircuitCase


def make_square_toric_cx_code(
        *,
        distance: int,
        basis: Union[Literal['X', 'Z'], str],
        rounds: int) -> CircuitCase:
    assert rounds >= 1

    patch = toric_patch(distance=distance)
    builder = gen.Builder.for_qubits(patch.used_set)

    builder.gate("R", patch.used_set)
    builder.tick()
    if basis == 'X':
        builder.gate("H", patch.data_set)
        builder.tick()

    do_measure_cycle(out=builder, patch=patch, save_layer='init')
    for tile in patch.tiles:
        if tile.basis == basis:
            builder.detector([gen.AtLayer(tile.measurement_qubit, 'init')], pos=tile.measurement_qubit)
    builder.shift_coords(dt=1)
    builder.tick()

    loop = builder.fork()
    do_measure_cycle(out=loop, patch=patch, save_layer='loop')
    for tile in patch.tiles:
        loop.detector([gen.AtLayer(tile.measurement_qubit, 'init'), gen.AtLayer(tile.measurement_qubit, 'loop')], pos=tile.measurement_qubit)
    loop.shift_coords(dt=1)
    loop.tick()
    builder.circuit += loop.circuit * (rounds - 1)

    if basis == 'X':
        builder.gate("H", patch.data_set)
        builder.tick()
    builder.measure(patch.data_set, save_layer='end')
    for tile in patch.tiles:
        if tile.basis == basis:
            builder.detector([gen.AtLayer(tile.measurement_qubit, 'loop')] + [gen.AtLayer(q, 'end') for q in tile.data_set], pos=tile.measurement_qubit)

    x_qs = {q for q in patch.data_set if q.imag == 0}
    z_qs = {q for q in patch.data_set if q.real == 0}
    assert len(x_qs & z_qs) % 2 == 1
    obs_qs = x_qs if basis == 'X' else z_qs
    builder.obs_include([gen.AtLayer(q, 'end') for q in obs_qs], obs_index=0)

    return CircuitCase(
        circuit=builder.circuit,
        patches=[patch],
        expected_interactions=frozenset(['CX']),
    )


def toric_patch(*, distance: int) -> gen.Patch:
    ur, ul, dl, dr = [(0.5 + 0.5j)*1j**d for d in range(4)]
    order_z = [ul, ur, dl, dr]
    order_ᴎ = [ul, dl, ur, dr]

    def wrap(q: complex) -> complex:
        r = q.real % distance
        i = q.imag % distance
        return r + 1j*i

    tiles = []
    for x in range(distance):
        for y in range(distance):
            m = x + 1j*y + 0.5 + 0.5j
            tile_basis = gen.checkerboard_basis(m)
            order = order_z if tile_basis == 'Z' else order_ᴎ
            tiles.append(gen.Tile(
                bases=tile_basis,
                measurement_qubit=m,
                ordered_data_qubits=[wrap(m + d) for d in order],
            ))

    return gen.Patch(tiles)


def do_measure_cycle(*,
                     patch: gen.Patch,
                     out: gen.Builder,
                     save_layer: str) -> None:
    out.gate("R", patch.measure_set)
    out.tick()
    out.gate("H", [tile.measurement_qubit for tile in patch.tiles if tile.basis == 'X'])
    out.tick()
    for k in range(4):
        for tile in patch.tiles:
            targets = (tile.measurement_qubit, tile.ordered_data_qubits[k])
            if tile.basis == 'Z':
                targets = targets[::-1]
            out.gate2("CX", [targets])
        out.tick()
    out.gate("H", [tile.measurement_qubit for tile in patch.tiles if tile.basis == 'X'])
    out.tick()
    out.measure(patch.measure_set, save_layer=save_layer)
