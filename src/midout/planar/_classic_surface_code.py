from typing import Union, Literal

from midout import gen
from midout._circuit_case import CircuitCase
from midout.gen._layer_translate import to_z_basis_interaction_circuit


def make_square_planar_cz_code(
        *,
        distance: int,
        basis: Union[Literal['X', 'Z'], str],
        rounds: int) -> CircuitCase:
    result = make_square_planar_cx_code(distance=distance, basis=basis, rounds=rounds)
    return CircuitCase(
        circuit=to_z_basis_interaction_circuit(result.circuit),
        patches=result.patches,
        expected_interactions=frozenset(['CZ']),
        show_patch_order=result.show_patch_order,
        show_patch_measure_qubits=result.show_patch_measure_qubits,
    )


def make_square_planar_cx_code(
        *,
        distance: int,
        basis: Union[Literal['X', 'Z'], str],
        rounds: int) -> CircuitCase:
    """Makes a standard surface code circuit."""
    assert rounds >= 1

    patch = gen.surface_code_patch(distance=distance)
    builder = gen.Builder.for_qubits(patch.used_set)

    rxs = {tile.measurement_qubit for tile in patch.tiles if tile.basis == 'X'}
    if basis == 'X':
        rxs |= patch.data_set
    rzs = patch.used_set - rxs
    builder.gate("R", rzs)
    builder.gate("RX", rxs)
    builder.tick()

    do_measure_cycle(out=builder, patch=patch, save_layer='init', measure_qubits_already_initialized=True)
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

    builder.measure(patch.data_set, basis=basis, save_layer='end')
    for tile in patch.tiles:
        if tile.basis == basis:
            builder.detector([gen.AtLayer(tile.measurement_qubit, 'loop')] + [gen.AtLayer(q, 'end') for q in tile.data_set], pos=tile.measurement_qubit)

    if basis == 'X':
        obs_qubits = {q for q in patch.data_set if q.imag == 0}
    else:
        obs_qubits = {q for q in patch.data_set if q.real == 0}
    builder.obs_include([gen.AtLayer(q, 'end') for q in obs_qubits], obs_index=0)

    return CircuitCase(
        circuit=builder.circuit,
        patches=[patch],
        expected_interactions=frozenset(['CX']),
        show_patch_order=True,
        show_patch_measure_qubits=True,
    )


def do_measure_cycle(*,
                     patch: gen.Patch,
                     out: gen.Builder,
                     measure_qubits_already_initialized: bool = False,
                     save_layer: str) -> None:
    if not measure_qubits_already_initialized:
        out.gate("R", [tile.measurement_qubit for tile in patch.tiles if tile.basis == 'Z'])
        out.gate("RX", [tile.measurement_qubit for tile in patch.tiles if tile.basis == 'X'])
        out.tick()
    for layer in range(4):
        cxs = []
        for tile in patch.tiles:
            q = tile.ordered_data_qubits[layer]
            if q is not None:
                cxs.append((q, tile.measurement_qubit)[::-1 if tile.basis == 'X' else +1])
        out.gate2('CX', cxs)
        out.tick()
    for basis in ['X', 'Z']:
        out.measure([tile.measurement_qubit for tile in patch.tiles if tile.basis == basis], basis=basis, save_layer=save_layer)
