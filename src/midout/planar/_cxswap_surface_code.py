from typing import Optional, Literal, Union, Any, List, Tuple

from midout import gen
from midout._circuit_case import CircuitCase


def make_cx_swap_surface_code_patch(*, distance: int) -> gen.Patch:
    """Creates the stabilizer tiles of the mid-cycle state when using CXSWAP."""
    data_qubits = {
        (x + 1j*y)*0.5
        for x in range(-1, distance*2)
        for y in range(-1, distance*2)
        if x % 2 == y % 2
        if x != -1 or y != -1
        if x != 2*distance-1 or y != -1
        if y != 2*distance-1 or x != -1
        if x != 2*distance-1 or y != 2*distance-1
    }
    tiles = []
    ds = [0, 0.5 + 0.5j, 1j, -0.5 + 0.5j]
    for q in data_qubits:
        qs = [q + d for d in ds if q + d in data_qubits]
        if len(qs) == 4:
            b = 'X' if q.real % 1 == 0 else 'Z'
            if q.real % 1 == 0 and (q.real - q.imag) % 2 == 0:
                m = q + -0.5 + 0.5j
            elif q.real % 1 == 0 and (q.real - q.imag) % 2 == 1:
                m = q + 0.5 + 0.5j
            elif q.real % 1 == 0.5 and (q.real - q.imag) % 2 == 1:
                m = q + 1j
            else:
                m = q
            tiles.append(gen.Tile(
                bases=b,
                measurement_qubit=m,
                ordered_data_qubits=qs,
            ))
    for x in range(distance - 1):
        for q in [x + 0.5 - 0.5j, x + 0.5 - 0.5j + distance*1j]:
            tiles.append(gen.Tile(
                bases='Z',
                measurement_qubit=q,
                ordered_data_qubits=[q],
            ))
            q2 = q.imag + q.real*1j
            tiles.append(gen.Tile(
                bases='X',
                measurement_qubit=q2,
                ordered_data_qubits=[q2],
            ))
    return gen.Patch(tiles)


def make_cx_swap_surface_code_chunk(
        *,
        distance: int,
        basis: str,
        round_parity: bool,
        is_first_round: bool,
) -> gen.Chunk:
    """Creates a chunk that measures half the stabilizers of the code.

    The start state and end state are identical, except that half the
    stabilizers were measured. To measure the other half, change `round_parity`.
    """

    patch = make_cx_swap_surface_code_patch(distance=distance)
    builder = gen.Builder.for_qubits(patch.used_set)

    active_qubits = set(patch.used_set)
    measured_tiles = []
    for tile in patch.tiles:
        d = tile.ordered_data_qubits[0]
        b = (d.real + d.imag) % 2 == round_parity
        if len(tile.ordered_data_qubits) == 4 and b:
            measured_tiles.append(tile)
        if len(tile.ordered_data_qubits) == 1:
            deactivate = False

            if d.real == distance-0.5 and b:
                measured_tiles.append(tile)
                deactivate = round_parity
            if d.imag == -0.5 and not b:
                measured_tiles.append(tile)
                deactivate = round_parity
            if d.imag == distance-0.5 and b:
                measured_tiles.append(tile)
                deactivate = not round_parity
            if d.real == -0.5 and not b:
                measured_tiles.append(tile)
                deactivate = not round_parity

            if deactivate:
                active_qubits.remove(d)

    def towards(source_basis: str, delta: complex, sign: int) -> List[Tuple[complex, complex]]:
        """Finds two qubit gate locations."""
        if round_parity:
            delta *= -1
        mx = {m for m in patch.measure_set if (m.real + m.imag) % 2 == 0}
        mz = {m for m in patch.measure_set if (m.real + m.imag) % 2 != 0}
        ms = mx if source_basis == 'X' else mz
        result = []
        for start in ms:
            partner = start + delta
            if start in active_qubits and partner in active_qubits:
                result.append((start, partner)[::sign])
        return result

    layer1 = towards('X', 0.5 + 0.5j, -1) + towards('Z', -0.5 - 0.5j, +1)
    layer2 = towards('X', 0.5 - 0.5j, -1) + towards('Z', 0.5 - 0.5j, +1)
    measure_xs = {tile.measurement_qubit for tile in measured_tiles if tile.basis == 'X'}
    measure_zs = {tile.measurement_qubit for tile in measured_tiles if tile.basis == 'Z'}
    if is_first_round:
        builder.gate('RX', measure_xs)
        builder.gate(f'R{basis}', patch.used_set - measure_xs - measure_zs)
        builder.gate('R', measure_zs)
    else:
        for layer in [layer1, layer2]:
            builder.gate2('CXSWAP', layer)
            builder.tick()
        builder.demolition_measure_with_feedback_passthrough(
            xs=measure_xs,
            zs=measure_zs,
            save_layer='solo',
        )

    for layer in [layer2, layer1]:
        builder.tick()
        builder.gate2('SWAPCX', layer)

    flows = []
    discarded_outputs = []
    for tile in patch.tiles:
        t = gen.PauliString.from_tile_data(tile)
        if tile in measured_tiles:
            measurement_indices = []
            if not is_first_round:
                measurement_indices = builder.tracker.measurement_indices([gen.AtLayer(tile.measurement_qubit, 'solo')])
                flows.append(gen.Flow(
                    start=t,
                    measurement_indices=measurement_indices,
                    center=tile.measurement_qubit,
                ))
            flows.append(gen.Flow(
                end=t,
                measurement_indices=measurement_indices,
                center=tile.measurement_qubit,
            ))
        else:
            if not is_first_round or tile.basis == basis:
                flows.append(gen.Flow(
                    start=None if is_first_round else t,
                    end=t,
                    center=tile.measurement_qubit,
                ))
            else:
                discarded_outputs.append(t)
    obs_x = gen.PauliString({q: 'X' for q in patch.used_set if q.imag == 0})
    obs_z = gen.PauliString({q: 'Z' for q in patch.used_set if q.real == 0})
    assert obs_x.anticommutes(obs_z)
    obs = obs_x if basis == 'X' else obs_z
    flows.append(gen.Flow(
        start=None if is_first_round else obs,
        end=obs,
        obs_index=0,
        center=0,
    ))

    return gen.Chunk(
        circuit=builder.circuit,
        q2i=builder.q2i,
        flows=flows,
        discarded_outputs=discarded_outputs,
    )


def make_square_planar_cxswap_code(
        *,
        distance: int,
        basis: Union[Literal['X', 'Z'], str],
        rounds: int,
        use_iswaps: bool = False,
) -> CircuitCase:
    """Creates a full 4-CXSWAP memory experiment from init to measure."""
    assert rounds >= 2
    chunks = []

    chunks.append(make_cx_swap_surface_code_chunk(
        distance=distance,
        basis=basis,
        round_parity=False,
        is_first_round=True,
    ))
    chunks.append(make_cx_swap_surface_code_chunk(
        distance=distance,
        basis=basis,
        round_parity=True,
        is_first_round=False,
    ))

    chunks.append(gen.ChunkLoop(
        chunks=[
            make_cx_swap_surface_code_chunk(
                distance=distance,
                basis=basis,
                round_parity=False,
                is_first_round=False,
            ),
            make_cx_swap_surface_code_chunk(
                distance=distance,
                basis=basis,
                round_parity=True,
                is_first_round=False,
            ),
        ],
        repetitions=(rounds - 2) // 2,
    ))

    if rounds % 2 == 1:
        chunks.append(make_cx_swap_surface_code_chunk(
            distance=distance,
            basis=basis,
            round_parity=False,
            is_first_round=False,
        ))
    chunks.append(make_cx_swap_surface_code_chunk(
        distance=distance,
        basis=basis,
        round_parity=rounds % 2 == 0,
        is_first_round=True,
    ).inverted())

    circuit = gen.compile_chunks_into_circuit(chunks).with_inlined_feedback()
    if use_iswaps:
        circuit = gen.to_z_basis_interaction_circuit(circuit)

    return CircuitCase(
        circuit=circuit,
        patches=[make_cx_swap_surface_code_patch(distance=distance)],
        expected_interactions=frozenset(['ISWAP']) if use_iswaps else frozenset(['CXSWAP']),
        show_patch_order=False,
        show_patch_measure_qubits=True,
    )