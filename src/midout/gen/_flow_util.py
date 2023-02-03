from typing import Union, List, Tuple, Any, Optional, Dict, Literal, Iterable

import stim

from midout.gen._builder import MeasurementTracker, Builder, AtLayer
from midout.gen._chunk import Chunk, ChunkLoop
from midout.gen._flow import PauliString, Flow
from midout.gen._patch import Patch
from midout.gen._util import sorted_complex


def magic_init_for_chunk(
        chunk: Chunk,
) -> Chunk:
    builder = Builder(
        q2i=chunk.q2i,
        circuit=stim.Circuit(),
        tracker=MeasurementTracker(),
    )
    index = 0
    flows = []
    for flow in chunk.flows:
        if flow.start:
            builder.measure_pauli_product(q2b=flow.start.qubits, key=AtLayer(index, 'solo'))
            flows.append(Flow(
                center=flow.center,
                end=flow.start,
                measurement_indices=[index],
                obs_index=flow.obs_index,
            ))
            index += 1

    return Chunk(
        circuit=builder.circuit,
        q2i=builder.q2i,
        flows=flows,
        magic=True,
    )


def magic_measure_for_chunk(
        chunk: Chunk,
) -> Chunk:
    builder = Builder(
        q2i=chunk.q2i,
        circuit=stim.Circuit(),
        tracker=MeasurementTracker(),
    )
    index = 0
    flows = []
    for flow in chunk.flows:
        if flow.end:
            key = AtLayer(index, 'solo')
            builder.measure_pauli_product(q2b=flow.end.qubits, key=key)
            flows.append(Flow(
                center=flow.center,
                start=flow.end,
                measurement_indices=[index],
                obs_index=flow.obs_index,
            ))
            index += 1

    return Chunk(
        circuit=builder.circuit,
        q2i=builder.q2i,
        flows=flows,
        magic=True,
    )


def build_surface_code_round_circuit(
        patch: Patch,
        *,
        init_data_basis: Union[None, str, Dict[complex, str]] = None,
        measure_data_basis: Union[None, str, Dict[complex, str]] = None,
        save_layer: Any,
        out: Builder,
):
    measure_xs = Patch([tile for tile in patch.tiles if tile.basis == 'X'])
    measure_zs = Patch([tile for tile in patch.tiles if tile.basis == 'Z'])
    if init_data_basis is None:
        init_data_basis = {}
    elif isinstance(init_data_basis, str):
        init_data_basis = {q: init_data_basis for q in patch.data_set}
    if measure_data_basis is None:
        measure_data_basis = {}
    elif isinstance(measure_data_basis, str):
        measure_data_basis = {q: measure_data_basis for q in patch.data_set}

    out.gate("RX", measure_xs.measure_set)
    for basis in 'XYZ':
        qs = [q for q in init_data_basis if init_data_basis[q] == basis]
        if qs:
            out.gate(f"R{basis}", qs)
    out.gate("R", measure_zs.measure_set)
    out.tick()

    num_layers, = {len(tile.ordered_data_qubits) for tile in patch.tiles}
    for k in range(num_layers):
        out.gate2('CX', [
            (tile.measurement_qubit, tile.ordered_data_qubits[k])[::-1 if tile.basis == 'Z' else +1]
            for tile in patch.tiles
            if tile.ordered_data_qubits[k] is not None
        ])
        out.tick()

    out.measure(measure_xs.measure_set, basis='X', save_layer=save_layer)
    for basis in 'XYZ':
        qs = [q for q in measure_data_basis if measure_data_basis[q] == basis]
        if qs:
            out.measure(qs, basis=basis, save_layer=save_layer)
    out.measure(measure_zs.measure_set, basis='Z', save_layer=save_layer)


def standard_surface_code_chunk(
        patch: Patch,
        *,
        init_data_basis: Union[None, str, Dict[complex, str]] = None,
        measure_data_basis: Union[None, str, Dict[complex, str]] = None,
        obs: Optional[PauliString] = None,
) -> Chunk:
    if init_data_basis is None:
        init_data_basis = {}
    elif isinstance(init_data_basis, str):
        init_data_basis = {q: init_data_basis for q in patch.data_set}
    if measure_data_basis is None:
        measure_data_basis = {}
    elif isinstance(measure_data_basis, str):
        measure_data_basis = {q: measure_data_basis for q in patch.data_set}

    out = Builder.for_qubits(patch.used_set)
    save_layer = 'solo'
    build_surface_code_round_circuit(
        patch=patch,
        init_data_basis=init_data_basis,
        measure_data_basis=measure_data_basis,
        save_layer=save_layer,
        out=out,
    )

    discarded_inputs = []
    discarded_outputs = []
    flows = []
    for tile in patch.tiles:
        from_prev = PauliString({
            q: b
            for q, b in zip(tile.ordered_data_qubits, tile.bases)
            if q is not None and q not in init_data_basis
        })
        if any(init_data_basis.get(q, b) != b for q, b in zip(tile.ordered_data_qubits, tile.bases)):
            # Stabilizer anticommutes with a reset. Not prepared.
            if from_prev:
                discarded_inputs.append(from_prev)
            continue
        flows.append(
            Flow(
                center=tile.measurement_qubit,
                start=from_prev,
                measurement_indices=out.tracker.measurement_indices([AtLayer(tile.measurement_qubit, save_layer)]),
            )
        )

    for tile in patch.tiles:
        to_next = PauliString({
            q: b
            for q, b in zip(tile.ordered_data_qubits, tile.bases)
            if q is not None and q not in measure_data_basis
        })
        if any(measure_data_basis.get(q, b) != b for q, b in zip(tile.ordered_data_qubits, tile.bases)):
            # Stabilizer anticommutes with a reset. Not prepared.
            if to_next:
                discarded_outputs.append(to_next)
            continue
        flows.append(
            Flow(
                center=tile.measurement_qubit,
                end=to_next,
                measurement_indices=out.tracker.measurement_indices([
                    AtLayer(q, save_layer)
                    for q in tile.used_set
                    if q in measure_data_basis or q == tile.measurement_qubit
                ]),
            )
        )

    if obs is not None:
        start_obs = dict(obs.qubits)
        end_obs = dict(obs.qubits)
        for q in init_data_basis:
            if q in start_obs:
                if start_obs.pop(q) != init_data_basis[q]:
                    raise ValueError("wrong init basis for obs")
        measure_indices = []
        for q in measure_data_basis:
            if q in end_obs:
                if end_obs.pop(q) != measure_data_basis[q]:
                    raise ValueError("wrong measure basis for obs")
                measure_indices.extend(out.tracker.measurement_indices([
                    AtLayer(q, save_layer)
                ]))

        flows.append(
            Flow(
                center=0,
                start=PauliString(start_obs),
                end=PauliString(end_obs),
                obs_index=0,
                measurement_indices=measure_indices,
            )
        )

    return Chunk(
        circuit=out.circuit,
        q2i=out.q2i,
        flows=flows,
        discarded_inputs=discarded_inputs,
        discarded_outputs=discarded_outputs,
    )


def relabel_circuit_into(*, circuit: stim.Circuit, old_q2i: Dict[complex, int], new_q2i: Dict[complex, int], out: stim.Circuit):
    i2i = {i: new_q2i[q] for q, i in old_q2i.items()}

    for inst in circuit:
        if inst.name == 'QUBIT_COORDS':
            continue
        targets = []
        for t in inst.targets_copy():
            if t.is_qubit_target:
                targets.append(i2i[t.value])
            elif t.is_x_target:
                targets.append(stim.target_x(i2i[t.value]))
            elif t.is_y_target:
                targets.append(stim.target_y(i2i[t.value]))
            elif t.is_z_target:
                targets.append(stim.target_z(i2i[t.value]))
            elif t.is_combiner:
                targets.append(t)
            elif t.is_measurement_record_target:
                targets.append(t)
            else:
                raise NotImplementedError(f'{inst=}')
        out.append(inst.name, targets, inst.gate_args_copy())


class ChunkCompileState:
    def __init__(self, *, open_flows: Dict[Tuple[PauliString, Any], Union[Flow, Literal["discard"]]], measure_offset: int):
        self.open_flows = open_flows
        self.measure_offset = measure_offset


def _compile_chunk_into_circuit_many_repetitions(
        *,
        chunk: Union[Chunk, ChunkLoop],
        state: ChunkCompileState,
        include_detectors: bool,
        ignore_errors: bool,
        out_circuit: stim.Circuit,
        q2i: Dict[complex, int],
) -> ChunkCompileState:
    assert chunk.repetitions > 1
    no_reps = chunk.with_repetitions(1)
    circuits = []
    measure_offset_start_of_loop = state.measure_offset
    while len(circuits) < chunk.repetitions:
        fully_in_loop = min([
            m
            for flow in state.open_flows.values()
            if isinstance(flow, Flow)
            for m in flow.measurement_indices
        ], default=measure_offset_start_of_loop) >= measure_offset_start_of_loop

        circuits.append(stim.Circuit())
        state = compile_chunk_into_circuit(
            chunk=no_reps,
            state=state,
            include_detectors=include_detectors,
            ignore_errors=ignore_errors,
            out_circuit=circuits[-1],
            q2i=q2i,
        )

        if fully_in_loop:
            # The circuit is guaranteed to repeat now. Don't do each iteration individually.
            finish_reps = chunk.repetitions - len(circuits) + 1
            while len(circuits) > 1 and circuits[-1] == circuits[-2]:
                finish_reps += 1
                circuits.pop()
            circuits[-1] *= finish_reps
            break

    # Fuse iterations that happened to be equal.
    k = 0
    while k < len(circuits):
        k2 = k + 1
        while k2 < len(circuits) and circuits[k2] == circuits[k]:
            k2 += 1
        out_circuit += circuits[k] * (k2 - k)
        k = k2

    return state


def _compile_chunk_into_circuit_sequence(
        *,
        chunks: Iterable[Union[Chunk, ChunkLoop]],
        state: ChunkCompileState,
        include_detectors: bool,
        ignore_errors: bool,
        out_circuit: stim.Circuit,
        q2i: Dict[complex, int],
) -> ChunkCompileState:
    for sub_chunk in chunks:
        state = compile_chunk_into_circuit(
            chunk=sub_chunk,
            state=state,
            include_detectors=include_detectors,
            ignore_errors=ignore_errors,
            out_circuit=out_circuit,
            q2i=q2i,
        )
    return state


def _compile_chunk_into_circuit_atomic(
        *,
        chunk: Chunk,
        state: ChunkCompileState,
        include_detectors: bool,
        ignore_errors: bool,
        out_circuit: stim.Circuit,
        q2i: Dict[complex, int],
) -> ChunkCompileState:
    prev_flows = dict(state.open_flows)
    next_flows: Dict[Tuple[PauliString, Any], Union[Flow, Literal['discard']]] = {}
    dumped_flows: List[Flow] = []
    if include_detectors:
        for flow in chunk.flows:
            flow = Flow(
                center=flow.center,
                start=flow.start,
                end=flow.end,
                obs_index=flow.obs_index,
                measurement_indices=[m + state.measure_offset for m in flow.measurement_indices],
                postselect=flow.postselect,
            )
            if flow.start:
                prev = prev_flows.pop((flow.start, flow.obs_index), None)
                if prev is None:
                    if ignore_errors:
                        continue
                    else:
                        raise ValueError(f"Missing prev {flow!r} have {prev_flows!r}")
                elif prev == 'discard':
                    if flow.end:
                        next_flows[(flow.end, flow.obs_index)] = 'discard'
                    continue
                flow = prev.concat(flow, 0)
            if flow.end:
                if flow.obs_index is not None and flow.measurement_indices:
                    dumped_flows.append(flow)
                    flow = Flow(start=flow.start, end=flow.end, obs_index=flow.obs_index, center=flow.center)
                next_flows[(flow.end, flow.obs_index)] = flow
            else:
                dumped_flows.append(flow)
        for discarded in chunk.discarded_inputs:
            prev_flows.pop((discarded, None), None)
        for discarded in chunk.discarded_outputs:
            assert (discarded, None) not in next_flows
            next_flows[(discarded, None)] = "discard"
        for flow, val in prev_flows.items():
            if val != "discard" and not ignore_errors:
                raise ValueError(f"Some flows were left over (not matched) when moving into chunk: {list(prev_flows.values())!r}")

    new_measure_offset = state.measure_offset + chunk.circuit.num_measurements
    relabel_circuit_into(circuit=chunk.circuit, out=out_circuit, old_q2i=chunk.q2i, new_q2i=q2i)
    if include_detectors:
        any_detectors = False
        for flow in dumped_flows:
            targets = []
            for m in flow.measurement_indices:
                targets.append(stim.target_rec(m - new_measure_offset))
            if flow.obs_index is None:
                coords = (flow.center.real, flow.center.imag, 0)
                if flow.postselect:
                    coords += (999,)
                out_circuit.append("DETECTOR", targets, coords)
                any_detectors = True
            else:
                out_circuit.append("OBSERVABLE_INCLUDE", targets, flow.obs_index)
        if any_detectors:
            out_circuit.append("SHIFT_COORDS", [], (0, 0, 1))
    out_circuit.append("TICK")

    return ChunkCompileState(
        measure_offset=new_measure_offset,
        open_flows=next_flows,
    )


def compile_chunk_into_circuit(
    *,
    chunk: Union[Chunk, ChunkLoop],
    state: ChunkCompileState,
    include_detectors: bool,
    ignore_errors: bool,
    out_circuit: stim.Circuit,
    q2i: Dict[complex, int],
) -> ChunkCompileState:
    if chunk.repetitions == 0:
        return state
    if chunk.repetitions > 1:
        return _compile_chunk_into_circuit_many_repetitions(
            chunk=chunk,
            state=state,
            include_detectors=include_detectors,
            ignore_errors=ignore_errors,
            out_circuit=out_circuit,
            q2i=q2i,
        )
    if isinstance(chunk, ChunkLoop):
        return _compile_chunk_into_circuit_sequence(
            chunks=chunk.chunks,
            state=state,
            include_detectors=include_detectors,
            ignore_errors=ignore_errors,
            out_circuit=out_circuit,
            q2i=q2i,
        )

    return _compile_chunk_into_circuit_atomic(
        chunk=chunk,
        state=state,
        include_detectors=include_detectors,
        ignore_errors=ignore_errors,
        out_circuit=out_circuit,
        q2i=q2i,
    )


def compile_chunks_into_circuit(
        chunks: List[Union[Chunk, ChunkLoop]],
        *,
        include_detectors: bool = True,
        ignore_errors: bool = False,
) -> stim.Circuit:
    all_qubits = set()

    def _process(c: Union[Chunk, ChunkLoop]):
        nonlocal all_qubits
        if isinstance(c, ChunkLoop):
            for c2 in c.chunks:
                _process(c2)
        elif isinstance(c, Chunk):
            all_qubits |= c.q2i.keys()
        else:
            raise NotImplementedError(f'{c=}')
    for c in chunks:
        _process(c)

    q2i = {q: i for i, q in enumerate(sorted_complex(set(all_qubits)))}
    full_circuit = stim.Circuit()
    for q, i in q2i.items():
        full_circuit.append('QUBIT_COORDS', i, [q.real, q.imag])

    state = ChunkCompileState(open_flows={}, measure_offset=0)
    for k, chunk in enumerate(chunks):
        state = compile_chunk_into_circuit(
            chunk=chunk,
            state=state,
            include_detectors=include_detectors,
            ignore_errors=ignore_errors,
            out_circuit=full_circuit,
            q2i=q2i,
        )
    if include_detectors:
        if state.open_flows:
            if not ignore_errors:
                raise ValueError("Unterminated")
    return full_circuit


def verify_circuit_has_all_possible_detectors(circuit: stim.Circuit):
    """Checks that the number of predictable measurements is equal to dets+obs.
    """

    num_declarations = circuit.num_detectors + circuit.num_observables
    num_determined_measurements = 0
    num_detectors_seen = 0
    sim = stim.TableauSimulator()

    tick = 0
    def run_block(block: stim.Circuit, reps: int):
        nonlocal num_determined_measurements, tick, num_detectors_seen
        for _ in range(reps):
            for inst in block:
                if isinstance(inst, stim.CircuitRepeatBlock):
                    run_block(inst.body_copy(), inst.repeat_count)
                elif inst.name == 'DETECTOR':
                    num_detectors_seen += 1
                elif inst.name == 'TICK':
                    tick += 1
                elif inst.name == 'M' or inst.name == 'MR':
                    args = inst.gate_args_copy()
                    for t in inst.targets_copy():
                        assert t.is_qubit_target
                        known = sim.peek_z(t.value) != 0
                        num_determined_measurements += known
                        sim.do(stim.CircuitInstruction(inst.name, [t.value], args))
                elif inst.name == 'MX' or inst.name == 'MRX':
                    args = inst.gate_args_copy()
                    for t in inst.targets_copy():
                        assert t.is_qubit_target
                        known = sim.peek_x(t.value) != 0
                        num_determined_measurements += known
                        sim.do(stim.CircuitInstruction(inst.name, [t.value], args))
                elif inst.name == 'MY' or inst.name == 'MRY':
                    args = inst.gate_args_copy()
                    for t in inst.targets_copy():
                        assert t.is_qubit_target
                        known = sim.peek_y(t.value) != 0
                        num_determined_measurements += known
                        sim.do(stim.CircuitInstruction(inst.name, [t.value], args))
                elif inst.name == 'MPP':
                    raise NotImplementedError(f'{inst=}')
                else:
                    sim.do(inst)
    run_block(circuit, 1)
    if num_declarations != num_determined_measurements:
        raise ValueError(f"{num_declarations=} != {num_determined_measurements=}")
