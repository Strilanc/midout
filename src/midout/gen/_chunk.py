from typing import Iterable, Dict, Callable, Union

import sinter
import stim

from midout.gen._util import stim_circuit_with_transformed_coords
from midout.gen._flow import Flow, PauliString
from midout.gen._patch import Patch
from midout.gen._tile import Tile


class Chunk:
    """A circuit chunk with accompanying stabilizer flow assertions."""
    def __init__(self,
                 circuit: stim.Circuit,
                 q2i: Dict[complex, int],
                 flows: Iterable[Flow],
                 magic: bool = False,
                 discarded_inputs: Iterable[PauliString] = (),
                 discarded_outputs: Iterable[PauliString] = (),
                 repetitions: int = 1):
        self.q2i = q2i
        self.magic = magic
        self.circuit = circuit
        self.flows = tuple(flows)
        self.discarded_inputs = discarded_inputs
        self.discarded_outputs = discarded_outputs
        self.repetitions = repetitions

    def __eq__(self, other):
        if not isinstance(other, Chunk):
            return NotImplemented
        return (self.q2i == other.q2i and
                self.magic == other.magic and
                self.circuit == other.circuit and
                self.flows == other.flows and
                self.discarded_inputs == other.discarded_inputs and
                self.discarded_outputs == other.discarded_outputs and
                self.repetitions == other.repetitions)

    def with_flows_postselected(self, flow_predicate: Callable[[Flow], bool]) -> 'Chunk':
        return Chunk(
            circuit=self.circuit,
            q2i=self.q2i,
            magic=self.magic,
            flows=[flow.postselected() if flow_predicate(flow) else flow for flow in self.flows],
            discarded_inputs=self.discarded_inputs,
            discarded_outputs=self.discarded_outputs,
            repetitions=self.repetitions,
        )

    def with_repetitions(self, new_repetitions: int) -> 'Chunk':
        assert new_repetitions >= 0
        return Chunk(
            circuit=self.circuit,
            q2i=self.q2i,
            magic=self.magic,
            flows=self.flows,
            repetitions=new_repetitions,
            discarded_inputs=self.discarded_inputs,
            discarded_outputs=self.discarded_outputs,
        )

    def __mul__(self, other: int) -> 'Chunk':
        return self.with_repetitions(other)

    def verify(self):
        """Checks that this chunk's circuit actually implements its flows."""
        for key, group in sinter.group_by(self.flows, key=lambda flow: (flow.start, flow.obs_index)).items():
            if key[0] and len(group) > 1:
                raise ValueError(f"Multiple flows with same non-empty end: {group}")
        for key, group in sinter.group_by(self.flows, key=lambda flow: (flow.end, flow.obs_index)).items():
            if key[0] and len(group) > 1:
                raise ValueError(f"Multiple flows with same non-empty end: {group}")

        from midout.gen._flow_verifier import FlowStabilizerVerifier
        FlowStabilizerVerifier.verify(self)

        starts = {}
        ends = {}
        if self.repetitions != 1:
            for flow in self.flows:
                if flow.start:
                    starts[flow.start] = flow.obs_index
                if flow.end:
                    ends[flow.end] = flow.obs_index
            if starts != ends:
                raise ValueError("Not an exact loop.")

    def inverted(self) -> 'Chunk':
        """Checks that this chunk's circuit actually implements its flows."""
        from midout.gen._flow_verifier import FlowStabilizerVerifier
        return FlowStabilizerVerifier.invert(self)

    def with_xz_flipped(self) -> 'Chunk':
        return Chunk(
            q2i=self.q2i,
            magic=self.magic,
            circuit=circuit_with_xz_flipped(self.circuit),
            flows=[flow.with_xz_flipped() for flow in self.flows],
            discarded_inputs=[p.with_xz_flipped() for p in self.discarded_inputs],
            discarded_outputs=[p.with_xz_flipped() for p in self.discarded_outputs],
            repetitions=self.repetitions,
        )

    def with_transformed_coords(self, transform: Callable[[complex], complex]) -> 'Chunk':
        return Chunk(
            q2i={transform(q): i for q, i in self.q2i.items()},
            magic=self.magic,
            circuit=stim_circuit_with_transformed_coords(self.circuit, transform),
            flows=[flow.with_transformed_coords(transform) for flow in self.flows],
            discarded_inputs=[p.with_transformed_coords(transform) for p in self.discarded_inputs],
            discarded_outputs=[p.with_transformed_coords(transform) for p in self.discarded_outputs],
            repetitions=self.repetitions,
        )

    def magic_init_chunk(self) -> 'Chunk':
        """Returns a chunk that initializes the stabilizers needed by this one.

        The stabilizers are initialized using direct measurement by MPP, with
        no care for connectivity or physical limitations of hardware.
        """
        from midout.gen._flow_util import magic_init_for_chunk
        return magic_init_for_chunk(self)

    def magic_end_chunk(self) -> 'Chunk':
        """Returns a chunk that terminates the stabilizers produced by this one.

        The stabilizers are initialized using direct measurement by MPP, with
        no care for connectivity or physical limitations of hardware.
        """
        from midout.gen._flow_util import magic_measure_for_chunk
        return magic_measure_for_chunk(self)

    def _boundary_patch(self, end: bool) -> Patch:
        tiles = []
        for flow in self.flows:
            r = flow.end if end else flow.start
            if r.qubits and flow.obs_index is None:
                tiles.append(Tile(
                    ordered_data_qubits=r.qubits.keys(),
                    bases=''.join(r.qubits.values()),
                    measurement_qubit=list(r.qubits.keys())[0],
                ))
        return Patch(tiles)

    def start_patch(self) -> Patch:
        return self._boundary_patch(False)

    def end_patch(self) -> Patch:
        return self._boundary_patch(True)


class ChunkLoop:
    def __init__(self, chunks: Iterable[Union[Chunk, 'ChunkLoop']], repetitions: int):
        self.chunks = tuple(chunks)
        self.repetitions = repetitions

    def with_repetitions(self, new_repetitions: int) -> 'ChunkLoop':
        return ChunkLoop(chunks=self.chunks, repetitions=new_repetitions)

    def magic_init_chunk(self) -> 'Chunk':
        return self.chunks[0].magic_init_chunk()

    def magic_end_chunk(self) -> 'Chunk':
        return self.chunks[-1].magic_end_chunk()


XZ_FLIPPED = {
    "I": "I",
    "X": "Z",
    "Y": "Y",
    "Z": "X",
    "C_XYZ": "C_ZYX",
    "C_ZYX": "C_XYZ",
    "H": "H",
    "H_XY": "H_YZ",
    "H_XZ": "H_XZ",
    "H_YZ": "H_XY",
    "S": "SQRT_X",
    "SQRT_X": "S",
    "SQRT_X_DAG": "S_DAG",
    "SQRT_Y": "SQRT_Y",
    "SQRT_Y_DAG": "SQRT_Y_DAG",
    "S_DAG": "SQRT_X_DAG",
    "CX": "XCZ",
    "CY": "XCY",
    "CZ": "XCX",
    "ISWAP": None,
    "ISWAP_DAG": None,
    "SQRT_XX": "SQRT_ZZ",
    "SQRT_XX_DAG": "SQRT_ZZ_DAG",
    "SQRT_YY": "SQRT_YY",
    "SQRT_YY_DAG": "SQRT_YY_DAG",
    "SQRT_ZZ": "SQRT_XX",
    "SQRT_ZZ_DAG": "SQRT_XX_DAG",
    "SWAP": "SWAP",
    "XCX": "CZ",
    "XCY": "CY",
    "XCZ": "CX",
    "YCX": "YCZ",
    "YCY": "YCY",
    "YCZ": "YCX",
    "DEPOLARIZE1": "DEPOLARIZE1",
    "DEPOLARIZE2": "DEPOLARIZE2",
    "E": None,
    "ELSE_CORRELATED_ERROR": None,
    "PAULI_CHANNEL_1": None,
    "PAULI_CHANNEL_2": None,
    "X_ERROR": "Z_ERROR",
    "Y_ERROR": "Y_ERROR",
    "Z_ERROR": "X_ERROR",
    "M": "MX",
    "MPP": None,
    "MR": "MRX",
    "MRX": "MRZ",
    "MRY": "MRY",
    "MX": "M",
    "MY": "MY",
    "R": "RX",
    "RX": "R",
    "RY": "RY",
    "DETECTOR": "DETECTOR",
    "OBSERVABLE_INCLUDE": "OBSERVABLE_INCLUDE",
    "QUBIT_COORDS": "QUBIT_COORDS",
    "SHIFT_COORDS": "SHIFT_COORDS",
    "TICK": "TICK",
}


def circuit_with_xz_flipped(circuit: stim.Circuit) -> stim.Circuit:
    result = stim.Circuit()
    for inst in circuit:
        if isinstance(inst, stim.CircuitRepeatBlock):
            result.append(stim.CircuitRepeatBlock(
                body=circuit_with_xz_flipped(inst.body_copy()),
                repeat_count=inst.repeat_count))
        else:
            other = XZ_FLIPPED.get(inst.name)
            if other is None:
                raise NotImplementedError(f'{inst=}')
            result.append(stim.CircuitInstruction(other, inst.targets_copy(), inst.gate_args_copy()))
    return result
