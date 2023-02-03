import dataclasses
import itertools
from typing import Dict, Iterator, Optional

import stim

from midout.walking.cnot_layer import CnotLayer
from midout.walking.meas_layer import MeasLayer
from midout.walking.reset_layer import ResetLayer
from midout.walking.util import Basis, Measurement, Qubit

STIM_RESET = {Basis.z: "RZ", Basis.x: "RX"}
STIM_MEASURE = {Basis.z: "MZ", Basis.x: "MX"}
STIM_NONCOMMUTING_ERRORS = {Basis.z: "X_ERROR", Basis.x: "Z_ERROR"}
STIM_OBSV_INDEX = {Basis.z: 0, Basis.x: 1}


@dataclasses.dataclass
class StimCircuitBuilder:
    """Builds a stim circuit, including tracking all the things you don't want to track yourself.

    You probably just want to use `build_circuit_from_layers`

    Otherwise, just make a StimCircuitBuilder() and then process layers as you like
        Everything internally should be initialised for you,
        Notice that TICKs are added in `process_layer`

    Will add very simple errors if you ask it too,
        2Q depolarizing immediately after every CNOT
        1Q flip after every reset (X for RZ and Z for RX)
        1Q flips immediately before measurement (X for MZ and Z for MX)
    This is really aimed at covering every part of the circuit with errors,
        so that shortest_graphlike_error can work nicely
    For any proper error model, you should add your own errors to the clean circuit
    """

    qubits: Dict[Qubit, int] = dataclasses.field(default_factory=dict)
    next_qubit_index: int = 0
    measurements: Dict[Measurement, int] = dataclasses.field(default_factory=dict)
    next_measurement_index: int = 0
    _stim_circuit: stim.Circuit = dataclasses.field(default_factory=stim.Circuit)
    _qubit_coords_circuit: stim.Circuit = dataclasses.field(default_factory=stim.Circuit)

    errors: Optional[float] = None
    duid_counter: Iterator = itertools.count()

    @property
    def stim_circuit(self):
        return self._qubit_coords_circuit + self._stim_circuit

    @staticmethod
    def build_circuit_from_layers(layers, errors=None) -> stim.Circuit:
        scb = StimCircuitBuilder(errors=errors)
        for l in layers:
            scb.process_layer(l)
        return scb._stim_circuit

    def add_qubit(self, qubit: Qubit):
        self.qubits[qubit] = self.next_qubit_index
        self.next_qubit_index += 1

    def add_measurement(self, measurement: 'Measurement'):
        self.measurements[measurement] = self.next_measurement_index
        self.next_measurement_index += 1

    def measurement_index(self, measurement):
        return self.measurements[measurement] - self.next_measurement_index

    def process_layer(self, layer):
        if isinstance(layer, ResetLayer):
            self.process_reset_layer(layer)
        elif isinstance(layer, CnotLayer):
            self.process_cnot_layer(layer)
        elif isinstance(layer, MeasLayer):
            self.process_meas_layer(layer)
        else:
            raise ValueError(f"Unrecognised Layer: {layer, type(layer)}")
        self._stim_circuit.append("TICK")

    def process_reset_layer(self, layer: ResetLayer):
        """process a reset layer, which comes down to finding new qubits and adding QUBIT_COORDS."""
        for b, qubits in sorted(layer.resets.items(), key=lambda kvpair: kvpair[0].value):
            if b is None:
                raise ValueError("Can't build stim circuit for non-basis experiments")
            qubit_indices = []
            for q in sorted(qubits, key=lambda q: (q.real, q.imag)):
                if q not in self.qubits:
                    self.add_qubit(q)
                    self._qubit_coords_circuit.append("QUBIT_COORDS", [self.qubits[q]], [q.real, q.imag])
                qubit_indices.append(self.qubits[q])
            self._stim_circuit.append(STIM_RESET[b], qubit_indices)
            if self.errors:
                self._stim_circuit.append(STIM_NONCOMMUTING_ERRORS[b], qubit_indices, self.errors)

    def process_cnot_layer(self, layer: CnotLayer):
        """process a cnot layer, which is really just adding the right cnots.

        first qubit is control, second qubit is target
        """
        qubit_indices = []
        for t, c in sorted(layer.targets.items(), key=lambda kvpair: self.qubits[kvpair[0]]):
            qubit_indices.append(self.qubits[t])
            qubit_indices.append(self.qubits[c])
        self._stim_circuit.append("CX", qubit_indices)
        if self.errors:
            self._stim_circuit.append("DEPOLARIZE2", qubit_indices, self.errors)

    def process_meas_layer(self, layer: MeasLayer):
        """process a measurement layer, which means adding M, DETECTOR and OBVSERVABLE_INCLUDE."""
        for b, meas in layer.measurements.items():
            qubit_indices = []
            for m in sorted(meas, key=lambda m: (m.index, self.qubits[m.qubit])):
                qubit_indices.append(self.qubits[m.qubit])
                self.add_measurement(m)
            if len(qubit_indices) != 0:
                if self.errors:
                    self._stim_circuit.append(
                        STIM_NONCOMMUTING_ERRORS[b], qubit_indices, self.errors
                    )
                self._stim_circuit.append(STIM_MEASURE[b], qubit_indices)

        det_measurements_by_coords = {}
        for d, duid in sorted(layer.detectors):
            meas_indices = []
            last_m = None
            for m in sorted(d, key=lambda m: (m.index, self.qubits[m.qubit])):
                last_m = m if last_m is None or m.index > last_m.index else last_m
                meas_indices.append(self.measurement_index(m))
            if len(meas_indices) != 0:
                coords = [last_m.qubit.real, last_m.qubit.imag, last_m.index]
                if duid:
                    coords += [duid]
                det_measurements_by_coords[tuple(coords)] = meas_indices

        for coords, meas_indices in sorted(det_measurements_by_coords.items()):
            self._stim_circuit.append(
                "DETECTOR", [stim.target_rec(i) for i in meas_indices], coords
            )

        for b, obvs_incs in sorted(layer.observable_includes.items(), key=lambda kvpair:(kvpair[0].value)):
            meas_indices = [self.measurement_index(m) for m in obvs_incs]
            if len(meas_indices) != 0:
                self._stim_circuit.append(
                    "OBSERVABLE_INCLUDE",
                    [stim.target_rec(i) for i in meas_indices],
                    STIM_OBSV_INDEX[b],
                )

    def next_duid(self):
        return next(self.duid_counter)
