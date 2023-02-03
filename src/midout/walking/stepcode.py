import itertools
from typing import List

import stim


def reset_layer(first_qubit, last_qubit):
    return "R " + " ".join(str(i) for i in range(first_qubit, last_qubit + 1))


def cnot_layer(first_qubit, last_qubit, direction: str):
    """make a cnot layer for the stepcode."""

    qubits = list(range(first_qubit, last_qubit + 1))

    # if direction is up, reverse the whole qubit list
    qubits = qubits[::-1] if direction == 'up' else qubits

    gate_args = [(a, b) for a, b in zip(qubits[::2], qubits[1::2])]

    flatten_gate_args = list(itertools.chain.from_iterable(gate_args))

    return ["CNOT " + " ".join(str(i) for i in flatten_gate_args), "TICK"]


def measure_layer(
    first_qubit, last_qubit, lost_qubit=None, lost_qubit_last_round=False, first_round=False
):
    return_strings = []

    extra_qubit_meas_string = None
    if lost_qubit is not None:
        return_strings.append(f"H {lost_qubit}")
        return_strings.append("TICK")
        extra_qubit_meas_string = f"M {lost_qubit}"

    return_strings.append("MR " + " ".join(str(i) for i in range(first_qubit + 1, last_qubit, 2)))
    if extra_qubit_meas_string:
        return_strings.append(extra_qubit_meas_string)

    num_measure_qubits = int((last_qubit - first_qubit) / 2)
    offset_last_round = 1 if lost_qubit_last_round else 0
    offset_all_recs = 1 if lost_qubit is not None else 0
    if first_round:
        detectors = [f"DETECTOR rec[{-i - 1 - offset_all_recs}]" for i in range(num_measure_qubits)]
    else:
        detectors = [
            f"DETECTOR "
            f"rec[{-i - 1 - num_measure_qubits - offset_last_round - offset_all_recs}] "
            f"rec[{-i - 1 - offset_all_recs}]"
            for i in range(num_measure_qubits)
        ]
    return_strings.extend(detectors)
    return_strings.append("TICK")

    return return_strings


def last_measure_layer(first_qubit, last_qubit, lost_qubit_last_round=False):
    return_strings = []
    # measure all the data qubits
    return_strings.append("M " + " ".join([str(q) for q in range(first_qubit, last_qubit + 1, 2)]))

    offset = 1 if lost_qubit_last_round else 0
    num_data_qubits = int((last_qubit - first_qubit) / 2 + 1)
    return_strings.extend(
        [
            f"DETECTOR rec[{-i - 1 - num_data_qubits - offset}] rec[{-i - 2}] rec[{-i - 1}]"
            for i in range(num_data_qubits - 1)
        ]
    )

    return_strings.append("OBSERVABLE_INCLUDE(0) rec[-1]")
    return return_strings


def build_stepcode_circuit(
    first_qubit: int, last_qubit: int, step_sequence: List[str], errors=False
):
    """

    Assume you're starting with an actual repcode state over alternating measure and data qubits
        if your initial qubits aren't sane, no promise this will work

    Args:
        first_qubit: index of the first qubit at the start of the code
        last_qubit: index of the last qubit at the start of the code
        step_sequence: sequence of strings indicating the cnot layers to build
            'up' 'down' means a normal rep code round
            'down' 'up' also makes a regular round, but with the CNOT layers switched,
            'down' 'down' makes a round that steps the measure qubits down by 1 qubit
            'up' 'up' makes a round that moves the measure qubits up by 1 qubit
        errors: if False or None, doesn't add any errors
            if number like, applies an X error to every qubit in the code right before measurement
            with the given probability

    """
    assert abs(first_qubit - last_qubit) % 2 == 0  # make sure there's an odd number of qubits

    full_circuit_list = [reset_layer(first_qubit, last_qubit), "TICK"]
    lost_qubit_last_round = False
    first_round = True

    for d0, d1 in zip(step_sequence[::2], step_sequence[1::2]):

        # figure out where we should put gates
        if (d0, d1) == ('up', 'up'):
            extra_qubit = first_qubit - 1
            full_circuit_list.extend([f"R {extra_qubit}", "TICK"])
            full_circuit_list.extend(cnot_layer(first_qubit - 1, last_qubit, direction=d0))
            full_circuit_list.extend(cnot_layer(first_qubit, last_qubit - 1, direction=d1))
            # now shift the qubits for the measurements
            lost_qubit = last_qubit
            first_qubit, last_qubit = first_qubit - 1, last_qubit - 1

        elif (d0, d1) == ('down', 'down'):
            extra_qubit = last_qubit + 1
            full_circuit_list.extend([f"R {extra_qubit}", "TICK"])
            full_circuit_list.extend(cnot_layer(first_qubit, last_qubit + 1, direction=d0))
            full_circuit_list.extend(cnot_layer(first_qubit + 1, last_qubit, direction=d1))

            lost_qubit = first_qubit
            first_qubit, last_qubit = first_qubit + 1, last_qubit + 1
        else:
            # regular rep code, no extra qubit, no lost qubit
            full_circuit_list.extend(cnot_layer(first_qubit, last_qubit, direction=d0))
            full_circuit_list.extend(cnot_layer(first_qubit, last_qubit, direction=d1))
            lost_qubit = None

        if errors:
            full_circuit_list.append(
                f"X_ERROR({errors}) "
                + " ".join(str(q) for q in range(first_qubit, last_qubit + 1))
                + ("" if lost_qubit is None else f" {lost_qubit}")
            )  # don't tick before measurement

        full_circuit_list.extend(
            measure_layer(first_qubit, last_qubit, lost_qubit, lost_qubit_last_round, first_round)
        )

        first_round = False
        lost_qubit_last_round = lost_qubit is not None

    full_circuit_list.extend(last_measure_layer(first_qubit, last_qubit, lost_qubit_last_round))

    return stim.Circuit("\n".join(full_circuit_list))


def normal_rep_code(distance: int, rounds: int, errors=None):
    num_qubits = distance * 2 - 1
    return build_stepcode_circuit(
        first_qubit=0,
        last_qubit=num_qubits - 1,
        step_sequence=['down', 'up'] * rounds,
        errors=errors,
    )


def other_normal_rep_code(distance: int, rounds: int, errors=None):
    num_qubits = distance * 2 - 1
    return build_stepcode_circuit(
        first_qubit=0,
        last_qubit=num_qubits - 1,
        step_sequence=['up', 'down'] * rounds,
        errors=errors,
    )


def strange_rep_code(distance: int, rounds: int, errors=None):
    num_qubits = distance * 2 - 1
    return build_stepcode_circuit(
        first_qubit=0,
        last_qubit=num_qubits - 1,
        step_sequence=(
            ['up', 'down', 'down', 'up'] * int(rounds / 2)
            + (['up', 'down'] if rounds % 2 == 1 else [])
        ),
        errors=errors,
    )


def gliding_step_code(distance, rounds, errors=None):
    num_qubits = distance * 2 - 1
    return build_stepcode_circuit(
        first_qubit=0,
        last_qubit=num_qubits - 1,
        step_sequence=(['down', 'down'] * rounds),
        errors=errors,
    )


def wiggling_step_code(distance, rounds, errors=None):
    num_qubits = distance * 2 - 1
    return build_stepcode_circuit(
        first_qubit=0,
        last_qubit=num_qubits - 1,
        step_sequence=(
            ['down', 'down', 'up', 'up'] * int(rounds / 2)
            + (['down', 'down'] if rounds % 2 == 1 else [])
        ),
        errors=errors,
    )
