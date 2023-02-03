import pytest

from midout.walking.stepcode import (
    cnot_layer,
    gliding_step_code,
    last_measure_layer,
    measure_layer,
    normal_rep_code,
    other_normal_rep_code,
    reset_layer,
    strange_rep_code,
    wiggling_step_code,
)


class TestStepCode:
    def test_reset_layer(self):
        assert reset_layer(0, 10) == 'R 0 1 2 3 4 5 6 7 8 9 10'
        assert reset_layer(5, 10) == 'R 5 6 7 8 9 10'

    def test_cnot_layer(self):
        assert cnot_layer(0, 10, 'down') == ["CNOT 0 1 2 3 4 5 6 7 8 9", "TICK"]
        assert cnot_layer(0, 10, 'up') == ["CNOT 10 9 8 7 6 5 4 3 2 1", "TICK"]
        assert cnot_layer(1, 10, 'down') == ["CNOT 1 2 3 4 5 6 7 8 9 10", "TICK"]

    def test_measure_layer(self):
        assert measure_layer(0, 10) == [
            'MR 1 3 5 7 9',
            'DETECTOR rec[-6] rec[-1]',
            'DETECTOR rec[-7] rec[-2]',
            'DETECTOR rec[-8] rec[-3]',
            'DETECTOR rec[-9] rec[-4]',
            'DETECTOR rec[-10] rec[-5]',
            'TICK',
        ]
        assert measure_layer(1, 11) == [
            'MR 2 4 6 8 10',
            'DETECTOR rec[-6] rec[-1]',
            'DETECTOR rec[-7] rec[-2]',
            'DETECTOR rec[-8] rec[-3]',
            'DETECTOR rec[-9] rec[-4]',
            'DETECTOR rec[-10] rec[-5]',
            'TICK',
        ]

        assert measure_layer(0, 10, lost_qubit=50) == [
            'H 50',
            'TICK',
            'MR 1 3 5 7 9',
            'M 50',
            'DETECTOR rec[-7] rec[-2]',
            'DETECTOR rec[-8] rec[-3]',
            'DETECTOR rec[-9] rec[-4]',
            'DETECTOR rec[-10] rec[-5]',
            'DETECTOR rec[-11] rec[-6]',
            'TICK',
        ]

        assert measure_layer(0, 6, lost_qubit_last_round=True) == [
            'MR 1 3 5',
            'DETECTOR rec[-5] rec[-1]',
            'DETECTOR rec[-6] rec[-2]',
            'DETECTOR rec[-7] rec[-3]',
            'TICK',
        ]

        assert measure_layer(0, 6, lost_qubit=10, lost_qubit_last_round=True) == [
            'H 10',
            'TICK',
            'MR 1 3 5',
            'M 10',
            'DETECTOR rec[-6] rec[-2]',
            'DETECTOR rec[-7] rec[-3]',
            'DETECTOR rec[-8] rec[-4]',
            'TICK',
        ]

        assert measure_layer(0, 6, first_round=True) == [
            'MR 1 3 5',
            'DETECTOR rec[-1]',
            'DETECTOR rec[-2]',
            'DETECTOR rec[-3]',
            'TICK',
        ]

    def test_last_measure_layer(self):
        assert last_measure_layer(0, 6) == [
            'M 0 2 4 6',
            'DETECTOR rec[-5] rec[-2] rec[-1]',
            'DETECTOR rec[-6] rec[-3] rec[-2]',
            'DETECTOR rec[-7] rec[-4] rec[-3]',
            'OBSERVABLE_INCLUDE(0) rec[-1]',
        ]
        assert last_measure_layer(2, 8) == [
            'M 2 4 6 8',
            'DETECTOR rec[-5] rec[-2] rec[-1]',
            'DETECTOR rec[-6] rec[-3] rec[-2]',
            'DETECTOR rec[-7] rec[-4] rec[-3]',
            'OBSERVABLE_INCLUDE(0) rec[-1]',
        ]
        assert last_measure_layer(1, 11, lost_qubit_last_round=True) == [
            'M 1 3 5 7 9 11',
            'DETECTOR rec[-8] rec[-2] rec[-1]',
            'DETECTOR rec[-9] rec[-3] rec[-2]',
            'DETECTOR rec[-10] rec[-4] rec[-3]',
            'DETECTOR rec[-11] rec[-5] rec[-4]',
            'DETECTOR rec[-12] rec[-6] rec[-5]',
            'OBSERVABLE_INCLUDE(0) rec[-1]',
        ]

    @pytest.mark.parametrize("distance", [3, 5, 11])
    @pytest.mark.parametrize(
        "code_func",
        [
            normal_rep_code,
            other_normal_rep_code,
            strange_rep_code,
            gliding_step_code,
            wiggling_step_code,
        ],
    )
    def test_circuit_distances(self, distance, code_func):
        circuit = code_func(distance=distance, rounds=5, errors=0.1)
        assert len(circuit.shortest_graphlike_error()) == distance
