import pytest
import stim

from midout.walking.circuit import Circuit
from midout.walking.util import Basis, DOWN_RIGHT, DOWN_LEFT, UP_LEFT


@pytest.mark.parametrize("distance", [3, 4, 5, 6, 7, 15])
@pytest.mark.parametrize("rounds", ['d', 4, 5, 6])
def test_regular_surface_code_distances(distance, rounds):
    if rounds is 'd':
        rounds = distance
    c = Circuit(distance=distance, rounds=[None] * rounds, basis=Basis.x)
    stim_circuit = c.build_stim_circuit(errors=1e-3)
    error_path = stim_circuit.shortest_graphlike_error()
    assert len(error_path) == distance


strategy_to_rounds = {
    'gliding': [DOWN_RIGHT] * 2,
    'sliding': [DOWN_RIGHT, DOWN_LEFT],
    'wiggling': [DOWN_RIGHT, UP_LEFT],
}


@pytest.mark.parametrize("distance", [3, 4, 5, 6, 7, 15])
@pytest.mark.parametrize("strategy", list(strategy_to_rounds.keys()))
@pytest.mark.parametrize("basis", [Basis.z, Basis.x])
@pytest.mark.parametrize("num_rounds", ['d', 4, 5, 6])
def test_walking_surface_code_distances(distance, strategy, basis, num_rounds):
    if num_rounds is 'd':
        num_rounds = distance
    rounds = strategy_to_rounds[strategy] * (num_rounds//2)
    if num_rounds%2:
        rounds += [None]
    assert len(rounds) == num_rounds
    c = Circuit(distance=distance, rounds=rounds, basis=basis, duids=False)
    stim_circuit = c.build_stim_circuit(errors=1e-3)
    error_path = stim_circuit.shortest_graphlike_error()
    assert len(error_path) == distance


def test_circuit_determinism():
    distance = 10
    rounds = [DOWN_RIGHT] * 6
    c0 = Circuit(distance=distance, rounds=rounds, basis=Basis.x)
    c1 = Circuit(distance=distance, rounds=rounds, basis=Basis.x)

    stim_circuit_0 = c0.build_stim_circuit()
    stim_circuit_1 = c1.build_stim_circuit()

    assert stim_circuit_0 == stim_circuit_1
