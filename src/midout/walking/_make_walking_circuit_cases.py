"""This file contains the interface for use of the walking circuit code by the rest of the code base"""

from midout._circuit_case import CircuitCase
from midout.gen._layer_translate import to_z_basis_interaction_circuit

from midout.walking.circuit import Circuit
from midout.walking.util import Basis, DOWN_RIGHT, DOWN_LEFT, UP_LEFT


strategy_to_rounds = {
    'gliding': [DOWN_RIGHT, DOWN_RIGHT],
    'sliding': [DOWN_RIGHT, DOWN_LEFT],
    'wiggling': [DOWN_RIGHT, UP_LEFT],
}


def make_walking_code(distance, basis, rounds, strategy, gate):

    rounds_kernel = strategy_to_rounds[strategy]
    rounds_directions = (rounds_kernel * (rounds // 2))
    if len(rounds_directions) < rounds:
        rounds_directions += [None]
    c = Circuit(distance=distance, rounds=rounds_directions, basis=Basis[basis.lower()], duids=False)
    stim_circuit = c.build_stim_circuit(errors=None)

    if gate == "CX":
        expected_interactions = frozenset(['CX'])
    elif gate == "CZ":
        expected_interactions = frozenset(['CZ'])
        stim_circuit = to_z_basis_interaction_circuit(stim_circuit)
    else:
        raise ValueError(f"Unrecognised Gate Type: {gate}")

    return CircuitCase(
        circuit=stim_circuit,
        patches=[],
        expected_interactions=expected_interactions
    )
