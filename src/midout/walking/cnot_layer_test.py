import itertools

from midout.walking.cnot_layer import _make_cnot_layer


def test_cnot_layer():
    cl = _make_cnot_layer(
        qubits={0, 1, 1j, 1 + 1j, 0.5 + 0.5j},
        center=0.5 + 0.5j,
        direction=0.5 + 0.5j,
        aligned=True,
        flip_cnots=False,
    )
    assert cl == {(0.5 + 0.5j, 1 + 1j)}

    cl = _make_cnot_layer(
        qubits={0, 1, 1j, 1 + 1j, 0.5 + 0.5j},
        center=0.5 + 0.5j,
        direction=-0.5 + 0.5j,  # different direction
        aligned=True,
        flip_cnots=True,  # also flip cnot
    )
    assert cl == {(1j, 0.5 + 0.5j)}

    cl = _make_cnot_layer(
        qubits={0, 1, 1j, 1 + 1j, 0.5 + 0.5j},
        center=0.5 - 0.5j,  # move center so 0.5+0.5j is on the odd sublattice
        direction=0.5 + 0.5j,
        aligned=True,
        flip_cnots=False,
    )
    assert cl == {(1 + 1j, 0.5 + 0.5j)}  # odd sublattice should also have reversed cnot

    # make a qubit set like
    # D   D   D
    #   M   M
    # D   D   D
    #   M   M
    # D   D   D
    data_qubits = [x + 1j * y for x, y in itertools.product(range(3), repeat=2)]
    measure_qubits = [x + 0.5 + 1j * (y + 0.5) for x, y in itertools.product(range(2), repeat=2)]
    qubits = set(data_qubits + measure_qubits)

    cl = _make_cnot_layer(
        qubits=qubits, center=0.5 + 0.5j, direction=0.5 + 0.5j, aligned=True, flip_cnots=False
    )
    assert cl == {
        (0.5 + 0.5j, 1 + 1j),
        (1.5 + 0.5j, 2 + 1j)[::-1],
        (0.5 + 1.5j, 1 + 2j)[::-1],
        (1.5 + 1.5j, 2 + 2j),
    }

    cl = _make_cnot_layer(
        qubits=qubits,
        center=0.5 + 0.5j,
        direction=0.5 + 0.5j,
        aligned=False,  # anti-aligned
        flip_cnots=False,
    )
    assert cl == {
        (0.5 + 0.5j, 1 + 1j),
        (1.5 + 0.5j, 1 + 0j)[::-1],
        (0.5 + 1.5j, 0 + 1j)[::-1],
        (1.5 + 1.5j, 2 + 2j),
    }

    cl = _make_cnot_layer(
        qubits=qubits,
        center=1.5 + 0.5j,  # change the sublattice
        direction=0.5 + 0.5j,
        aligned=True,
        flip_cnots=False,
    )
    assert cl == {
        (0.5 + 0.5j, 1 + 1j)[::-1],
        (1.5 + 0.5j, 2 + 1j),
        (0.5 + 1.5j, 1 + 2j),
        (1.5 + 1.5j, 2 + 2j)[::-1],
    }
