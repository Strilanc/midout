from midout import gen


def checkerboard_basis(q: complex) -> str:
    """Classifies a coordinate as X type or Z type according to a checkerboard.
    """
    is_x = int(q.real + q.imag) & 1 == 0
    return 'X' if is_x else 'Z'


def surface_code_patch(*, distance: int) -> gen.Patch:
    top_bot_basis = 'Z'
    left_right_basis = 'X'

    ur, ul, dl, dr = [(0.5 + 0.5j)*1j**d for d in range(4)]
    order_z = [ul, ur, dl, dr]
    order_ᴎ = [ul, dl, ur, dr]

    data_qubits = {
        x + 1j*y
        for x in range(distance)
        for y in range(distance)
    }
    potential_measure_qubits = {
        q + d
        for q in data_qubits
        for d in order_z
    }

    tiles = []
    for m in gen.sorted_complex(potential_measure_qubits):
        on_top = m.imag < 0
        on_bottom = m.imag > distance - 1
        on_left = m.real < 0
        on_right = m.real > distance - 1
        tile_basis = checkerboard_basis(m)
        if (on_top or on_bottom) and tile_basis != top_bot_basis:
            continue
        if (on_left or on_right) and tile_basis != left_right_basis:
            continue
        order = order_z if tile_basis == 'Z' else order_ᴎ
        tiles.append(gen.Tile(
            bases=tile_basis,
            measurement_qubit=m,
            ordered_data_qubits=[m + d if m + d in data_qubits else None for d in order],
        ))

    return gen.Patch(tiles)
