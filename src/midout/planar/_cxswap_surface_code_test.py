import itertools

import pytest
import stim

from midout import gen
from midout.planar._cxswap_surface_code import make_cx_swap_surface_code_patch, \
    make_cx_swap_surface_code_chunk, make_square_planar_cxswap_code


def test_make_cx_swap_surface_code_patch():
    patch = make_cx_swap_surface_code_patch(distance=2)
    assert patch == gen.Patch(tiles=[
        gen.Tile(
            ordered_data_qubits=(0j, (0.5+0.5j), 1j, (-0.5+0.5j)),
            measurement_qubit=(-0.5+0.5j),
            bases='X',
        ),
        gen.Tile(
            ordered_data_qubits=((-0.5+0.5j),),
            measurement_qubit=(-0.5+0.5j),
            bases='X',
        ),
        gen.Tile(
            ordered_data_qubits=((0.5-0.5j),),
            measurement_qubit=(0.5-0.5j),
            bases='Z',
        ),
        gen.Tile(
            ordered_data_qubits=((0.5-0.5j), (1+0j), (0.5+0.5j), 0j),
            measurement_qubit=(0.5+0.5j),
            bases='Z',
        ),
        gen.Tile(
            ordered_data_qubits=((0.5+0.5j), (1+1j), (0.5+1.5j), 1j),
            measurement_qubit=(0.5+0.5j),
            bases='Z',
        ),
        gen.Tile(
            ordered_data_qubits=((0.5+1.5j),),
            measurement_qubit=(0.5+1.5j),
            bases='Z',
        ),
        gen.Tile(
            ordered_data_qubits=((1+0j), (1.5+0.5j), (1+1j), (0.5+0.5j)),
            measurement_qubit=(1.5+0.5j),
            bases='X',
        ),
        gen.Tile(
            ordered_data_qubits=((1.5+0.5j),),
            measurement_qubit=(1.5+0.5j),
            bases='X',
        ),
    ])


@pytest.mark.parametrize('distance,basis,round_parity,is_first_round', itertools.product(
    [2, 5, 6],
    'XZ',
    [False, True],
    [False, True]
))
def test_make_cx_swap_surface_code_chunk(distance: int, basis: str, round_parity: bool, is_first_round: bool):
    r1 = make_cx_swap_surface_code_chunk(
        distance=distance,
        round_parity=round_parity,
        basis=basis,
        is_first_round=is_first_round,
    )
    r1.verify()

    r2 = make_cx_swap_surface_code_chunk(
        distance=distance,
        round_parity=not round_parity,
        basis=basis,
        is_first_round=False,
    )
    r2.verify()

    gen.compile_chunks_into_circuit([r1.magic_init_chunk(), r1, r2, r2.magic_end_chunk()])


def test_make_square_planar_cxswap_code():
    assert make_square_planar_cxswap_code(
        distance=3,
        basis='X',
        rounds=100,
    ).circuit == stim.Circuit('''
        QUBIT_COORDS(0, 0) 0
        QUBIT_COORDS(0, 1) 1
        QUBIT_COORDS(0, 2) 2
        QUBIT_COORDS(1, 0) 3
        QUBIT_COORDS(1, 1) 4
        QUBIT_COORDS(1, 2) 5
        QUBIT_COORDS(2, 0) 6
        QUBIT_COORDS(2, 1) 7
        QUBIT_COORDS(2, 2) 8
        QUBIT_COORDS(-0.5, 0.5) 9
        QUBIT_COORDS(-0.5, 1.5) 10
        QUBIT_COORDS(0.5, -0.5) 11
        QUBIT_COORDS(0.5, 0.5) 12
        QUBIT_COORDS(0.5, 1.5) 13
        QUBIT_COORDS(0.5, 2.5) 14
        QUBIT_COORDS(1.5, -0.5) 15
        QUBIT_COORDS(1.5, 0.5) 16
        QUBIT_COORDS(1.5, 1.5) 17
        QUBIT_COORDS(1.5, 2.5) 18
        QUBIT_COORDS(2.5, 0.5) 19
        QUBIT_COORDS(2.5, 1.5) 20
        RX 9 10 13 16 20 0 1 2 3 4 5 6 7 8 11 19
        R 12 14 15 17 18
        TICK
        CXSWAP 9 0 13 4 16 6 3 12 5 14 7 17
        TICK
        CXSWAP 9 1 11 3 13 5 16 7 0 12 2 14 4 17 6 19
        TICK
        CXSWAP 1 13 3 16 5 18 7 20 10 2 12 4 15 6 17 8
        TICK
        CXSWAP 2 13 4 16 8 20 12 1 15 3 17 5
        TICK
        MX 9 13 16 19 20
        M 11 12 14 15 17
        TICK
        RX 9 13 16 19 20
        R 11 12 14 15 17
        TICK
        CXSWAP 13 2 16 4 20 8 1 12 3 15 5 17
        TICK
        CXSWAP 13 1 16 3 18 5 20 7 2 10 4 12 6 15 8 17
        DETECTOR(-0.5, 0.5, 0) rec[-10]
        DETECTOR(0.5, 1.5, 0) rec[-9]
        DETECTOR(1.5, 0.5, 0) rec[-8]
        DETECTOR(2.5, 0.5, 0) rec[-7]
        DETECTOR(2.5, 1.5, 0) rec[-6]
        SHIFT_COORDS(0, 0, 1)
        TICK
        CXSWAP 1 9 3 11 5 13 7 16 12 0 14 2 17 4 19 6
        TICK
        CXSWAP 0 9 4 13 6 16 12 3 14 5 17 7
        TICK
        MX 9 10 13 16 20
        M 12 14 15 17 18
        TICK
        RX 9 10 13 16 20
        R 12 14 15 17 18
        TICK
        CXSWAP 9 0 13 4 16 6 3 12 5 14 7 17
        TICK
        CXSWAP 9 1 11 3 13 5 16 7 0 12 2 14 4 17 6 19
        DETECTOR(-0.5, 0.5, 0) rec[-20] rec[-18] rec[-10]
        DETECTOR(-0.5, 1.5, 0) rec[-19] rec[-9]
        DETECTOR(0.5, 0.5, 0) rec[-15] rec[-5]
        DETECTOR(0.5, 1.5, 0) rec[-16] rec[-8]
        DETECTOR(0.5, 2.5, 0) rec[-14] rec[-13] rec[-4]
        DETECTOR(1.5, -0.5, 0) rec[-3]
        DETECTOR(1.5, 0.5, 0) rec[-17] rec[-7]
        DETECTOR(1.5, 1.5, 0) rec[-12] rec[-2]
        DETECTOR(1.5, 2.5, 0) rec[-11] rec[-1]
        DETECTOR(2.5, 1.5, 0) rec[-6]
        SHIFT_COORDS(0, 0, 1)
        TICK
        CXSWAP 1 13 3 16 5 18 7 20 10 2 12 4 15 6 17 8
        TICK
        CXSWAP 2 13 4 16 8 20 12 1 15 3 17 5
        TICK
        MX 9 13 16 19 20
        M 11 12 14 15 17
        TICK
        RX 9 13 16 19 20
        R 11 12 14 15 17
        TICK
        CXSWAP 13 2 16 4 20 8 1 12 3 15 5 17
        TICK
        CXSWAP 13 1 16 3 18 5 20 7 2 10 4 12 6 15 8 17
        DETECTOR(-0.5, 0.5, 0) rec[-10]
        DETECTOR(0.5, -0.5, 0) rec[-15] rec[-5]
        DETECTOR(0.5, 0.5, 0) rec[-14] rec[-4]
        DETECTOR(0.5, 1.5, 0) rec[-19] rec[-9]
        DETECTOR(0.5, 2.5, 0) rec[-3]
        DETECTOR(1.5, -0.5, 0) rec[-13] rec[-12] rec[-2]
        DETECTOR(1.5, 0.5, 0) rec[-20] rec[-8]
        DETECTOR(1.5, 1.5, 0) rec[-11] rec[-1]
        DETECTOR(2.5, 0.5, 0) rec[-17] rec[-7]
        DETECTOR(2.5, 1.5, 0) rec[-18] rec[-16] rec[-6]
        SHIFT_COORDS(0, 0, 1)
        TICK
        REPEAT 48 {
            CXSWAP 1 9 3 11 5 13 7 16 12 0 14 2 17 4 19 6
            TICK
            CXSWAP 0 9 4 13 6 16 12 3 14 5 17 7
            TICK
            MX 9 10 13 16 20
            M 12 14 15 17 18
            TICK
            RX 9 10 13 16 20
            R 12 14 15 17 18
            TICK
            CXSWAP 9 0 13 4 16 6 3 12 5 14 7 17
            TICK
            CXSWAP 9 1 11 3 13 5 16 7 0 12 2 14 4 17 6 19
            DETECTOR(-0.5, 0.5, 0) rec[-20] rec[-18] rec[-10]
            DETECTOR(-0.5, 1.5, 0) rec[-19] rec[-9]
            DETECTOR(0.5, 0.5, 0) rec[-15] rec[-5]
            DETECTOR(0.5, 1.5, 0) rec[-16] rec[-8]
            DETECTOR(0.5, 2.5, 0) rec[-14] rec[-13] rec[-4]
            DETECTOR(1.5, -0.5, 0) rec[-3]
            DETECTOR(1.5, 0.5, 0) rec[-17] rec[-7]
            DETECTOR(1.5, 1.5, 0) rec[-12] rec[-2]
            DETECTOR(1.5, 2.5, 0) rec[-11] rec[-1]
            DETECTOR(2.5, 1.5, 0) rec[-6]
            SHIFT_COORDS(0, 0, 1)
            TICK
            CXSWAP 1 13 3 16 5 18 7 20 10 2 12 4 15 6 17 8
            TICK
            CXSWAP 2 13 4 16 8 20 12 1 15 3 17 5
            TICK
            MX 9 13 16 19 20
            M 11 12 14 15 17
            TICK
            RX 9 13 16 19 20
            R 11 12 14 15 17
            TICK
            CXSWAP 13 2 16 4 20 8 1 12 3 15 5 17
            TICK
            CXSWAP 13 1 16 3 18 5 20 7 2 10 4 12 6 15 8 17
            DETECTOR(-0.5, 0.5, 0) rec[-10]
            DETECTOR(0.5, -0.5, 0) rec[-15] rec[-5]
            DETECTOR(0.5, 0.5, 0) rec[-14] rec[-4]
            DETECTOR(0.5, 1.5, 0) rec[-19] rec[-9]
            DETECTOR(0.5, 2.5, 0) rec[-3]
            DETECTOR(1.5, -0.5, 0) rec[-13] rec[-12] rec[-2]
            DETECTOR(1.5, 0.5, 0) rec[-20] rec[-8]
            DETECTOR(1.5, 1.5, 0) rec[-11] rec[-1]
            DETECTOR(2.5, 0.5, 0) rec[-17] rec[-7]
            DETECTOR(2.5, 1.5, 0) rec[-18] rec[-16] rec[-6]
            SHIFT_COORDS(0, 0, 1)
            TICK
        }
        CXSWAP 17 8 15 6 12 4 10 2 7 20 5 18 3 16 1 13
        TICK
        CXSWAP 17 5 15 3 12 1 8 20 4 16 2 13
        TICK
        M 17 15 14 12 11
        MX 18 10 8 7 6 5 4 3 2 1 0 20 19 16 13 9
        DETECTOR(-0.5, 0.5, 0) rec[-31] rec[-29] rec[-10] rec[-8] rec[-7] rec[-6] rec[-3] rec[-1]
        DETECTOR(-0.5, 0.5, 0) rec[-1]
        DETECTOR(-0.5, 1.5, 0) rec[-30] rec[-15] rec[-8] rec[-2]
        DETECTOR(0.5, -0.5, 0) rec[-17]
        DETECTOR(0.5, 0.5, 0) rec[-18]
        DETECTOR(0.5, 1.5, 0) rec[-2]
        DETECTOR(0.5, 1.5, 0) rec[-27] rec[-16] rec[-14] rec[-11] rec[-7] rec[-5]
        DETECTOR(0.5, 2.5, 0) rec[-19]
        DETECTOR(1.5, -0.5, 0) rec[-20]
        DETECTOR(1.5, 0.5, 0) rec[-3]
        DETECTOR(1.5, 0.5, 0) rec[-28] rec[-14] rec[-13] rec[-9] rec[-4]
        DETECTOR(1.5, 1.5, 0) rec[-21]
        DETECTOR(2.5, 0.5, 0) rec[-4]
        DETECTOR(2.5, 1.5, 0) rec[-5]
        DETECTOR(2.5, 1.5, 0) rec[-13]
        OBSERVABLE_INCLUDE(0) rec[-10] rec[-9] rec[-6]
        SHIFT_COORDS(0, 0, 1)
        TICK
    ''')
