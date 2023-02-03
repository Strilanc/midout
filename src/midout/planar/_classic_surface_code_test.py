import stim

from midout.planar._classic_surface_code import make_square_planar_cz_code, make_square_planar_cx_code


def test_make_square_planar_cx_code():
    circuit = make_square_planar_cx_code(
        basis='X',
        distance=3,
        rounds=100,
    ).circuit
    assert circuit == stim.Circuit("""
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
        QUBIT_COORDS(0.5, 0.5) 10
        QUBIT_COORDS(0.5, 1.5) 11
        QUBIT_COORDS(0.5, 2.5) 12
        QUBIT_COORDS(1.5, -0.5) 13
        QUBIT_COORDS(1.5, 0.5) 14
        QUBIT_COORDS(1.5, 1.5) 15
        QUBIT_COORDS(2.5, 1.5) 16
        R 10 12 13 15
        RX 0 1 2 3 4 5 6 7 8 9 11 14 16
        TICK
        CX 1 10 3 13 5 15 11 2 14 4 16 8
        TICK
        CX 4 10 6 13 8 15 11 1 14 3 16 7
        TICK
        CX 0 10 2 12 4 15 9 1 11 5 14 7
        TICK
        CX 3 10 5 12 7 15 9 0 11 4 14 6
        TICK
        MX 9 11 14 16
        M 10 12 13 15
        DETECTOR(-0.5, 0.5, 0) rec[-8]
        DETECTOR(0.5, 1.5, 0) rec[-7]
        DETECTOR(1.5, 0.5, 0) rec[-6]
        DETECTOR(2.5, 1.5, 0) rec[-5]
        SHIFT_COORDS(0, 0, 1)
        TICK
        REPEAT 99 {
            R 10 12 13 15
            RX 9 11 14 16
            TICK
            CX 1 10 3 13 5 15 11 2 14 4 16 8
            TICK
            CX 4 10 6 13 8 15 11 1 14 3 16 7
            TICK
            CX 0 10 2 12 4 15 9 1 11 5 14 7
            TICK
            CX 3 10 5 12 7 15 9 0 11 4 14 6
            TICK
            MX 9 11 14 16
            M 10 12 13 15
            DETECTOR(-0.5, 0.5, 0) rec[-16] rec[-8]
            DETECTOR(0.5, 0.5, 0) rec[-12] rec[-4]
            DETECTOR(0.5, 1.5, 0) rec[-15] rec[-7]
            DETECTOR(0.5, 2.5, 0) rec[-11] rec[-3]
            DETECTOR(1.5, -0.5, 0) rec[-10] rec[-2]
            DETECTOR(1.5, 0.5, 0) rec[-14] rec[-6]
            DETECTOR(1.5, 1.5, 0) rec[-9] rec[-1]
            DETECTOR(2.5, 1.5, 0) rec[-13] rec[-5]
            SHIFT_COORDS(0, 0, 1)
            TICK
        }
        MX 0 1 2 3 4 5 6 7 8
        DETECTOR(-0.5, 0.5, 0) rec[-17] rec[-9] rec[-8]
        DETECTOR(0.5, 1.5, 0) rec[-16] rec[-8] rec[-7] rec[-5] rec[-4]
        DETECTOR(1.5, 0.5, 0) rec[-15] rec[-6] rec[-5] rec[-3] rec[-2]
        DETECTOR(2.5, 1.5, 0) rec[-14] rec[-2] rec[-1]
        OBSERVABLE_INCLUDE(0) rec[-9] rec[-6] rec[-3]
    """)


def test_make_square_planar_cz_code():
    circuit = make_square_planar_cz_code(
        basis='X',
        distance=3,
        rounds=100,
    ).circuit
    assert circuit == stim.Circuit("""
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
        QUBIT_COORDS(0.5, 0.5) 10
        QUBIT_COORDS(0.5, 1.5) 11
        QUBIT_COORDS(0.5, 2.5) 12
        QUBIT_COORDS(1.5, -0.5) 13
        QUBIT_COORDS(1.5, 0.5) 14
        QUBIT_COORDS(1.5, 1.5) 15
        QUBIT_COORDS(2.5, 1.5) 16
        R 10 12 13 15 0 1 2 3 4 5 6 7 8 9 11 14 16
        TICK
        H 0 1 3 5 6 9 10 11 12 13 14 15 16
        TICK
        CZ 1 10 2 11 3 13 4 14 5 15 8 16
        TICK
        H 1 2 3 4 5 8
        TICK
        CZ 1 11 3 14 4 10 6 13 7 16 8 15
        TICK
        CZ 0 10 1 9 2 12 4 15 5 11 7 14
        TICK
        H 0 1 3 4 5 6 7 13
        TICK
        CZ 0 9 3 10 4 11 5 12 6 14 7 15
        TICK
        H 0 4 6 9 10 11 12 14 15 16
        TICK
        M 9 11 14 16 10 12 13 15
        DETECTOR(-0.5, 0.5, 0) rec[-8]
        DETECTOR(0.5, 1.5, 0) rec[-7]
        DETECTOR(1.5, 0.5, 0) rec[-6]
        DETECTOR(2.5, 1.5, 0) rec[-5]
        SHIFT_COORDS(0, 0, 1)
        TICK
        REPEAT 99 {
            R 10 12 13 15 9 11 14 16
            TICK
            H 2 4 8 9 10 11 12 13 14 15 16
            TICK
            CZ 1 10 2 11 3 13 4 14 5 15 8 16
            TICK
            H 1 2 3 4 5 7 8
            TICK
            CZ 1 11 3 14 4 10 6 13 7 16 8 15
            TICK
            CZ 0 10 1 9 2 12 4 15 5 11 7 14
            TICK
            H 0 1 3 4 5 6 7 13
            TICK
            CZ 0 9 3 10 4 11 5 12 6 14 7 15
            TICK
            H 0 4 6 9 10 11 12 14 15 16
            TICK
            M 9 11 14 16 10 12 13 15
            DETECTOR(-0.5, 0.5, 0) rec[-16] rec[-8]
            DETECTOR(0.5, 0.5, 0) rec[-12] rec[-4]
            DETECTOR(0.5, 1.5, 0) rec[-15] rec[-7]
            DETECTOR(0.5, 2.5, 0) rec[-11] rec[-3]
            DETECTOR(1.5, -0.5, 0) rec[-10] rec[-2]
            DETECTOR(1.5, 0.5, 0) rec[-14] rec[-6]
            DETECTOR(1.5, 1.5, 0) rec[-9] rec[-1]
            DETECTOR(2.5, 1.5, 0) rec[-13] rec[-5]
            SHIFT_COORDS(0, 0, 1)
            TICK
        }
        H 0 1 2 3 4 5 6 7 8 9 11 14 16
        TICK
        M 0 1 2 3 4 5 6 7 8
        DETECTOR(-0.5, 0.5, 0) rec[-17] rec[-9] rec[-8]
        DETECTOR(0.5, 1.5, 0) rec[-16] rec[-8] rec[-7] rec[-5] rec[-4]
        DETECTOR(1.5, 0.5, 0) rec[-15] rec[-6] rec[-5] rec[-3] rec[-2]
        DETECTOR(2.5, 1.5, 0) rec[-14] rec[-2] rec[-1]
        OBSERVABLE_INCLUDE(0) rec[-9] rec[-6] rec[-3]
    """)
