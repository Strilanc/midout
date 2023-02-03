import stim

from midout.planar._iswap_surface_code import make_square_planar_iswap_code


def test_make_square_planar_iswap_code():
    circuit = make_square_planar_iswap_code(
        basis='X',
        distance=3,
        rounds=100,
    ).circuit
    assert circuit == stim.Circuit("""
        QUBIT_COORDS(0, 0) 0
        QUBIT_COORDS(0, 1) 1
        QUBIT_COORDS(0, 2) 2
        QUBIT_COORDS(1, -1) 3
        QUBIT_COORDS(1, 0) 4
        QUBIT_COORDS(1, 1) 5
        QUBIT_COORDS(1, 2) 6
        QUBIT_COORDS(2, -1) 7
        QUBIT_COORDS(2, 0) 8
        QUBIT_COORDS(2, 1) 9
        QUBIT_COORDS(2, 2) 10
        QUBIT_COORDS(3, -1) 11
        QUBIT_COORDS(3, 0) 12
        QUBIT_COORDS(3, 1) 13
        QUBIT_COORDS(0.5, -0.5) 14
        QUBIT_COORDS(0.5, 0.5) 15
        QUBIT_COORDS(0.5, 1.5) 16
        QUBIT_COORDS(1.5, -0.5) 17
        QUBIT_COORDS(1.5, 0.5) 18
        QUBIT_COORDS(1.5, 1.5) 19
        QUBIT_COORDS(2.5, -0.5) 20
        QUBIT_COORDS(2.5, 0.5) 21
        QUBIT_COORDS(2.5, 1.5) 22
        R 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
        TICK
        H 11 13 17 19 21
        H_YZ 1 3 4 5 8 9
        TICK
        ISWAP 1 16 3 17 4 18 5 19 8 21 9 22
        TICK
        C_XYZ 1 3 4 5 8 9
        H 6 11 13
        H_YZ 0 7
        TICK
        ISWAP 4 17 5 16 6 19 8 18 9 21 13 22 3 14 1 15 11 20
        TICK
        ISWAP 2 16 10 22 12 21 0 14 4 15 7 17 5 18 9 19 8 20
        TICK
        C_ZYX 4 5 8 9 10 12
        H 0 2 7
        S 14 15 16 17 18 19 20 21 22
        TICK
        ISWAP 4 14 5 15 8 17 9 18 10 19 12 20
        TICK
        H 4 5 7 8 9 10 12 14 18 20
        H_YZ 6 13
        TICK
        M 1 3 4 5 6 8 9 10 11 12 13 16 21 22
        DETECTOR(0.5, -0.5, 0) rec[-12]
        DETECTOR(1.5, 0.5, 0) rec[-8]
        DETECTOR(2.5, -0.5, 0) rec[-5]
        DETECTOR(3.5, 0.5, 0) rec[-4]
        DETECTOR(0, 1, 0) rec[-14]
        DETECTOR(1, -1, 0) rec[-13]
        DETECTOR(3, -1, 0) rec[-6]
        DETECTOR(0.5, 1.5, 0) rec[-3]
        DETECTOR(2.5, 0.5, 0) rec[-2]
        DETECTOR(2.5, 1.5, 0) rec[-1]
        SHIFT_COORDS(0, 0, 1)
        TICK
        R 1 3 4 5 6 8 9 10 11 12 13 16 21 22
        TICK
        H 4 5 7 8 9 10 12 14 18 20
        TICK
        ISWAP 4 14 5 15 8 17 9 18 10 19 12 20
        TICK
        C_ZYX 0 2 7
        H 6 13
        S 16 21 22
        H_YZ 4 5 8 9 10 12
        TICK
        ISWAP 2 16 10 22 12 21 0 14 4 15 7 17 5 18 9 19 8 20
        TICK
        ISWAP 4 17 5 16 6 19 8 18 9 21 13 22 3 14 1 15 11 20
        TICK
        C_XYZ 1 3 4 5 6 8 9 11 13
        S 16 17 18 19 21 22
        TICK
        ISWAP 1 16 3 17 4 18 5 19 8 21 9 22
        TICK
        C_XYZ 1 3 4 5 8 9 16 18 22
        H 0 6 7
        S 17 19 21
        TICK
        M 0 1 2 3 4 5 7 8 9 10 12 14 15 20
        DETECTOR(-0.5, 0.5, 0) rec[-26] rec[-14]
        DETECTOR(0.5, 1.5, 0) rec[-22] rec[-13]
        DETECTOR(1.5, -0.5, 0) rec[-25] rec[-11]
        DETECTOR(1.5, 0.5, 0) rec[-19] rec[-10]
        DETECTOR(1.5, 1.5, 0) rec[-24] rec[-9]
        DETECTOR(2.5, -1.5, 0) rec[-23] rec[-8]
        DETECTOR(2.5, 0.5, 0) rec[-21] rec[-7]
        DETECTOR(2.5, 1.5, 0) rec[-18] rec[-6]
        DETECTOR(0, 2, 0) rec[-12]
        DETECTOR(2, 2, 0) rec[-5]
        DETECTOR(3, 0, 0) rec[-4]
        DETECTOR(0.5, -0.5, 0) rec[-3]
        DETECTOR(0.5, 0.5, 0) rec[-2]
        DETECTOR(2.5, -0.5, 0) rec[-1]
        TICK
        REPEAT 49 {
            R 0 1 2 3 4 5 7 8 9 10 12 14 15 20
            TICK
            H 6 16 18 22
            H_YZ 1 3 4 5 8 9
            TICK
            ISWAP 1 16 3 17 4 18 5 19 8 21 9 22
            TICK
            C_XYZ 1 3 4 5 8 9
            H 6 11 13
            H_YZ 0 7
            TICK
            ISWAP 4 17 5 16 6 19 8 18 9 21 13 22 3 14 1 15 11 20
            TICK
            ISWAP 2 16 10 22 12 21 0 14 4 15 7 17 5 18 9 19 8 20
            TICK
            C_ZYX 4 5 8 9 10 12
            H 0 2 7
            S 14 15 16 17 18 19 20 21 22
            TICK
            ISWAP 4 14 5 15 8 17 9 18 10 19 12 20
            TICK
            H 4 5 7 8 9 10 12 14 18 20
            H_YZ 6 13
            TICK
            M 1 3 4 5 6 8 9 10 11 12 13 16 21 22
            DETECTOR(0.5, -0.5, 0) rec[-28] rec[-12]
            DETECTOR(0.5, 0.5, 0) rec[-25] rec[-11]
            DETECTOR(0.5, 2.5, 0) rec[-23] rec[-10]
            DETECTOR(1.5, -0.5, 0) rec[-22] rec[-9]
            DETECTOR(1.5, 0.5, 0) rec[-27] rec[-8]
            DETECTOR(1.5, 1.5, 0) rec[-21] rec[-7]
            DETECTOR(2.5, -0.5, 0) rec[-24] rec[-5]
            DETECTOR(3.5, 0.5, 0) rec[-20] rec[-4]
            DETECTOR(0, 1, 0) rec[-14]
            DETECTOR(1, -1, 0) rec[-13]
            DETECTOR(3, -1, 0) rec[-6]
            DETECTOR(0.5, 1.5, 0) rec[-3]
            DETECTOR(2.5, 0.5, 0) rec[-2]
            DETECTOR(2.5, 1.5, 0) rec[-1]
            SHIFT_COORDS(0, 0, 1)
            TICK
            R 1 3 4 5 6 8 9 10 11 12 13 16 21 22
            TICK
            H 4 5 7 8 9 10 12 14 18 20
            TICK
            ISWAP 4 14 5 15 8 17 9 18 10 19 12 20
            TICK
            C_ZYX 0 2 7
            H 6 13
            S 16 21 22
            H_YZ 4 5 8 9 10 12
            TICK
            ISWAP 2 16 10 22 12 21 0 14 4 15 7 17 5 18 9 19 8 20
            TICK
            ISWAP 4 17 5 16 6 19 8 18 9 21 13 22 3 14 1 15 11 20
            TICK
            C_XYZ 1 3 4 5 6 8 9 11 13
            S 16 17 18 19 21 22
            TICK
            ISWAP 1 16 3 17 4 18 5 19 8 21 9 22
            TICK
            C_XYZ 1 3 4 5 8 9 16 18 22
            H 0 6 7
            S 17 19 21
            TICK
            M 0 1 2 3 4 5 7 8 9 10 12 14 15 20
            DETECTOR(-0.5, 0.5, 0) rec[-26] rec[-14]
            DETECTOR(0.5, 1.5, 0) rec[-22] rec[-13]
            DETECTOR(1.5, -0.5, 0) rec[-25] rec[-11]
            DETECTOR(1.5, 0.5, 0) rec[-19] rec[-10]
            DETECTOR(1.5, 1.5, 0) rec[-24] rec[-9]
            DETECTOR(2.5, -1.5, 0) rec[-23] rec[-8]
            DETECTOR(2.5, 0.5, 0) rec[-21] rec[-7]
            DETECTOR(2.5, 1.5, 0) rec[-18] rec[-6]
            DETECTOR(0, 2, 0) rec[-12]
            DETECTOR(2, 2, 0) rec[-5]
            DETECTOR(3, 0, 0) rec[-4]
            DETECTOR(0.5, -0.5, 0) rec[-3]
            DETECTOR(0.5, 0.5, 0) rec[-2]
            DETECTOR(2.5, -0.5, 0) rec[-1]
            TICK
        }
        H 6 11 13 16 17 18 19 21 22
        TICK
        M 6 11 13 16 17 18 19 21 22
        OBSERVABLE_INCLUDE(0) rec[-8] rec[-5] rec[-4]
        DETECTOR(-0.5, 0.5, 0) rec[-23] rec[-6] rec[-5]
        DETECTOR(0.5, 1.5, 0) rec[-22] rec[-9] rec[-6] rec[-3] rec[-2]
        DETECTOR(1.5, 0.5, 0) rec[-19] rec[-8] rec[-4] rec[-2] rec[-1]
        DETECTOR(2.5, 1.5, 0) rec[-15] rec[-7] rec[-1]
        SHIFT_COORDS(0, 0, 1)
        TICK
    """)
