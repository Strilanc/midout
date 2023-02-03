import stim

from midout.planar._three_coupler_surface_codes import make_hex_planar_surface_code


def test_cz_circuit_details():
    circuit = make_hex_planar_surface_code(
        basis='X',
        distance=3,
        gate='CZ',
        wiggle=False,
        rounds=100,
    ).circuit
    assert circuit == stim.Circuit("""
        QUBIT_COORDS(0, 0) 0
        QUBIT_COORDS(0, 1) 1
        QUBIT_COORDS(0, 2) 2
        QUBIT_COORDS(0.5, 0.5) 3
        QUBIT_COORDS(0.5, 1.5) 4
        QUBIT_COORDS(0.5, 2.5) 5
        QUBIT_COORDS(1, 0) 6
        QUBIT_COORDS(1, 1) 7
        QUBIT_COORDS(1, 2) 8
        QUBIT_COORDS(1.5, 0.5) 9
        QUBIT_COORDS(1.5, 1.5) 10
        QUBIT_COORDS(1.5, 2.5) 11
        QUBIT_COORDS(2, 0) 12
        QUBIT_COORDS(2, 1) 13
        QUBIT_COORDS(2, 2) 14
        QUBIT_COORDS(2.5, 0.5) 15
        QUBIT_COORDS(2.5, 1.5) 16
        R 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16
        TICK
        H 1 4 6 9 13 14 15 16
        TICK
        M 3 5 10 11 4 9 15 16
        DETECTOR(0.5, 0.5, 0) rec[-8]
        DETECTOR(0.5, 2.5, 0) rec[-7]
        DETECTOR(1.5, 1.5, 0) rec[-6]
        DETECTOR(1.5, 2.5, 0) rec[-5]
        SHIFT_COORDS(0, 0, 1)
        TICK
        R 4 9 15 16 3 5 10 11
        TICK
        H 3 4 5 9 10 11 15 16
        TICK
        CZ 0 3 1 4 2 5 6 9 7 10 13 16
        TICK
        H 0 2 6 7 13
        TICK
        CZ 2 4 3 6 5 8 7 9 10 13 14 16
        TICK
        H 3 4 5 7 8 9 10 13 14 16
        TICK
        CZ 1 3 4 7 8 10 9 12 11 14 13 15
        TICK
        H 1 7 8 12 14
        TICK
        CZ 0 3 1 4 6 9 7 10 8 11 12 15
        TICK
        H 1 3 4 6 8 9 10 11 15
        TICK
        M 4 5 9 11 3 10 15 16
        DETECTOR(0.5, 1.5, 0) rec[-16] rec[-8]
        DETECTOR(0.5, 2.5, 0) rec[-15] rec[-7]
        DETECTOR(1.5, 0.5, 0) rec[-6]
        DETECTOR(1.5, 2.5, 0) rec[-14] rec[-13] rec[-5]
        SHIFT_COORDS(0, 0, 1)
        TICK
        R 3 10 15 16 4 5 9 11
        TICK
        REPEAT 49 {
            H 1 3 4 5 6 8 9 10 11 15 16
            TICK
            CZ 0 3 1 4 6 9 7 10 8 11 12 15
            TICK
            H 0 1 7 8 12 14
            TICK
            CZ 1 3 4 7 8 10 9 12 11 14 13 15
            TICK
            H 3 4 7 8 9 10 11 13 14 15
            TICK
            CZ 2 4 3 6 5 8 7 9 10 13 14 16
            TICK
            H 2 6 7 13
            TICK
            CZ 0 3 1 4 2 5 6 9 7 10 13 16
            TICK
            H 3 4 5 9 10 16
            TICK
            M 3 5 10 11 4 9 15 16
            DETECTOR(0.5, 0.5, 0) rec[-8]
            DETECTOR(0.5, 1.5, 0) rec[-4]
            DETECTOR(0.5, 2.5, 0) rec[-16] rec[-15] rec[-7]
            DETECTOR(1.5, 0.5, 0) rec[-12] rec[-3]
            DETECTOR(1.5, 1.5, 0) rec[-14] rec[-6]
            DETECTOR(1.5, 2.5, 0) rec[-13] rec[-5]
            DETECTOR(2.5, 0.5, 0) rec[-10] rec[-2]
            DETECTOR(2.5, 1.5, 0) rec[-11] rec[-9] rec[-1]
            SHIFT_COORDS(0, 0, 1)
            TICK
            R 4 9 15 16 3 5 10 11
            TICK
            H 3 4 5 9 10 11 15 16
            TICK
            CZ 0 3 1 4 2 5 6 9 7 10 13 16
            TICK
            H 0 2 6 7 13
            TICK
            CZ 2 4 3 6 5 8 7 9 10 13 14 16
            TICK
            H 3 4 5 7 8 9 10 13 14 16
            TICK
            CZ 1 3 4 7 8 10 9 12 11 14 13 15
            TICK
            H 1 7 8 12 14
            TICK
            CZ 0 3 1 4 6 9 7 10 8 11 12 15
            TICK
            H 1 3 4 6 8 9 10 11 15
            TICK
            M 4 5 9 11 3 10 15 16
            DETECTOR(0.5, 0.5, 0) rec[-4]
            DETECTOR(0.5, 1.5, 0) rec[-16] rec[-8]
            DETECTOR(0.5, 2.5, 0) rec[-15] rec[-7]
            DETECTOR(1.5, 0.5, 0) rec[-6]
            DETECTOR(1.5, 1.5, 0) rec[-12] rec[-3]
            DETECTOR(1.5, 2.5, 0) rec[-14] rec[-13] rec[-5]
            DETECTOR(2.5, 0.5, 0) rec[-11] rec[-10] rec[-2]
            DETECTOR(2.5, 1.5, 0) rec[-9] rec[-1]
            SHIFT_COORDS(0, 0, 1)
            TICK
            R 3 10 15 16 4 5 9 11
            TICK
        }
        H 1 3 4 6 8 9 10 11 15 16
        TICK
        CZ 0 3 1 4 6 9 7 10 8 11 12 15
        TICK
        H 0 1 7 8 12 14
        TICK
        CZ 1 3 4 7 8 10 9 12 11 14 13 15
        TICK
        H 1 2 4 8 9 11 13
        TICK
        M 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16
        DETECTOR(0.5, 0, 0) rec[-17] rec[-14] rec[-11]
        DETECTOR(0.5, 1, 0) rec[-16] rec[-14] rec[-13] rec[-10]
        DETECTOR(0.5, 2, 0) rec[-25] rec[-24] rec[-15] rec[-13] rec[-12] rec[-9]
        DETECTOR(0.5, 2.5, 0) rec[-12]
        DETECTOR(1.5, 0, 0) rec[-11] rec[-8] rec[-5]
        DETECTOR(1.5, 1, 0) rec[-23] rec[-10] rec[-8] rec[-7] rec[-4]
        DETECTOR(1.5, 2, 0) rec[-9] rec[-7] rec[-6] rec[-3]
        DETECTOR(1.5, 2.5, 0) rec[-22] rec[-6]
        OBSERVABLE_INCLUDE(0) rec[-17] rec[-16] rec[-15]
    """)


def test_iswap_circuit_details():
    circuit = make_hex_planar_surface_code(
        basis='X',
        distance=3,
        gate='ISWAP',
        wiggle=False,
        rounds=100,
    ).circuit
    assert circuit == stim.Circuit("""
        QUBIT_COORDS(-0.5, 0.5) 0
        QUBIT_COORDS(-0.5, 1.5) 1
        QUBIT_COORDS(0, 0) 2
        QUBIT_COORDS(0, 1) 3
        QUBIT_COORDS(0, 2) 4
        QUBIT_COORDS(0.5, -0.5) 5
        QUBIT_COORDS(0.5, 0.5) 6
        QUBIT_COORDS(0.5, 1.5) 7
        QUBIT_COORDS(0.5, 2.5) 8
        QUBIT_COORDS(1, 0) 9
        QUBIT_COORDS(1, 1) 10
        QUBIT_COORDS(1, 2) 11
        QUBIT_COORDS(1.5, -0.5) 12
        QUBIT_COORDS(1.5, 0.5) 13
        QUBIT_COORDS(1.5, 1.5) 14
        QUBIT_COORDS(1.5, 2.5) 15
        QUBIT_COORDS(2, 0) 16
        QUBIT_COORDS(2, 1) 17
        QUBIT_COORDS(2, 2) 18
        QUBIT_COORDS(2.5, 0.5) 19
        QUBIT_COORDS(2.5, 1.5) 20
        R 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20
        TICK
        H 0 1 3 7 9 13 16 17 18 19 20
        TICK
        M 6 8 12 14 15 0 7 13 19 20
        DETECTOR(0.5, 0.5, 0) rec[-10]
        DETECTOR(0.5, 2.5, 0) rec[-9]
        DETECTOR(1.5, -0.5, 0) rec[-8]
        DETECTOR(1.5, 1.5, 0) rec[-7]
        DETECTOR(1.5, 2.5, 0) rec[-6]
        SHIFT_COORDS(0, 0, 1)
        TICK
        R 0 7 13 19 20 6 8 12 14 15
        TICK
        H 0 6 7 8 12 13 14 20
        TICK
        ISWAP 2 6 3 7 4 8 9 13 10 14 17 20
        TICK
        C_XYZ 6 7 8 13 14 20
        S 2 3 4 9 10 17
        TICK
        ISWAP 1 3 2 5 4 7 6 9 8 11 10 13 14 17 18 20
        TICK
        C_XYZ 1 2 3 4 5 6 7 9 10 11 13 14 17 18
        S 8 20
        TICK
        ISWAP 0 2 3 6 7 10 9 12 11 14 13 16 15 18 17 19
        TICK
        C_XYZ 6 7 12 13 14 15 18 19
        S 0 2 3 9 10 11 16 17
        TICK
        ISWAP 2 6 3 7 9 13 10 14 11 15 16 19
        TICK
        C_XYZ 3 6 7 9 11 13 14 15 19
        S 2 10 16
        TICK
        M 5 7 8 13 15 1 6 14 19 20
        DETECTOR(0.5, -0.5, 0) rec[-20] rec[-10]
        DETECTOR(0.5, 1.5, 0) rec[-19] rec[-9]
        DETECTOR(0.5, 2.5, 0) rec[-8]
        DETECTOR(1.5, 0.5, 0) rec[-18] rec[-17] rec[-7]
        DETECTOR(1.5, 2.5, 0) rec[-16] rec[-6]
        SHIFT_COORDS(0, 0, 1)
        TICK
        R 1 6 14 19 20 5 7 8 13 15
        TICK
        REPEAT 49 {
            H 1 3 5 6 7 9 11 13 14 15 19
            TICK
            ISWAP 2 6 3 7 9 13 10 14 11 15 16 19
            TICK
            C_XYZ 6 7 13 14 15 19
            H 12 18
            S 2 3 9 10 11 16
            TICK
            ISWAP 0 2 3 6 7 10 9 12 11 14 13 16 15 18 17 19
            TICK
            C_XYZ 0 2 3 6 7 9 10 11 12 13 14 17 18
            H 4
            S 15 16 19
            TICK
            ISWAP 1 3 2 5 4 7 6 9 8 11 10 13 14 17 18 20
            TICK
            C_XYZ 6 7 8 13 14 20
            S 1 2 3 4 5 9 10 11 17 18
            TICK
            ISWAP 2 6 3 7 4 8 9 13 10 14 17 20
            TICK
            C_XYZ 6 7 8 13 14 20
            S 2 3 4 9 10 17
            TICK
            M 6 8 12 14 15 0 7 13 19 20
            DETECTOR(-0.5, 0.5, 0) rec[-14] rec[-5]
            DETECTOR(0.5, 0.5, 0) rec[-20] rec[-19] rec[-10]
            DETECTOR(0.5, 1.5, 0) rec[-15] rec[-13] rec[-4]
            DETECTOR(0.5, 2.5, 0) rec[-18] rec[-9]
            DETECTOR(1.5, -0.5, 0) rec[-17] rec[-8]
            DETECTOR(1.5, 0.5, 0) rec[-12] rec[-3]
            DETECTOR(1.5, 1.5, 0) rec[-16] rec[-7]
            DETECTOR(1.5, 2.5, 0) rec[-6]
            DETECTOR(2.5, 0.5, 0) rec[-2]
            DETECTOR(2.5, 1.5, 0) rec[-11] rec[-1]
            SHIFT_COORDS(0, 0, 1)
            TICK
            R 0 7 13 19 20 6 8 12 14 15
            TICK
            H 0 6 7 8 12 13 14 20
            TICK
            ISWAP 2 6 3 7 4 8 9 13 10 14 17 20
            TICK
            C_XYZ 6 7 8 13 14 20
            S 2 3 4 9 10 17
            TICK
            ISWAP 1 3 2 5 4 7 6 9 8 11 10 13 14 17 18 20
            TICK
            C_XYZ 1 2 3 4 5 6 7 9 10 11 13 14 17 18
            S 8 20
            TICK
            ISWAP 0 2 3 6 7 10 9 12 11 14 13 16 15 18 17 19
            TICK
            C_XYZ 6 7 12 13 14 15 18 19
            S 0 2 3 9 10 11 16 17
            TICK
            ISWAP 2 6 3 7 9 13 10 14 11 15 16 19
            TICK
            C_XYZ 3 6 7 9 11 13 14 15 19
            S 2 10 16
            TICK
            M 5 7 8 13 15 1 6 14 19 20
            DETECTOR(-0.5, 1.5, 0) rec[-14] rec[-5]
            DETECTOR(0.5, -0.5, 0) rec[-20] rec[-10]
            DETECTOR(0.5, 0.5, 0) rec[-15] rec[-13] rec[-4]
            DETECTOR(0.5, 1.5, 0) rec[-19] rec[-9]
            DETECTOR(0.5, 2.5, 0) rec[-8]
            DETECTOR(1.5, 0.5, 0) rec[-18] rec[-17] rec[-7]
            DETECTOR(1.5, 1.5, 0) rec[-11] rec[-3]
            DETECTOR(1.5, 2.5, 0) rec[-16] rec[-6]
            DETECTOR(2.5, 0.5, 0) rec[-12] rec[-2]
            DETECTOR(2.5, 1.5, 0) rec[-1]
            SHIFT_COORDS(0, 0, 1)
            TICK
            R 1 6 14 19 20 5 7 8 13 15
            TICK
        }
        H 1 3 6 7 9 11 13 14 15 19 20
        TICK
        ISWAP 2 6 3 7 9 13 10 14 11 15 16 19
        TICK
        C_XYZ 6 7 13 14 15 19
        H 12 18
        S 2 3 9 10 11 16
        TICK
        ISWAP 0 2 3 6 7 10 9 12 11 14 13 16 15 18 17 19
        TICK
        C_XYZ 2 6 10 12 14 16 18 19
        H 4
        S 0 3 7 9 11 13 15 17
        TICK
        M 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20
        DETECTOR(0.5, -0.5, 0) rec[-16]
        DETECTOR(0.5, 0, 0) rec[-31] rec[-30] rec[-19] rec[-16] rec[-15] rec[-12]
        DETECTOR(0.5, 1, 0) rec[-18] rec[-15] rec[-14] rec[-11]
        DETECTOR(0.5, 2, 0) rec[-29] rec[-17] rec[-14] rec[-13] rec[-10]
        DETECTOR(0.5, 2.5, 0) rec[-13]
        DETECTOR(1.5, -0.5, 0) rec[-28] rec[-9]
        DETECTOR(1.5, 0, 0) rec[-12] rec[-9] rec[-8] rec[-5]
        DETECTOR(1.5, 1, 0) rec[-27] rec[-11] rec[-8] rec[-7] rec[-4]
        DETECTOR(1.5, 2, 0) rec[-10] rec[-7] rec[-6] rec[-3]
        DETECTOR(1.5, 2.5, 0) rec[-6]
        OBSERVABLE_INCLUDE(0) rec[-19] rec[-18] rec[-17]
    """)


def test_cx_mpp2_circuit_details():
    circuit = make_hex_planar_surface_code(
        basis='X',
        distance=3,
        gate='CZ_MZZ',
        wiggle=False,
        rounds=100,
    ).circuit
    assert circuit == stim.Circuit("""
        QUBIT_COORDS(0, 0) 0
        QUBIT_COORDS(0, 1) 1
        QUBIT_COORDS(0, 2) 2
        QUBIT_COORDS(0.5, 0.5) 3
        QUBIT_COORDS(0.5, 1.5) 4
        QUBIT_COORDS(0.5, 2.5) 5
        QUBIT_COORDS(1, 0) 6
        QUBIT_COORDS(1, 1) 7
        QUBIT_COORDS(1, 2) 8
        QUBIT_COORDS(1.5, 0.5) 9
        QUBIT_COORDS(1.5, 1.5) 10
        QUBIT_COORDS(1.5, 2.5) 11
        QUBIT_COORDS(2, 0) 12
        QUBIT_COORDS(2, 1) 13
        QUBIT_COORDS(2, 2) 14
        QUBIT_COORDS(2.5, 0.5) 15
        QUBIT_COORDS(2.5, 1.5) 16
        R 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16
        TICK
        H 1 4 6 9 13 14 15 16
        TICK
        MPP Z0*Z3 Z6*Z9 Z15 Z1*Z4 Z7*Z10 Z13*Z16 Z2*Z5 Z11
        DETECTOR(0.5, 0.5, 0) rec[-8]
        DETECTOR(0.5, 2.5, 0) rec[-2]
        DETECTOR(1.5, 1.5, 0) rec[-4]
        DETECTOR(1.5, 2.5, 0) rec[-1]
        SHIFT_COORDS(0, 0, 1)
        TICK
        H 0 2 3 4 5 6 7 9 10 11 13 16
        TICK
        CZ 2 4 3 6 5 8 7 9 10 13 14 16
        TICK
        H 3 4 5 7 8 9 10 13 14 15 16
        TICK
        CZ 1 3 4 7 8 10 9 12 11 14 13 15
        TICK
        H 1 3 4 7 8 9 10 11 12 14 15
        TICK
        MPP Z0*Z3 Z6*Z9 Z12*Z15 Z1*Z4 Z7*Z10 Z16 Z5 Z8*Z11
        DETECTOR(0.5, 1.5, 0) rec[-5]
        DETECTOR(0.5, 2.5, 0) rec[-2]
        DETECTOR(1.5, 0.5, 0) rec[-7]
        DETECTOR(1.5, 2.5, 0) rec[-1]
        SHIFT_COORDS(0, 0, 1)
        TICK
        REPEAT 49 {
            H 1 3 4 5 7 8 9 10 11 12 14 15
            TICK
            CZ 1 3 4 7 8 10 9 12 11 14 13 15
            TICK
            H 3 4 7 8 9 10 11 13 14 15 16
            TICK
            CZ 2 4 3 6 5 8 7 9 10 13 14 16
            TICK
            H 0 2 3 4 5 6 7 9 10 13 16
            TICK
            MPP Z0*Z3 Z6*Z9 Z15 Z1*Z4 Z7*Z10 Z13*Z16 Z2*Z5 Z11
            DETECTOR(0.5, 0.5, 0) rec[-24] rec[-8]
            DETECTOR(0.5, 1.5, 0) rec[-21] rec[-5]
            DETECTOR(0.5, 2.5, 0) rec[-18] rec[-2]
            DETECTOR(1.5, 0.5, 0) rec[-23] rec[-7]
            DETECTOR(1.5, 1.5, 0) rec[-20] rec[-4]
            DETECTOR(1.5, 2.5, 0) rec[-17] rec[-1]
            DETECTOR(2.5, 0.5, 0) rec[-22] rec[-6]
            DETECTOR(2.5, 1.5, 0) rec[-19] rec[-3]
            SHIFT_COORDS(0, 0, 1)
            TICK
            H 0 2 3 4 5 6 7 9 10 11 13 16
            TICK
            CZ 2 4 3 6 5 8 7 9 10 13 14 16
            TICK
            H 3 4 5 7 8 9 10 13 14 15 16
            TICK
            CZ 1 3 4 7 8 10 9 12 11 14 13 15
            TICK
            H 1 3 4 7 8 9 10 11 12 14 15
            TICK
            MPP Z0*Z3 Z6*Z9 Z12*Z15 Z1*Z4 Z7*Z10 Z16 Z5 Z8*Z11
            DETECTOR(0.5, 0.5, 0) rec[-24] rec[-8]
            DETECTOR(0.5, 1.5, 0) rec[-21] rec[-5]
            DETECTOR(0.5, 2.5, 0) rec[-18] rec[-2]
            DETECTOR(1.5, 0.5, 0) rec[-23] rec[-7]
            DETECTOR(1.5, 1.5, 0) rec[-20] rec[-4]
            DETECTOR(1.5, 2.5, 0) rec[-17] rec[-1]
            DETECTOR(2.5, 0.5, 0) rec[-22] rec[-6]
            DETECTOR(2.5, 1.5, 0) rec[-19] rec[-3]
            SHIFT_COORDS(0, 0, 1)
            TICK
        }
        H 1 3 4 7 8 9 10 11 12 14 15
        TICK
        CZ 1 3 4 7 8 10 9 12 11 14 13 15
        TICK
        H 0 1 2 4 8 9 11 13 16
        TICK
        M 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16
        DETECTOR(0.5, 0, 0) rec[-33] rec[-17] rec[-14] rec[-11]
        DETECTOR(0.5, 1, 0) rec[-22] rec[-16] rec[-14] rec[-13] rec[-10]
        DETECTOR(0.5, 2, 0) rec[-27] rec[-15] rec[-13] rec[-12] rec[-9]
        DETECTOR(0.5, 2.5, 0) rec[-19] rec[-12]
        DETECTOR(1.5, 0, 0) rec[-24] rec[-11] rec[-8] rec[-5]
        DETECTOR(1.5, 1, 0) rec[-29] rec[-10] rec[-8] rec[-7] rec[-4]
        DETECTOR(1.5, 2, 0) rec[-18] rec[-9] rec[-7] rec[-6] rec[-3]
        DETECTOR(1.5, 2.5, 0) rec[-26] rec[-6]
        OBSERVABLE_INCLUDE(0) rec[-17] rec[-16] rec[-15]
    """)
