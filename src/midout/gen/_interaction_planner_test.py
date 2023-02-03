import pytest
import stim

from midout.gen._interaction_planner import XZ_BASE_TRANSITION_MAP
from midout.gen._builder import Builder


def test_interaction_helper_single():
    builder = Builder.for_qubits([0, 1])
    builder.circuit.clear()

    with builder.plan_interactions(layer_count=1) as planner:
        planner.pcp('X', 'Y', 0, 1, layer=0)
        with pytest.raises(ValueError, match="Collision"):
            planner.pcp('X', 'Y', 0, 3, layer=0)
    assert builder.circuit == stim.Circuit("""
        C_XYZ 1
        H 0
        TICK
        CZ 0 1
        TICK
        C_ZYX 1
        H 0
    """)


def test_interaction_helper_rep_code():
    builder = Builder.for_qubits(range(7))
    builder.circuit.clear()
    with builder.plan_interactions(layer_count=2,
                                   end_orientations={3: 'ZY'}) as planner:
        planner.pcp('Z', 'X', 0, 1, layer=0)
        planner.pcp('X', 'Z', 5, 4, layer=0)
        planner.cx(2, 3, layer=0)
        planner.pcp('Z', 'X', 2, 1, layer=1)
        planner.pcp('Z', 'X', 6, 5, layer=1)
        planner.pcp('Z', 'X', 4, 3, layer=1)
    assert builder.circuit == stim.Circuit("""
        H 1 3 5
        TICK
        CZ 0 1 2 3 4 5
        TICK
        CZ 1 2 3 4 5 6
        TICK
        H 1 5
        H_XY 3
    """)


def test_interaction_helper_merge_layers():
    builder = Builder.for_qubits(range(7))
    builder.circuit.clear()
    with builder.plan_interactions(layer_count=2) as planner:
        planner.pcp('X', 'X', 1, 2, layer=0)
        planner.pcp('X', 'X', 0, 1, layer=1)
    assert builder.circuit == stim.Circuit("""
        H 0 1 2
        TICK
        CZ 1 2
        TICK
        CZ 0 1
        TICK
        H 0 1 2
    """)

    builder.circuit.clear()
    with builder.plan_interactions(
            layer_count=2,
            end_orientations={0: 'ZX', 1: 'ZX', 2: 'ZX'},
    ) as planner:
        planner.pcp('X', 'X', 1, 2, layer=0)
        planner.pcp('X', 'X', 0, 1, layer=1)
    assert builder.circuit == stim.Circuit("""
        H 0 1 2
        TICK
        CZ 1 2
        TICK
        CZ 0 1
    """)

    builder.circuit.clear()
    with builder.plan_interactions(
            layer_count=2,
            start_orientations={0: 'ZX', 1: 'ZX', 2: 'ZX'},
    ) as planner:
        planner.pcp('X', 'X', 1, 2, layer=0)
        planner.pcp('X', 'X', 0, 1, layer=1)
    assert builder.circuit == stim.Circuit("""
        CZ 1 2
        TICK
        CZ 0 1
        TICK
        H 0 1 2
    """)


def test_basis_transition_map():
    for k, v in XZ_BASE_TRANSITION_MAP.items():
        b0, b1 = k.split(' -> ')
        for p0, p1 in zip(b0, b1):
            # Conversion fails if measurement is not deterministic, verifying basis transformation.
            try:
                stim.Circuit(f"""
                    R{p0} 0
                    {v} 0
                    M{p1} 0
                    DETECTOR rec[-1]
                """).detector_error_model()
                success = True
            except ValueError:
                success = False
            assert success, f"{k} not achieved by {v}"
