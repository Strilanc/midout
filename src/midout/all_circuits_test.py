import itertools
from typing import Set

import pytest
import stim

from midout import gen
from midout.all_circuits import make_requested_surface_code, CONSTRUCTIONS, \
    xz_piece_error_rate


@pytest.mark.parametrize("distance,basis,style,round_offset", itertools.product(
    range(3, 10),
    "XZ",
    sorted(CONSTRUCTIONS.keys()),
    [0, 1],
))
def test_all_surface_code_constructions(distance: int, basis: str, style: str, round_offset: int):
    if "TORIC" in style and (distance % 2 == 1 or round_offset == 1):
        return
    if style == '4-ISWAP' and round_offset == 1:
        return

    details, circuit = make_requested_surface_code(
        distance=distance,
        basis=basis,
        noise=gen.NoiseModel.uniform_depolarizing(1e-3),
        style=style,
        rounds=distance * 2 + round_offset,
    )

    dem = circuit.detector_error_model(decompose_errors=True)
    assert dem is not None

    actual_distance = len(circuit.shortest_graphlike_error())
    expected_distance = distance
    assert actual_distance == expected_distance

    expected_gates = details.expected_interactions | {
        'R',
        'RX',
        'X',
        'M',
        'MX',
        'DETECTOR',
        'OBSERVABLE_INCLUDE',
        'QUBIT_COORDS',
        'SHIFT_COORDS',
        'TICK',
        'Z_ERROR',
        'X_ERROR',
        'DEPOLARIZE1',
        'DEPOLARIZE2',
    }
    single_qubit_rotations = {
        'SQRT_X',
        'S',
        'C_XYZ',
        'C_ZYX',
        'H',
        'H_YZ',
    }
    if 'CX' not in style or 'TORIC' in style:
        expected_gates |= single_qubit_rotations
    assert gates_used_by_circuit(circuit) <= expected_gates


def gates_used_by_circuit(circuit: stim.Circuit) -> Set[str]:
    out = set()
    for instruction in circuit:
        if isinstance(instruction, stim.CircuitRepeatBlock):
            out |= gates_used_by_circuit(instruction.body_copy())
        elif instruction.name in ['CX', 'CY', 'CZ', 'XCZ', 'YCZ']:
            targets = instruction.targets_copy()
            for k in range(0, len(targets), 2):
                if targets[k].is_measurement_record_target or targets[k + 1].is_measurement_record_target:
                    out.add('feedback')
                elif targets[k].is_sweep_bit_target or targets[k + 1].is_sweep_bit_target:
                    out.add('sweep')
                else:
                    out.add(instruction.name)
        elif instruction.name == 'MPP':
            op = 'M'
            targets = instruction.targets_copy()
            is_continuing = True
            for t in targets:
                if t.is_combiner:
                    is_continuing = True
                    continue
                p = 'X' if t.is_x_target else 'Y' if t.is_y_target else 'Z' if t.is_z_target else '?'
                if is_continuing:
                    op += p
                    is_continuing = False
                else:
                    if op == 'MZ':
                        op = 'M'
                    out.add(op)
                    op = 'M' + p
            if op:
                if op == 'MZ':
                    op = 'M'
                out.add(op)
        else:
            out.add(instruction.name)
    return out


def test_xz_piece_error_rate():
    assert xz_piece_error_rate(0.74, pieces=10, combo=False) > 0.74
    assert xz_piece_error_rate(0.74, pieces=10, combo=True) < 0.74
    assert xz_piece_error_rate(0.75, pieces=10, combo=True) == 0.75
    assert xz_piece_error_rate(0.76, pieces=10, combo=True) > 0.76
