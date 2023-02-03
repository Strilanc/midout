import functools
import pathlib
from typing import Optional, Callable, Dict, Tuple

import sinter
import stim

from midout import gen
from midout._circuit_case import CircuitCase
from midout.planar._classic_surface_code import make_square_planar_cx_code, \
    make_square_planar_cz_code
from midout.planar._cxswap_surface_code import make_square_planar_cxswap_code
from midout.planar._iswap_surface_code import make_square_planar_iswap_code
from midout.planar._three_coupler_surface_codes import \
    make_hex_planar_surface_code
from midout.toric._hex_toric_cxmpp import make_hex_toric_cxmpp_code
from midout.toric._toric_code import make_square_toric_cx_code
from midout.toric._toric_heavyhex import make_heavyhex_toric_cx_code
from midout.toric._toric_semiheavyhex_3cycle import \
    make_semiheavyhex_toric_cx_3cycle_code
from midout.walking._make_walking_circuit_cases import make_walking_code


def make_construction_dict() -> Dict[str, Callable]:
    return {
        '4-ISWAP': make_square_planar_iswap_code,
        '4-ISWAP-ALT': functools.partial(make_square_planar_cxswap_code, use_iswaps=True),
        '4-CXSWAP': make_square_planar_cxswap_code,
        '4-CZ': make_square_planar_cz_code,
        '4-CX': make_square_planar_cx_code,
        '3-CX': functools.partial(make_hex_planar_surface_code, wiggle=False, gate='CX'),
        '3-CZ_MZZ': functools.partial(make_hex_planar_surface_code, wiggle=False, gate='CZ_MZZ'),
        '3-CX_MXX_MZZ': functools.partial(make_hex_planar_surface_code, wiggle=False, gate='CX_MXX_MZZ'),
        '3-CXSWAP': functools.partial(make_hex_planar_surface_code, wiggle=False, gate='CXSWAP'),
        '3-CZ': functools.partial(make_hex_planar_surface_code, wiggle=False, gate='CZ'),
        '3-ISWAP': functools.partial(make_hex_planar_surface_code, wiggle=False, gate='ISWAP'),
        '3-CX-wiggle': functools.partial(make_hex_planar_surface_code, wiggle=True, gate='CX'),
        '3-CXSWAP-wiggle': functools.partial(make_hex_planar_surface_code, wiggle=True, gate='CXSWAP'),
        '3-CZ-wiggle': functools.partial(make_hex_planar_surface_code, wiggle=True, gate='CZ'),
        '3-ISWAP-wiggle': functools.partial(make_hex_planar_surface_code, wiggle=True, gate='ISWAP'),
        'TORIC-4-CX': make_square_toric_cx_code,
        'TORIC-3-CX_MXX_MZZ': make_hex_toric_cxmpp_code,
        'TORIC-3_HEAVY-CX': make_heavyhex_toric_cx_code,
        'TORIC-3_SEMI_HEAVY-CX': make_semiheavyhex_toric_cx_3cycle_code,
        'WIGGLING-CX':  functools.partial(make_walking_code, strategy='wiggling', gate='CX'),
        'WIGGLING-CZ': functools.partial(make_walking_code, strategy='wiggling', gate='CZ'),
        'GLIDING-CX': functools.partial(make_walking_code, strategy='gliding', gate='CX'),
        'GLIDING-CZ': functools.partial(make_walking_code, strategy='gliding', gate='CZ'),
        'SLIDING-CX': functools.partial(make_walking_code, strategy='sliding', gate='CX'),
        'SLIDING-CZ': functools.partial(make_walking_code, strategy='sliding', gate='CZ'),
    }


CONSTRUCTIONS = make_construction_dict()


def hide_long_range_tiles(patch: gen.Patch, remove_entirely: bool = False) -> gen.Patch:
    back_tiles = []
    tiles = []
    for tile in patch.tiles:
        max_distance = max(abs(a - b) for a in tile.ordered_data_qubits for b in tile.ordered_data_qubits if a is not None and b is not None)
        if max_distance > 4:
            if not remove_entirely:
                back_tiles.append(tile)
        else:
            tiles.append(tile)
    return gen.Patch(back_tiles + tiles, do_not_sort=True)


def make_requested_surface_code(
        *,
        basis: str,
        distance: int,
        noise: gen.NoiseModel,
        style: str,
        rounds: int,
        debug_out_dir: Optional[pathlib.Path] = None,
) -> Tuple[CircuitCase, stim.Circuit]:
    if style not in CONSTRUCTIONS:
        raise NotImplementedError(f'{style}=')
    result: CircuitCase = CONSTRUCTIONS[style](
        distance=distance,
        basis=basis,
        rounds=rounds,
    )
    assert isinstance(result, CircuitCase)

    main_patch = None
    if debug_out_dir is not None:
        patches = [hide_long_range_tiles(patch) for patch in result.patches]
        main_patch = hide_long_range_tiles(patches[0], True)

        path = debug_out_dir / "tiles.svg"
        with open(path, "w") as f:
            print(gen.patch_svg_viewer(
                patches,
                show_order=result.show_patch_order,
                show_measure_qubits=result.show_patch_measure_qubits,
            ), file=f)
        print(f'wrote file://{path.absolute()}')

        path = debug_out_dir / "ideal_circuit.html"
        with open(debug_out_dir / "ideal_circuit.html", "w") as f:
            print(gen.stim_circuit_html_viewer(result.circuit, patch=main_patch), file=f)
        print(f'wrote file://{path.absolute()}')

    noisy_circuit = noise.noisy_circuit(result.circuit)

    if debug_out_dir is not None:
        path = debug_out_dir / "circuit.html"
        with open(debug_out_dir / "circuit.html", "w") as f:
            print(gen.stim_circuit_html_viewer(noisy_circuit, patch=main_patch), file=f)
        print(f'wrote file://{path.absolute()}')

    return result, noisy_circuit


def xz_piece_error_rate(p_combo: float, *, pieces: float, combo: bool) -> float:
    if not combo:
        return sinter.shot_error_rate_to_piece_error_rate(p_combo, pieces=pieces)
    p_solo = 1 - (1 - p_combo)**0.5
    unit_solo = sinter.shot_error_rate_to_piece_error_rate(p_solo, pieces=pieces)
    unit_combo = unit_solo*(2 - unit_solo)
    return unit_combo
