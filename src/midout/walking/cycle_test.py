import pytest

from midout.walking.cycle import Cycle
from midout.walking.tiles import TileState
from midout.walking.util import Basis, DOWN_LEFT, DOWN_RIGHT, UP_LEFT, UP_RIGHT


@pytest.mark.parametrize("top_basis", [Basis.x, Basis.z])
@pytest.mark.parametrize("observable_side", [True, False])
@pytest.mark.parametrize("reverse", [True, False])
@pytest.mark.parametrize("direction", [UP_LEFT, UP_RIGHT, DOWN_LEFT, DOWN_RIGHT])
def test_cycle(top_basis, observable_side, reverse, direction):
    """make sure the cycle builds and the measurement layer checks pass.

    Note does not check if you made hook error mistakes explicitly
    """
    ts = TileState.make_surface_code(
        distance=7, top_boundary_basis=top_basis, observables_at_bottom=observable_side
    )
    flip_left_right = top_basis == Basis.x and not reverse
    _ = Cycle(
        index=0,
        input_tile_state=ts,
        input_contracting_diamond_state=ts.make_init_diamonds(),
        direction=direction,
        dont_make_measurements=False,
        flip_left_right=flip_left_right,
        reverse_expansion=reverse,
    )
