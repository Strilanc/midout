import dataclasses
from typing import Dict, FrozenSet, Optional, Set, TYPE_CHECKING

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

from midout.walking.diamonds import Diamond, DiamondState
from midout.walking.util import (
    Basis,
    Direction,
    get_compass_directions,
    get_default_color,
    MarkerType,
    Measurement,
    Qubit,
    qubit_coords,
)

if TYPE_CHECKING:
    from midout.walking.tiles import TileState

RESET_CIRCLE_RADIUS = 0.09


@dataclasses.dataclass(frozen=True)
class ResetLayer:
    """a reset layer, which adds expanding diamonds to a diamond state.

    which generally takes an existing set of only contracting diamonds
    and a set of corresponding tiles, and constructs new expanding diamonds,
    adds the measure qubits back into contracting diamonds,
    and records the reset gates
    """

    resets: Dict[Optional[Basis], FrozenSet[Qubit]]

    @staticmethod
    def from_tiles_and_diamonds(
        tile_state: 'TileState',
        contracting_diamond_state: DiamondState,
        direction: Optional[Direction] = None,
        reverse_expansion=False,
    ):
        new_diamonds = []
        resets = {Basis.z: set(), Basis.x: set()}

        for t in tile_state.tiles:
            # each of tile looks
            # and produce a new expanding diamond

            # first, lets get the corresponding contracting diamond for this tile
            possibly_multiple_diamonds = (
                contracting_diamond_state.get_diamonds_with_reinclude_qubit(t.measure_qubit)
            )
            if len(possibly_multiple_diamonds) > 1:
                raise ValueError
            elif len(possibly_multiple_diamonds) == 0:
                # this is normal for the first reset layer when only 1 basis it initialized
                # we just don't make a contracting diamond for the 1st cycle for this basis
                cont_diamond = None
                contracting_qubits = {t.measure_qubit}
            else:
                cont_diamond = possibly_multiple_diamonds.pop()
                contracting_qubits = set(cont_diamond.qubits)
            qubit_to_reinclude = None
            override_color = None

            # work out where we are:
            # boundary or bulk
            is_boundary = len(t.data_qubits) == 2
            is_trailing_boundary = False
            if direction is not None:
                check_qubit = t.measure_qubit + direction
                is_trailing_boundary = is_boundary and check_qubit in t.data_qubits

            # work out whether the expanding / contracting qubits need more qubits involved
            if direction is not None and is_boundary and not is_trailing_boundary:
                override_color = 'orange' if t.basis == Basis.x else 'plum'

                candidate_expansion_qubit = t.measure_qubit + direction
                # first work out which way along the boundary to expand
                perpendicular_directions = get_compass_directions(direction)[1:3]
                # this is the direction back toward the boundary in the given direction
                # such that direction + back_toward_boundary_direction should move you
                # straight along the boundary
                back_toward_boundary_direction = None
                for pd in perpendicular_directions:
                    if t.measure_qubit + pd in t.data_qubits:
                        back_toward_boundary_direction = pd
                if back_toward_boundary_direction is None:
                    raise ValueError

                if not reverse_expansion:
                    # work out the two relevant extra qubits
                    extra_expanding_qubit = (
                        t.measure_qubit + direction + back_toward_boundary_direction
                    )

                    qubit_to_reinclude = candidate_expansion_qubit
                    contracting_qubits.add(t.measure_qubit)
                    resets[t.basis].add(t.measure_qubit)

                    if extra_expanding_qubit - direction in tile_state.all_qubits:
                        expanding_qubits = {t.measure_qubit, extra_expanding_qubit}
                        resets[t.basis].add(t.measure_qubit)
                        resets[t.basis].add(extra_expanding_qubit)

                else:  # expansion type is backwards
                    contracting_qubits.add(t.measure_qubit)
                    resets[t.basis].add(t.measure_qubit)
                    contracting_qubits.add(candidate_expansion_qubit)
                    expanding_qubits = {candidate_expansion_qubit}
                    resets[t.basis].add(candidate_expansion_qubit)
                    # candidate_obvs_qubits[t.basis].add(candidate_expansion_qubit)

            else:  # not leading edge
                if direction is None or not reverse_expansion:
                    contracting_qubits = contracting_qubits | {t.measure_qubit}
                    expanding_qubits = {t.measure_qubit}
                    resets[t.basis].add(t.measure_qubit)

                    # if we're forward expanding and we're a bulk stabilizer on the leading edge,
                    # we have to extend our expanding diamond out into the new qubits as well
                    #   the if cases we've already passed through mean we only have to check
                    #   if we're expanding and not a boundary
                    #   to know we're a bulk stab and we're fwd expanding
                    if direction is not None and not is_boundary:
                        # now find out if we're near the leading edge
                        # check both perpendicular directions - if either are out of the tile,
                        # then we are near the leading edge, and that's the qubit to include
                        perpendicular_directions = get_compass_directions(direction)[1:3]
                        for pd in perpendicular_directions:
                            relevant_qubit = t.measure_qubit + direction + pd
                            if relevant_qubit not in tile_state.all_qubits:
                                expanding_qubits.add(relevant_qubit)
                                resets[t.basis].add(relevant_qubit)

                else:  # expansion happening and is reversed
                    exp_offset_qubit = t.measure_qubit + direction + direction
                    expanding_qubits = {exp_offset_qubit}
                    resets[t.basis].add(exp_offset_qubit)
                    if not is_trailing_boundary:
                        contracting_qubits = contracting_qubits | {t.measure_qubit}
                        resets[t.basis].add(t.measure_qubit)

            # actually make the new diamonds we found the qubits for
            exp_diamond = Diamond(
                qubits=frozenset(expanding_qubits),
                basis=t.basis,
                marker_type=MarkerType.expanding,
                override_color=override_color,
                qubit_to_reinclude=qubit_to_reinclude,
            )
            new_diamonds.append(exp_diamond)
            if cont_diamond is not None:
                cont_diamond = Diamond(
                    qubits=frozenset(contracting_qubits),
                    basis=cont_diamond.basis,
                    marker_type=MarkerType.contracting,
                    override_color=override_color,
                    measurements=cont_diamond.measurements,
                    duid=cont_diamond.duid,
                )
                new_diamonds.append(cont_diamond)

        # now sort out expanding the logical observable
        # the candidates are any qubit that's been reset in the right basis
        # if they're not covered by 2 detectors, and closer than 1 away from the LO, add them
        new_obvs = contracting_diamond_state.observables.copy()
        for b, obvs in contracting_diamond_state.observables.items():
            for q in resets[b]:
                # check if it's near the observable
                if any([np.abs(qo - q) < 1 for qo in obvs]):
                    diamonds_that_have_this_qubit = [d for d in new_diamonds if q in d.qubits]
                    if len(diamonds_that_have_this_qubit) == 0:
                        raise ValueError("a qubit was reset but not included in any diamonds")
                    elif len(diamonds_that_have_this_qubit) == 1:
                        new_obvs[b] = frozenset(new_obvs[b] | {q})

        new_ds = DiamondState(diamonds=frozenset(new_diamonds), observables=new_obvs)
        return new_ds, ResetLayer(resets={b: frozenset(s) for b, s in resets.items()})

    @staticmethod
    def make_init_layer(tile_state: "TileState", init_basis: Optional[Basis]):
        """make the resets and contracting diamonds for the initialization of the code."""
        diamonds = []
        qubits_to_reset: Set[Qubit] = set()
        for t in tile_state.tiles:
            if init_basis is None:
                # we're doing a dummy init
                # include the measurements you'd 'expect' to see on each diamond
                # basically, measure the measure qubit in the previous round
                measurements = frozenset([Measurement(qubit=t.measure_qubit, index=-1)])
            elif init_basis == t.basis:
                # we're doing a real initialization,
                # don't include any measurements so the
                # 1st round detector singletons are correct
                measurements = None
            else:
                continue
            cont_diamond = Diamond(
                qubits=frozenset(t.data_qubits),
                basis=t.basis,
                marker_type=MarkerType.contracting,
                override_color=None,
                measurements=measurements,
                qubit_to_reinclude=t.measure_qubit,
            )
            diamonds.append(cont_diamond)
            qubits_to_reset.update(t.data_qubits)

        if init_basis is None:
            observables = tile_state.observables
        else:
            observables = {init_basis: tile_state.observables[init_basis]}
        contracting_ds = DiamondState(diamonds=frozenset(diamonds), observables=observables)
        return contracting_ds, ResetLayer(resets={init_basis: frozenset(qubits_to_reset)})

    def plot(self, ax=None):
        """plot this layer's resets."""
        if ax is None:
            _, ax = plt.subplots()

        for b, qubits in self.resets.items():
            for q in qubits:
                circle = mpl.patches.Circle(
                    qubit_coords(q),
                    radius=RESET_CIRCLE_RADIUS,
                    facecolor=get_default_color(basis=b),
                    edgecolor="grey",
                    linewidth=2,
                    zorder=51,
                )
                ax.add_patch(circle)

        if not ax.yaxis_inverted():
            ax.invert_yaxis()
        ax.set_ylabel("Y = Imag = 2nd coord")
        ax.xaxis.tick_top()
        ax.xaxis.set_label_position('top')
        ax.set_xlabel("X = Real = 1st coord")

        ax.autoscale_view()
        ax.grid(True)
        ax.set_aspect('equal')

        return ax
