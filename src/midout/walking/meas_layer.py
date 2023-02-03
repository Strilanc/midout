import dataclasses
from typing import Dict, FrozenSet, List, Optional, Set, Tuple

import matplotlib as mpl
import matplotlib.pyplot as plt

from midout.walking.diamonds import Diamond, DiamondState
from midout.walking.tiles import TileState
from midout.walking.util import (
    Basis,
    get_default_color,
    MarkerType,
    Measurement,
    Qubit,
    qubit_coords,
)

MEAS_CIRCLE_RADIUS = 0.09


@dataclasses.dataclass(frozen=True)
class MeasLayer:

    measurements: Dict[Basis, FrozenSet[Measurement]]
    detectors: FrozenSet[Tuple[FrozenSet[Measurement], int]]
    observable_includes: Dict[Basis, FrozenSet[Measurement]]

    @staticmethod
    def make_measure_layer(tile_state: TileState, diamond_state: DiamondState, index: int):
        """Given a set of measurements and a DiamondState, work out the detectors and obvs includes.

        Basically,

        All the contracting diamonds should be fully measured, complain if not
        All the expanding diamonds should probably be partially measured,
            record in those diamonds which qubit is measured in this meas round,
            if none are measured, leave this diamond alone,
            we assume it already knows how to handle the next reset layer

        Also, convert all expanding diamonds into contracting diamonds, and return those
        """
        qubits_to_keep = tile_state.all_qubits
        qubits_to_reinclude = tile_state.measure_qubits

        qubits_to_dropout = diamond_state.all_qubits.difference(qubits_to_keep)
        qubits_to_measure = qubits_to_reinclude | qubits_to_dropout

        measurements: Dict[Basis, Set[Measurement]] = {Basis.x: set(), Basis.z: set()}
        detectors: List[Tuple[FrozenSet[Measurement], int]] = []
        observable_includes: Dict[Basis, Set[Measurement]] = {}
        contracting_diamonds = []

        for d in diamond_state.diamonds:
            if d.marker_type == MarkerType.contracting:
                # make sure all qubits are in qubits to measure, else complain
                this_detector = set(d.measurements) if d.measurements is not None else set()
                for q in d.qubits:
                    if q not in qubits_to_measure:
                        raise ValueError(
                            "A contracting diamond did not have all its qubits measured: "
                            f"Qubit {q} in Diamond {d} wasn't included in qubits_to_measure "
                            f"{qubits_to_measure}"
                        )
                    measurements[d.basis].add(Measurement(index, q))
                    this_detector.add(Measurement(index, q))
                detectors.append((frozenset(this_detector), d.duid))

            elif d.marker_type == MarkerType.expanding:
                # measure the qubits you're supposed to measure
                this_diamond_measurements = set(
                    Measurement(index, q) for q in d.qubits if q in qubits_to_measure
                )
                measurements[d.basis].update(this_diamond_measurements)
                # work out if we are measuring a qubit to re-include later
                possible_measure_qubits = [q for q in d.qubits if q in qubits_to_reinclude]
                if len(possible_measure_qubits) == 0:
                    # this diamond isn't getting any of its qubits measured,
                    # if it already thinks it has a measure qubit to include next round, keep that
                    measure_qubit = d.qubit_to_reinclude
                elif len(possible_measure_qubits) != 1:
                    raise ValueError(
                        "Something Went Wrong: measure layer found multiple "
                        "overlaps between the tile measure qubits and a single diamond. "
                        "This probably means that the diamonds for this layer were messed up. "
                        "Try plotting the diamond right before measurement."
                    )
                else:
                    measure_qubit = possible_measure_qubits[0]
                # make a new contracting diamond for the rest
                contracting_diamond = Diamond(
                    qubits=frozenset([q for q in d.qubits if q not in qubits_to_measure]),
                    basis=d.basis,
                    marker_type=MarkerType.contracting,
                    measurements=frozenset(this_diamond_measurements),
                    qubit_to_reinclude=measure_qubit,
                    duid=d.duid,
                )
                contracting_diamonds.append(contracting_diamond)
            elif d.marker_type == MarkerType.error:
                pass

        new_observables: Dict[Basis, Set[Qubit]] = {}
        for b, obv in diamond_state.observables.items():
            new_observables[b] = set()
            observable_includes[b] = set()
            for q in obv:
                # each obvs qubit is either being measured, or surviving
                for m in measurements[b]:
                    if m.qubit == q:
                        observable_includes[b].add(m)
                        break
                else:
                    new_observables[b].add(q)
            # check we didn't find any 'survivors' that are measured in the other basis
            if len(new_observables[b].intersection({m.qubit for m in measurements[b.flip()]})) != 0:
                raise ValueError("a qubit in an observable was measured in a non-commuting basis")

        contracting_ds = DiamondState(
            diamonds=frozenset(contracting_diamonds),
            observables={b: frozenset(s) for b, s in new_observables.items()},
        )
        m_layer = MeasLayer(
            measurements={b: frozenset(s) for b, s in measurements.items()},
            detectors=frozenset(detectors),
            observable_includes={b: frozenset(s) for b, s in observable_includes.items()},
        )
        return contracting_ds, m_layer

    @staticmethod
    def make_final_measure_layer(
        tile_state: TileState,
        diamond_state: DiamondState,
        index: int,
        basis: Optional[Basis] = None,
    ):
        """perform the final measurements for a memory experiment.

        measure all the remaining qubits in the given basis
        and construct the appropriate detectors and observable includes

        the input diamond state should be a just-measured contracting-only set of diamonds
        """

        # measurements are easy - measure everyone in the given basis
        measurements = {basis: frozenset(Measurement(index, q) for q in diamond_state.all_qubits)}
        # detectors are also easy - for diamonds of the basis being measured,
        # form a detector out of the prev measurement and all the current qubits
        # if we're doing a dodgy final measurement, make no detectors
        detectors: List[Tuple[FrozenSet[Measurement], int]] = []
        for d in diamond_state.diamonds:
            if d.basis == basis:
                this_detector = set(d.measurements) if d.measurements is not None else set()
                this_detector.update({Measurement(index, q) for q in d.qubits})
                detectors.append((frozenset(this_detector), d.duid))

        # observable includes are just include the qubits on the observable of matching basis.
        # if there's a non-commuting observable, complain.
        observable_includes = {}
        if basis is not None:
            if basis.flip() in diamond_state.observables:
                raise ValueError(
                    f"Non-commuting observable {basis.flip()} "
                    f"in final measurement round in basis {basis}"
                )
            observable_includes[basis] = frozenset(
                Measurement(index, q) for q in diamond_state.observables[basis]
            )
        else:
            observable_includes[basis]: FrozenSet[Measurement] = frozenset()

        return MeasLayer(
            measurements=measurements,
            detectors=frozenset(detectors),
            observable_includes=observable_includes,
        )

    def plot(self, ax=None):
        """plot this layer's measurements."""
        if ax is None:
            _, ax = plt.subplots()

        for b, ms in self.measurements.items():
            for m in ms:
                circle = mpl.patches.Circle(
                    qubit_coords(m.qubit),
                    radius=MEAS_CIRCLE_RADIUS,
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
