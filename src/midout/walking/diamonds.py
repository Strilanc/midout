import dataclasses
import itertools
from typing import ClassVar, Dict, FrozenSet, Iterable, Optional, TYPE_CHECKING, Union

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

from midout.walking.util import (
    Basis,
    centroid,
    get_default_color,
    get_observable_marker,
    line_ordered_qubit_coords,
    MarkerType,
    Measurement,
    Qubit,
    rot_order_qubit_coords,
    shoelace_area,
)

if TYPE_CHECKING:
    from midout.walking.cnot_layer import CnotLayer
    from midout.walking.stim_circuit_builder import StimCircuitBuilder

DIAMOND_LINE_HALF_WIDTH = 0.09
DIAMOND_ERROR_HALF_WIDTH = 0.05
DIAMOND_CIRCLE_RADIUS = DIAMOND_LINE_HALF_WIDTH * 2.0 / 3.0


@dataclasses.dataclass(frozen=True)
class Diamond:
    """A diamond, defining a time slice of a detecting region during a surface code cycle.

    Attributes:
        qubits: a tuple defining the qubits included in the diamond
        basis: the basis of the stabilizer this diamons is representing
        expanding: whether this diamond is
            expanding: part of a new detecting region, and should be grown into a stabilizer
            contracting: part of an old detecting region, and should be shrunk into the meas qubit
            This is important because following contracting diamonds indicates
            how we should define detectors
        override_color: an absolute override on the color of the diamond
        measurements: an optional series of previously included measurements
            helpful because a fully contracted diamond locally knows what detector to create
        qubit_to_reinclude: a special extra field for DiamondStates between measurements and resets
            a reset layer generally needs to re-include the measure qubit for the tile it is reping
            helpful because now the diamond knows locally when it hits the reset layer what qubit
            should be reset and reincluded
        duid: a unique identifier for this detecting region
            this should be preserved when diamonds are transformed by gate layers
            and included when the diamond is finally measured and forms a detector
    """

    qubits: FrozenSet[Qubit]
    basis: Basis
    marker_type: MarkerType
    override_color: Optional[str] = None
    measurements: Optional[FrozenSet[Measurement]] = None
    qubit_to_reinclude: Optional[Qubit] = None
    duid: Optional[int] = None

    duid_counter: ClassVar[Optional['StimCircuitBuilder']] = None

    @classmethod
    def reset_duid_counter(cls):
        cls.duid_counter = itertools.count(start=1)

    @classmethod
    def delete_duid_counter(cls):
        cls.duid_counter = None

    def __post_init__(self):
        if self.duid_counter is not None and self.duid is None:
            object.__setattr__(self, 'duid', next(self.duid_counter))

    @property
    def color(self):
        return self.override_color or get_default_color(self.basis, self.marker_type)

    def __len__(self):
        return len(self.qubits)

    def poly_coords(self):
        """returns coords to plot a polygon representing this diamond."""
        qubit_coords = np.array(rot_order_qubit_coords(self.qubits))
        if np.isclose(shoelace_area(qubit_coords), 0):
            # we have a number of points that are co-linear
            edge_spacing = (
                DIAMOND_ERROR_HALF_WIDTH
                if self.marker_type == MarkerType.error
                else DIAMOND_LINE_HALF_WIDTH
            )
            corner_spacing = edge_spacing * np.sqrt(2)  # distance from qubit center to corner

            x_then_y_sorted_coords = qubit_coords[
                np.lexsort([qubit_coords[:, 1], qubit_coords[:, 0]])
            ]
            left_q_coords = x_then_y_sorted_coords[0]
            right_q_coords = x_then_y_sorted_coords[-1]
            if left_q_coords[0] == right_q_coords[0] and len(qubit_coords) != 1:
                # they're aligned along the x axis
                coords = [
                    [left_q_coords[0] - edge_spacing, left_q_coords[1] - edge_spacing],
                    [left_q_coords[0] + edge_spacing, left_q_coords[1] - edge_spacing],
                    [right_q_coords[0] + edge_spacing, right_q_coords[1] + edge_spacing],
                    [right_q_coords[0] - edge_spacing, right_q_coords[1] + edge_spacing],
                ]
            elif left_q_coords[1] == right_q_coords[1] and len(qubit_coords) != 1:
                # they're aligned along the y axis
                # because we did lex sort, left means bottom and right means top
                coords = [
                    [left_q_coords[0] - edge_spacing, left_q_coords[1] + edge_spacing],
                    [left_q_coords[0] - edge_spacing, left_q_coords[1] - edge_spacing],
                    [right_q_coords[0] + edge_spacing, right_q_coords[1] - edge_spacing],
                    [right_q_coords[0] + edge_spacing, right_q_coords[1] + edge_spacing],
                ]
            elif right_q_coords[1] > left_q_coords[1]:
                # they're diagonal and right qubit is above the left qubit
                coords = [
                    [left_q_coords[0] - corner_spacing, left_q_coords[1]],
                    [left_q_coords[0], left_q_coords[1] - corner_spacing],
                    [right_q_coords[0] + corner_spacing, right_q_coords[1]],
                    [right_q_coords[0], right_q_coords[1] + corner_spacing],
                ]
            else:
                # they're diagonal and right qubit is below the left qubit
                coords = [
                    [left_q_coords[0] - corner_spacing, left_q_coords[1]],
                    [left_q_coords[0], left_q_coords[1] + corner_spacing],
                    [right_q_coords[0] + corner_spacing, right_q_coords[1]],
                    [right_q_coords[0], right_q_coords[1] - corner_spacing],
                ]
        else:
            # the points are not co-linear, so we can just use the qubit coordinates
            coords = qubit_coords
        return coords, centroid(coords)

    def interpolate_coords_through_cnot_layer(self, cnot_layer: 'CnotLayer', p: float):
        """returns the coords of a polygon representing this diamond part way through a not layer.

        Args:
            cnot_layer: the cnot layer to apply to this diamond
            p: float between 0 adn 1, how far through the layer to draw the poly

        TODO: something about the rot sorting, which changes randomly
        """
        qubit_cnot_map = cnot_layer.track_stabilizers_dict(self.qubits, self.basis)
        interp_coords = []
        for q in self.qubits:
            final_qubits = qubit_cnot_map[q]
            for qf in final_qubits:
                # add a point between q and qf
                interp_coords.append((1 - p) * q + p * qf)
        return rot_order_qubit_coords(interp_coords), centroid(interp_coords)


@dataclasses.dataclass(frozen=True)
class DiamondState:
    """A set of diamonds indicating the state of a code during the cycle."""

    diamonds: FrozenSet[Diamond]
    observables: Dict[Basis, FrozenSet[Qubit]]

    @property
    def all_qubits(self):
        qubits = set()
        for d in self.diamonds:
            qubits.update(d.qubits)
        for obv in self.observables.values():
            qubits.update(obv)
        return qubits

    def add_diamond(self, diamond):
        return DiamondState(
            diamonds=frozenset(self.diamonds | {diamond}), observables=self.observables
        )

    def replace_diamond(self, old_diamond: Diamond, *args, **kwargs):
        """lets you evolve a given diamond in this diamond state."""
        if old_diamond not in self.diamonds:
            raise ValueError
        new_diamond = dataclasses.replace(old_diamond, *args, **kwargs)
        return DiamondState(
            diamonds=frozenset(self.diamonds.difference({old_diamond}) | {new_diamond}),
            observables=self.observables,
        )

    def get_diamonds_with_qubit(self, q: Qubit):
        return set([d for d in self.diamonds if q in d.qubits])

    def get_diamonds_with_reinclude_qubit(self, q: Qubit):
        return set([d for d in self.diamonds if q == d.qubit_to_reinclude])

    def get_diamond_with_duid(self, duid: int):
        possible_diamonds = [d for d in self.diamonds if d.duid == duid]
        if len(possible_diamonds) > 1:
            raise ValueError(
                f"Got multiple diamonds with the same duid, "
                f"which shouldn't happen: {possible_diamonds}"
            )
        elif len(possible_diamonds) == 0:
            return None
        return possible_diamonds[0]

    def get_diamonds(
        self,
        filter_basis: Optional[Basis] = None,
        filter_marker: Optional[Union[MarkerType, Iterable[MarkerType]]] = None,
        filter_size: Optional[int] = None,
        filter_qubit: Optional[Qubit] = None,
    ):
        """return diamonds with given characteristics. None doesn't filter on that argument.

        Args:
            filter_basis: if not None, returns only diamonds with a matching basis
            filter_marker: if not None, returns only diamonds with a matching marker type
            filter_size: if not None, returns only diamonds with fewer than the given num of qubits
            filter_qubit: if not None, returns only diamonds with this qubit
        """
        if filter_marker is not None and isinstance(filter_marker, MarkerType):
            filter_marker = [filter_marker]
        return set(
            d
            for d in self.diamonds
            if (filter_basis is None or d.basis == filter_basis)
            and (filter_marker is None or d.marker_type in filter_marker)
            and (filter_size is None or len(d) < filter_size)
            and (filter_qubit is None or d in self.get_diamonds_with_qubit(filter_qubit))
        )

    def mark_error_diamonds_by_duid(self, duids: Iterable[int]):
        new_ds = self
        for duid in duids:
            d = new_ds.get_diamond_with_duid(duid)
            if d is not None:
                new_ds = new_ds.replace_diamond(
                    d, override_color=get_default_color(d.basis, MarkerType.error)
                )
        return new_ds

    def get_boundary_diamonds(self):
        return self.get_diamonds(filter_size=4)

    def has_any_error_diamonds(self):
        return any([d.marker_type == MarkerType.error for d in self.diamonds])

    def plot(
        self,
        ax=None,
        plot_observables=True,
        filter_kwargs=None,
        interpolate_kwargs=None,
        plot_expanding_marker=True,
    ):
        if ax is None:
            _, ax = plt.subplots()

        filter_kwargs = filter_kwargs or {}
        for diamond in self.get_diamonds(**filter_kwargs):
            if interpolate_kwargs is None:
                coords, center_coord = diamond.poly_coords()
            else:
                coords, center_coord = diamond.interpolate_coords_through_cnot_layer(
                    **interpolate_kwargs
                )

            z = 15 - len(diamond.qubits)  # order them with largest polys at the back
            if diamond.marker_type is MarkerType.error:
                z += 10

            patch = mpl.patches.Polygon(
                coords, facecolor=diamond.color, edgecolor='k', linewidth=2, zorder=z
            )
            ax.add_patch(patch)

            # plot an indicator for markertype
            # do nothing for contracting, or for errors and observables
            if diamond.marker_type == MarkerType.expanding and plot_expanding_marker:
                patch = mpl.patches.Circle(
                    xy=center_coord,
                    radius=DIAMOND_CIRCLE_RADIUS,
                    facecolor=diamond.color,
                    edgecolor='k',
                    linewidth=2,
                    zorder=z,
                )
                ax.add_patch(patch)

        z_patch_legend = mpl.patches.Patch(facecolor=get_default_color(Basis.z))
        x_patch_legend = mpl.patches.Patch(facecolor=get_default_color(Basis.x))
        legend_handles = [z_patch_legend, x_patch_legend]
        legend_labels = ['Z Diamond', 'X Diamond']

        if self.has_any_error_diamonds():
            z_error_patch_legend = mpl.patches.Patch(
                facecolor=get_default_color(Basis.z, marker_type=MarkerType.error)
            )
            x_error_patch_legend = mpl.patches.Patch(
                facecolor=get_default_color(Basis.x, marker_type=MarkerType.error)
            )
            legend_handles.extend([z_error_patch_legend, x_error_patch_legend])
            legend_labels.extend(['Z Error', 'X Error'])

        if plot_observables:
            for basis, obv in self.observables.items():
                coords = line_ordered_qubit_coords(obv, flip=(basis == Basis.z))
                ax.plot(
                    coords[:, 0],
                    coords[:, 1],
                    marker=get_observable_marker(basis),
                    markersize=20,
                    linewidth=5,
                    color=get_default_color(basis, marker_type=MarkerType.observable),
                    zorder=15,
                )
            z_obvs = mpl.lines.Line2D(
                [],
                [],
                marker=get_observable_marker(Basis.z),
                color=get_default_color(Basis.z, marker_type=MarkerType.observable),
            )
            x_obvs = mpl.lines.Line2D(
                [],
                [],
                marker=get_observable_marker(Basis.x),
                color=get_default_color(Basis.x, marker_type=MarkerType.observable),
            )
            legend_handles.extend([z_obvs, x_obvs])
            legend_labels.extend(['Z Observable', 'X Observable'])

        exp = mpl.lines.Line2D([], [], marker='o', color='k', fillstyle='none', linestyle='none')
        legend_handles.append(exp)
        legend_labels.append("Expanding")

        ax.legend(legend_handles, legend_labels)

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
