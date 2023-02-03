import dataclasses
import itertools
from typing import cast, Dict, FrozenSet, List, Set

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

from midout.walking.reset_layer import ResetLayer
from midout.walking.util import (
    Basis,
    get_default_color,
    get_observable_marker,
    line_ordered_qubit_coords,
    MarkerType,
    Qubit,
    qubit_coords,
    rot_order_qubit_coords,
)


@dataclasses.dataclass(frozen=True)
class Tile:
    """A tile for defining the state at measurements of a normal surface code.

    More strictly, tiles here define detecting regions (DRs) at the Measure Reset layer,

    All tiles feature a measure qubit,
    and the corresponding DR that starts from the reset on that qbit in this round

    Tiles can also have data qubits, and if so also represent a fully expanded det region
    that includes the measurement and reset on the measure qubit,
    and all the data qubits

    Attributes:
        data_qubits: a set defining the data qubits for this stabilizer
        measure_qubit: a single measurement ancilla qubit used for measuring this stabilizer
        basis: the basis of the stabilizer this tile is representing
    """

    data_qubits: FrozenSet[Qubit]
    measure_qubit: Qubit
    basis: Basis

    @property
    def all_qubits(self) -> FrozenSet[Qubit]:
        return frozenset(self.data_qubits | {self.measure_qubit})

    def __len__(self):
        return len(self.all_qubits)

    def has_data_qubits(self):
        return len(self.data_qubits) != 0

    def translate(self, direction):
        return Tile(
            data_qubits=frozenset([q + direction for q in self.data_qubits]),
            measure_qubit=self.measure_qubit + direction,
            basis=self.basis,
        )


@dataclasses.dataclass()
class TileState:
    tiles: FrozenSet[Tile]
    top_boundary_basis: Basis
    observables: Dict[Basis, FrozenSet[Qubit]]

    def tiles_with_qubit(self, q: Qubit):
        return set([t for t in self.tiles if q in t.all_qubits])

    def get_center_qubit(self):
        """for cnot layers, we often want the 'center' to be the MQ of a Z type tile."""
        for t in self.tiles:
            if t.basis == Basis.z:
                return t.measure_qubit

    def data_qubit_on_boundary(self, q, basis):
        """returns true if this data qubit is on the boundary of the given basis.

        boundary is given by the stabilizer flavor that has leafs checks
        so the z logical observable terminates on the z boundary
        where there is only 1x type check for each data qubit
        """
        return len([t for t in self.tiles_with_qubit(q) if t.basis == basis.flip()]) == 1

    @property
    def measure_qubits(self) -> Set[Qubit]:
        return set([t.measure_qubit for t in self.tiles])

    @property
    def all_qubits(self) -> Set[Qubit]:
        return set(q for t in self.tiles for q in t.all_qubits)

    @property
    def all_bulk_qubits(self):
        return set(q for t in self.tiles for q in t.all_qubits if len(t.data_qubits) == 4)

    @staticmethod
    def make_surface_code(
        distance,
        flip_stabilizers: bool = False,
        top_boundary_basis=Basis.z,
        observables_at_bottom=False,
    ):
        tiles: Set[Tile] = set()
        h_observable: Set[Qubit] = set()
        v_observable: Set[Qubit] = set()

        for x, y in itertools.product(range(distance + 1), repeat=2):
            # we're going to go through each 'unit cell'
            # which has a mq on the origin coordinate of the unit cell,
            # and decide whether to include it as a tile
            # and which data qubits are associated with the tile
            mq: Qubit = x + 1j * y

            # which kind of mq we have is just set by the parity
            # if flip_stabilizers is false, main diagonal stabilizers, (0,0), (1,1) etc is Z type
            is_x_stabilizer = (x % 2 != y % 2) != flip_stabilizers
            basis: Basis = Basis.x if is_x_stabilizer else Basis.z
            # whether we include the mq is basically yes unless it's on the 'wrong' boundary
            if (x, y) in [(0, 0), (0, distance), (distance, 0), (distance, distance)]:
                # it's the intersection of both boundaries, off the patch, ignore it
                continue
            elif y == 0:  # top boundary
                dqs = cast(List[Qubit], [mq + (+0.5 + 0.5j), mq + (-0.5 + 0.5j)])
                if observables_at_bottom is False:
                    h_observable.update(dqs)
                if top_boundary_basis != basis:
                    continue
            elif y == distance:  # bottom boundary
                dqs = cast(List[Qubit], [mq + (+0.5 - 0.5j), mq + (-0.5 - 0.5j)])
                if observables_at_bottom is True:
                    h_observable.update(dqs)
                if top_boundary_basis != basis:
                    continue
            elif x == 0:  # left boundary
                dqs = cast(List[Qubit], [mq + (+0.5 + 0.5j), mq + (+0.5 - 0.5j)])
                if observables_at_bottom is False:
                    v_observable.update(dqs)
                if top_boundary_basis == basis:
                    continue
            elif x == distance:  # right boundary
                dqs = cast(List[Qubit], [mq + (-0.5 + 0.5j), mq + (-0.5 - 0.5j)])
                if observables_at_bottom is True:
                    v_observable.update(dqs)
                if top_boundary_basis == basis:
                    continue
            else:  # bulk stabilizer,
                dqs = cast(
                    List[Qubit],
                    [
                        mq + (+0.5 + 0.5j),
                        mq + (+0.5 - 0.5j),
                        mq + (-0.5 + 0.5j),
                        mq + (-0.5 - 0.5j),
                    ],
                )
            tiles.add(Tile(data_qubits=frozenset(dqs), measure_qubit=mq, basis=basis))

        # let's make the observables too
        # we're going to use across the top and down the left hand side
        # if the top boundary is x basis (has x leaf checks),
        # then the X observable is the vertical one
        # because it'll terminate on those leaf checks and not get noticed bc it commutes

        z_observable, x_observable = (
            (h_observable, v_observable)
            if top_boundary_basis == Basis.x
            else (v_observable, h_observable)
        )
        return TileState(
            tiles=frozenset(tiles),
            top_boundary_basis=Basis.x,
            observables={Basis.x: frozenset(x_observable), Basis.z: frozenset(z_observable)},
        )

    def plot(self, ax=None, plot_observables=True, plot_measure_qubits=True):
        if ax is None:
            _, ax = plt.subplots()

        patches = []
        for tile in self.tiles:
            mq_coord = qubit_coords(tile.measure_qubit)
            if tile.has_data_qubits():
                coords = rot_order_qubit_coords(tile.data_qubits)
                if len(tile.data_qubits) == 2:
                    # plot a nice wedge instead of a square
                    center = np.mean(coords, axis=0)
                    center_to_mq = mq_coord - center
                    angle_center_to_mq = (
                        np.arctan2(center_to_mq[1], center_to_mq[0]) / 2 / np.pi * 360
                    )
                    patch = mpl.patches.Wedge(
                        center=center,
                        r=np.sqrt(center_to_mq[0] ** 2 + center_to_mq[1] ** 2),
                        theta1=angle_center_to_mq - 90,
                        theta2=angle_center_to_mq + 90,
                        facecolor=get_default_color(tile.basis),
                        edgecolor='k',
                        linewidth=2,
                        zorder=4,
                    )
                else:
                    patch = mpl.patches.Polygon(
                        coords,
                        facecolor=get_default_color(tile.basis),
                        edgecolor='k',
                        linewidth=2,
                        zorder=3,
                    )
                patches.append(patch)
            # add circle for mq
            if plot_measure_qubits:
                circle = mpl.patches.Circle(
                    mq_coord,
                    radius=0.1,
                    facecolor=get_default_color(tile.basis),
                    edgecolor="k",
                    linewidth=2,
                    zorder=5,
                )
                patches.append(circle)
        x_patch = mpl.patches.Patch(facecolor=get_default_color(Basis.x))
        z_patch = mpl.patches.Patch(facecolor=get_default_color(Basis.z))
        legend_handles = [z_patch, x_patch]
        legend_labels = ['Z Stabilizer', 'X Stabilizer']

        for p in patches:
            ax.add_patch(p)

        if plot_observables:
            for basis, obv in self.observables.items():
                coords = line_ordered_qubit_coords(obv)
                ax.plot(
                    coords[:, 0],
                    coords[:, 1],
                    marker=get_observable_marker(basis),
                    markersize=20,
                    linewidth=5,
                    color=get_default_color(basis, marker_type=MarkerType.observable),
                    zorder=10,
                    label=f"{basis.value} Observable",
                )

            z_obvs = mpl.lines.Line2D(
                [],
                [],
                marker='D',
                color=get_default_color(Basis.z, marker_type=MarkerType.observable),
            )
            x_obvs = mpl.lines.Line2D(
                [],
                [],
                marker='s',
                color=get_default_color(Basis.x, marker_type=MarkerType.observable),
            )
            legend_handles.extend([z_obvs, x_obvs])
            legend_labels.extend(['Z Observable', 'X Observable'])
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

    def make_init_diamonds(self):
        contracting_diamond_state, _ = ResetLayer.make_init_layer(tile_state=self, init_basis=None)
        return contracting_diamond_state

    def make_diamonds(self, expand_direction=None):
        full_diamond_state, _ = ResetLayer.from_tiles_and_diamonds(
            contracting_diamond_state=self.make_init_diamonds(),
            tile_state=self,
            direction=expand_direction,
        )
        return full_diamond_state

    def translate(self, direction) -> 'TileState':
        """returns the equiv tilestate translated by direction."""
        return (
            TileState(
                tiles=frozenset([t.translate(direction) for t in self.tiles]),
                top_boundary_basis=self.top_boundary_basis,
                observables={
                    b: frozenset([q + direction for q in obv])
                    for b, obv in self.observables.items()
                },
            )
            if direction is not None
            else self
        )
