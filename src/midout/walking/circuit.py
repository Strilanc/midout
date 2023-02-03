import dataclasses
from typing import Iterable, List, Optional, Union

import matplotlib.pyplot as plt

from midout.walking.cycle import Cycle
from midout.walking.diamonds import Diamond
from midout.walking.meas_layer import MeasLayer
from midout.walking.reset_layer import ResetLayer
from midout.walking.stim_circuit_builder import StimCircuitBuilder
from midout.walking.tiles import TileState
from midout.walking.util import Basis, Direction


@dataclasses.dataclass
class Circuit:
    """represents the entire surface code circuit."""

    # given a sequence of directions and an initial set of tiles,
    # build out all the cycles and the init and final meas rounds

    def __init__(
        self, distance, rounds: Union[int, List[Optional[Direction]]], basis: Basis, duids=False
    ):
        if duids:
            Diamond.reset_duid_counter()

        self.basis = basis
        self.init_tiles = TileState.make_surface_code(distance=distance)
        self.init_diamonds, self.init_layer = ResetLayer.make_init_layer(
            tile_state=self.init_tiles, init_basis=basis
        )

        if isinstance(rounds, int):
            directions: List[Optional[Direction]] = [None] * rounds
        else:
            directions: List[Optional[Direction]] = rounds

        if len(directions) == 0:
            self.cycles: List[Cycle] = []
            self.final_meas_layer = MeasLayer.make_final_measure_layer(
                tile_state=self.init_tiles,
                diamond_state=self.init_diamonds,
                index=0,
                basis=self.basis,
            )
            return

        c0 = Cycle(
            index=0,
            input_tile_state=self.init_tiles,
            input_contracting_diamond_state=self.init_diamonds,
            direction=directions[0],
        )
        self.cycles = [c0]
        for d in directions[1:]:
            c = self.cycles[-1].make_next_cycle(d)
            self.cycles.append(c)

        self.final_meas_layer = MeasLayer.make_final_measure_layer(
            tile_state=self.cycles[-1].output_tile_state,
            diamond_state=self.cycles[-1].diamond_states[-1],
            index=self.cycles[-1].index + 1,
            basis=self.basis,
        )
        if duids:
            Diamond.delete_duid_counter()
        return

    def plot(self, axes=None, plot_gates=False):
        """plot out all the diamond states in this circuit."""
        num = 1 + 6 * len(self.cycles)  # number of diamonds to plot
        if plot_gates:
            num += 1  # if we're plotting gates as well, plot the final measurements at the end
        if axes is None:
            _, axes = plt.subplots(num, figsize=(15, 15 * num))
        self.init_diamonds.plot(ax=axes[0])
        if plot_gates:
            self.init_layer.plot(ax=axes[0])
        axes[0].set_title(
            "Initialization Resets"
            + (f": {self.basis.value} Basis" if self.basis is not None else "")
        )

        for i, c in enumerate(self.cycles):
            c.plot(axes=axes[1 + i * 6 : 1 + (i + 1) * 6], plot_gates=plot_gates)
            for ax in axes[1 + i * 6 : 1 + (i + 1) * 6]:
                ax.set_title(f"Cycle {i}: " + ax.get_title())

        if plot_gates:
            self.final_meas_layer.plot(ax=axes[-1])
            axes[-1].set_title(
                "Final Measurements"
                + (f": {self.basis.value} Basis" if self.basis is not None else "")
            )
        for i, ax in enumerate(axes):
            ax.set_title(ax.get_title() + f" : AFTER TICK {i}")
        return axes

    def plot_with_and_without_gates(self):
        num = 2 + 6 * len(self.cycles)
        _, axes = plt.subplots(num, 2, figsize=(15 * 2, 15 * num))
        self.plot(axes=axes[:, 0], plot_gates=False)
        self.plot(axes=axes[:, 1], plot_gates=True)
        return axes

    @property
    def layers(self):
        return (
            [self.init_layer] + [l for c in self.cycles for l in c.layers] + [self.final_meas_layer]
        )

    def build_stim_circuit(self, *, errors=None):
        builder = StimCircuitBuilder(errors=errors)
        for l in self.layers:
            builder.process_layer(l)
        return builder.stim_circuit

    def mark_error_diamonds_by_duid(self, duids: Iterable[int]):
        self.init_diamonds = self.init_diamonds.mark_error_diamonds_by_duid(duids)
        for c in self.cycles:
            c.mark_error_diamonds_by_duid(duids)
