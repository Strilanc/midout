from typing import Iterable

import matplotlib.pyplot as plt

from midout.walking.cnot_layer import CnotLayer
from midout.walking.diamonds import DiamondState
from midout.walking.meas_layer import MeasLayer
from midout.walking.reset_layer import ResetLayer
from midout.walking.tiles import TileState
from midout.walking.util import DOWN_RIGHT, get_compass_directions, UP_LEFT


class Cycle:
    """a full cycle for the surface code.

    on construction, takes a direction and a tile set,

    it then constructs the reset, cnot and measurement layers,
    and the set of tiles that should be used in the next cycle
    """

    def __init__(
        self,
        index,
        input_tile_state: TileState,
        input_contracting_diamond_state: DiamondState,
        direction=None,
        dont_make_measurements=False,
        flip_left_right=False,
        reverse_expansion=False,
    ):
        self.index = index
        self.input_tile_state = input_tile_state
        self.output_tile_state = input_tile_state.translate(direction)
        self._dont_make_measurements = dont_make_measurements

        self.post_reset_diamonds, self.reset_layer = ResetLayer.from_tiles_and_diamonds(
            input_tile_state,
            input_contracting_diamond_state,
            direction=direction,
            reverse_expansion=reverse_expansion,
        )
        # origin for the cnot layers should be the MQ of a Z tile, any Z tile will do
        center_qubit = input_tile_state.get_center_qubit()
        if direction is None:
            self.cnot_layers = list(
                CnotLayer.for_standard_surface_code_cycle(
                    qubits=self.post_reset_diamonds.all_qubits,
                    center=center_qubit,
                    arrangement='standard',
                )
            )
        else:
            just_original_qubits = self.input_tile_state.all_qubits
            just_new_qubits = self.output_tile_state.all_qubits
            original_and_new_bulk_qubits = (
                just_original_qubits | self.output_tile_state.all_bulk_qubits
            )
            original_bulk_and_new_qubits = self.input_tile_state.all_bulk_qubits | just_new_qubits
            left_then_right = (direction not in [DOWN_RIGHT, UP_LEFT]) != flip_left_right
            compass_directions = get_compass_directions(direction, left_then_right=left_then_right)
            if not reverse_expansion:  # normal expansion
                self.cnot_layers = [
                    CnotLayer.from_layer_spec(
                        just_original_qubits,
                        center_qubit,
                        direction=compass_directions[0],
                        aligned=True,
                    ),
                    CnotLayer.from_layer_spec(
                        just_original_qubits,
                        center_qubit,
                        direction=compass_directions[1],
                        aligned=False,
                    ),
                    CnotLayer.from_layer_spec(
                        original_and_new_bulk_qubits,
                        center_qubit + direction,
                        direction=compass_directions[2],
                        aligned=False,
                    ),
                    CnotLayer.from_layer_spec(
                        original_and_new_bulk_qubits,
                        center_qubit + direction,
                        direction=compass_directions[0],
                        aligned=True,
                    ),
                ]
            else:  # reverse expansion
                self.cnot_layers = [
                    CnotLayer.from_layer_spec(
                        original_bulk_and_new_qubits,
                        center_qubit,
                        direction=compass_directions[3],
                        aligned=True,
                    ),
                    CnotLayer.from_layer_spec(
                        original_bulk_and_new_qubits,
                        center_qubit,
                        direction=compass_directions[1],
                        aligned=False,
                    ),
                    CnotLayer.from_layer_spec(
                        just_new_qubits,
                        center_qubit + direction,
                        direction=compass_directions[2],
                        aligned=False,
                    ),
                    CnotLayer.from_layer_spec(
                        just_new_qubits,
                        center_qubit + direction,
                        direction=compass_directions[3],
                        aligned=True,
                    ),
                ]

        self.measure_layer = None
        self.post_measure_cont_diamonds = None
        self.layers = []
        self.diamond_states = []
        self.build_diamonds_states_and_layers()

    def build_diamonds_states_and_layers(self):
        while self.build_next_layer():
            pass
        return

    def build_next_layer(self):
        n = len(self.diamond_states)
        if n == 0:
            self.diamond_states = [self.post_reset_diamonds]
            self.layers = [self.reset_layer]
        elif n <= 4:
            cnot_layer = self.cnot_layers[n - 1]
            self.diamond_states.append(
                cnot_layer.track_diamond_state(diamond_state=self.diamond_states[-1])
            )
            self.layers.append(cnot_layer)
        elif n == 5 and not self._dont_make_measurements:
            # _dont_make_measurements is helpful for debugging the cnot layers,
            # because the measure layer is where we start making consistency checks
            # eg if contracting diamonds aren't all measured
            # or if expanding diamonds cover multiple measure qubits
            self.post_measure_cont_diamonds, self.measure_layer = MeasLayer.make_measure_layer(
                self.output_tile_state, self.diamond_states[-1], index=self.index
            )
            self.diamond_states.append(self.post_measure_cont_diamonds)
            self.layers.append(self.measure_layer)
        else:
            return False
        return True

    def make_next_cycle(self, direction):
        if self._dont_make_measurements:
            raise ValueError("can't make next cycle correctly if we didn't make measurements")
        return Cycle(
            index=self.index + 1,
            input_tile_state=self.output_tile_state,
            input_contracting_diamond_state=self.diamond_states[-1],
            direction=direction,
        )

    def plot(self, axes=None, plot_gates=False):
        """plot out the 6 diamonds from this cycle, optionally with cnots."""
        if axes is None:
            p = len(self.diamond_states)
            _, axes = plt.subplots(p, figsize=(15, 15 * p))
        titles = ['Resets'] + [f'CNOT layer {i+1}' for i in range(4)] + ['Measurements']
        for i, (ax, layer, d) in enumerate(zip(axes.flatten(), self.layers, self.diamond_states)):
            d.plot(ax=ax, plot_observables=True, filter_kwargs=dict(filter_marker=None))
            if layer is not None and plot_gates:
                layer.plot(ax=ax)
                ax.get_legend().remove()
            ax.set_title(titles[i])

    def mark_error_diamonds_by_duid(self, duids: Iterable[int]):
        self.diamond_states = [ds.mark_error_diamonds_by_duid(duids) for ds in self.diamond_states]
