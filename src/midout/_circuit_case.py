import dataclasses
from typing import List, FrozenSet

import stim

from midout import gen


@dataclasses.dataclass(frozen=True)
class CircuitCase:
    circuit: stim.Circuit
    patches: List[gen.Patch]
    expected_interactions: FrozenSet[str]
    show_patch_order: bool = False
    show_patch_measure_qubits: bool = False
