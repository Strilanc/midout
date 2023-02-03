import functools
from typing import Iterable, Optional, FrozenSet, Callable


class Tile:
    """A stabilizer with additional information related to how it is measured.

    Annotates the order in which data qubits are touched, the relevant basis of
    each data qubit, and also the measurement ancilla.
    """

    def __init__(self,
                 *,
                 bases: str,
                 measurement_qubit: complex,
                 ordered_data_qubits: Iterable[Optional[complex]]):
        """
        Args:
            bases: Basis of the stabilizer. A string of XYZ characters the same
                length as the ordered_data_qubits argument. It is permitted to
                give a single-character string, which will automatically be
                expanded to the full length. For example, "X" will become "XXXX"
                if there are four data qubits.
            measurement_qubit: The ancilla qubit used to measure the stabilizer.
            ordered_data_qubits: The data qubits in the stabilizer, in the order
                that they are interacted with. Some entries may be None,
                indicating that no data qubit is interacted with during the
                corresponding interaction layer.
        """
        self.ordered_data_qubits = tuple(ordered_data_qubits)
        self.measurement_qubit = measurement_qubit
        if len(bases) == 1:
            bases *= len(self.ordered_data_qubits)
        self.bases: str = bases
        if len(self.bases) != len(self.ordered_data_qubits):
            raise ValueError('len(self.bases_2) != len(self.data_qubits_order)')

    def __eq__(self, other):
        if not isinstance(other, Tile):
            return False
        return self.ordered_data_qubits == other.ordered_data_qubits and self.measurement_qubit == other.measurement_qubit and self.bases == other.bases

    def __ne__(self, other):
        return not (self == other)

    def __hash__(self):
        return hash((Tile, self.ordered_data_qubits, self.measurement_qubit, self.bases))

    def __repr__(self):
        b = self.basis or self.bases
        return f"""gen.Tile(
    ordered_data_qubits={self.ordered_data_qubits!r},
    measurement_qubit={self.measurement_qubit!r},
    bases={b!r},
)"""

    def after_coordinate_transform(self, coord_transform: Callable[[complex], complex]) -> 'Tile':
        return Tile(
            bases=self.bases,
            ordered_data_qubits=[None if d is None else coord_transform(d) for d in self.ordered_data_qubits],
            measurement_qubit=coord_transform(self.measurement_qubit),
        )

    def after_basis_transform(self, basis_transform: Callable[[str], str]) -> 'Tile':
        return Tile(
            bases=''.join(basis_transform(e) for e in self.bases),
            ordered_data_qubits=self.ordered_data_qubits,
            measurement_qubit=self.measurement_qubit,
        )

    @functools.cached_property
    def data_set(self) -> FrozenSet[complex]:
        return frozenset(e for e in self.ordered_data_qubits if e is not None)

    @functools.cached_property
    def used_set(self) -> FrozenSet[complex]:
        return self.data_set | frozenset([self.measurement_qubit])

    @property
    def basis(self) -> Optional[str]:
        if self.bases == self.bases[0] * len(self.bases):
            return self.bases[0]
        return None
