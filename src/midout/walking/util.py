import dataclasses
from enum import Enum
from typing import cast, NewType, Optional

import numpy as np

# I'm aware this is bad, but there you go, I'm lazy
Qubit = NewType("Qubit", complex)
Direction = NewType("Direction", complex)

DOWN_RIGHT: Direction = cast(Direction, 0.5 + 0.5j)
DOWN_LEFT: Direction = cast(Direction, -0.5 + 0.5j)
UP_RIGHT: Direction = cast(Direction, 0.5 - 0.5j)
UP_LEFT: Direction = cast(Direction, -0.5 - 0.5j)


def get_compass_directions(direction, left_then_right=True):
    """returns the original direction, perpendicular directions, and the opposite direction.

    Returns:
        direction: the input direction
        left_direction: the direction 90 to the left when facing in direction
        right_direction:  the direction 90 to the right when facing in direction
        anti_direction: the opposite direction to direction
    """
    left_right_correction = 1 if left_then_right else -1
    return [
        direction,
        (direction.imag - direction.real * 1j) * left_right_correction,
        (-direction.imag + direction.real * 1j) * left_right_correction,
        -1 * direction,
    ]


def qubit_coords(qubit):
    return (qubit.real, qubit.imag)


def qubit_coords_array(qubits):
    return np.array([[q.real, q.imag] for q in qubits])


def rot_order_qubit_coords(qubits):
    """return the data qubits in 'rotation order', convenient for drawing closed polygons."""
    qubits = np.array([q for q in qubits])
    center = np.mean(qubits)
    center_angles = np.angle(qubits - center)
    return [
        [x.real, x.imag] for _, x in sorted(zip(center_angles, qubits), key=lambda pair: pair[0])
    ]


def line_ordered_qubit_coords(qubits, flip=False):
    if flip:
        qubits = [q.imag + 1j * q.real for q in qubits]
    coords = qubit_coords_array(np.sort([q for q in qubits]))
    if flip:
        coords = np.flip(coords, axis=1)
    return coords


def shoelace_area(points):
    """return the signed shoelace area
    assumes points are around the perimeter of a non-overlapping polygon
    """
    return 0.5 * np.sum(
        [
            points[i - 1][0] * points[i][1] - points[i][0] * points[i - 1][1]
            for i in range(len(points))
        ]
    )


def centroid(points):
    """returns the centroid of a 2D polygon.

    assumes the points are ordered around the perimeter of the shape,
    and that the polygon is not self-overlapping
    """
    # so, this is a bit batshit, but:
    # https://en.wikipedia.org/wiki/Centroid#Of_a_polygon
    if len(points) == 1:
        return points[0]
    A = shoelace_area(points)
    if np.isclose(A, 0):
        # the area is 0 if there are 2 points, or the points are colinear
        return np.mean(points, axis=0)
    lx = [
        (points[i - 1][0] + points[i][0])
        * (points[i - 1][0] * points[i][1] - points[i][0] * points[i - 1][1])
        for i in range(len(points))
    ]
    cx = 1 / 6.0 / A * np.sum(lx)
    ly = [
        (points[i - 1][1] + points[i][1])
        * (points[i - 1][0] * points[i][1] - points[i][0] * points[i - 1][1])
        for i in range(len(points))
    ]
    cy = 1 / 6.0 / A * np.sum(ly)
    return [cx, cy]


class Basis(Enum):
    x = 'x'
    z = 'z'

    def flip(self):
        return {Basis.x: Basis.z, Basis.z: Basis.x}[self]


class MarkerType(Enum):
    expanding = 'expanding'
    contracting = 'contracting'
    error = 'error'
    observable = 'observable'


def get_default_color(basis: Optional[Basis], marker_type: MarkerType = MarkerType.expanding):
    if basis is None:
        return 'grey'
    colors = {
        MarkerType.expanding: {Basis.x: 'lightcoral', Basis.z: 'skyblue'},
        MarkerType.contracting: {Basis.x: 'lightcoral', Basis.z: 'skyblue'},
        MarkerType.error: {Basis.x: 'yellow', Basis.z: 'springgreen'},
        MarkerType.observable: {Basis.x: 'crimson', Basis.z: 'darkblue'},
    }
    return colors[marker_type][basis]


def get_observable_marker(basis: Basis):
    return 'X' if basis == Basis.x else 'P'


@dataclasses.dataclass(frozen=True)
class Measurement:
    index: int
    qubit: Qubit
