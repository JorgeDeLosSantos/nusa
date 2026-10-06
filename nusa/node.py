"""Finite-element node entity."""

import numpy as np


class Node:
    """Geometry/topology point used by finite-element models."""

    __slots__ = ("coordinates", "_label")

    def __init__(self, coordinates):
        try:
            coordinates = np.asarray(coordinates, dtype=float)
        except (TypeError, ValueError) as exc:
            raise ValueError(
                "Node coordinates must contain two finite numbers"
            ) from exc
        if coordinates.shape != (2,):
            raise ValueError("Node coordinates must contain exactly two values")
        if not np.isfinite(coordinates).all():
            raise ValueError("Node coordinates must be finite")

        self.coordinates = coordinates.copy()
        self._label = None

    @property
    def x(self):
        return self.coordinates[0]

    @property
    def y(self):
        return self.coordinates[1]

    @property
    def label(self):
        return self._label

    @label.setter
    def label(self, value):
        self._label = value

    def __str__(self):
        return f"Node {self.label}: ({self.x},{self.y})"

    def __repr__(self):
        return f"<Node {self.label}: ({self.x},{self.y})>"
