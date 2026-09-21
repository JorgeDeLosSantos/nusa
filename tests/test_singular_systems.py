"""Regression tests for singular stiffness-system handling."""

import numpy as np
import pytest

from nusa.core import Node
from nusa.element import Bar, LinearTriangle
from nusa.model import BarModel, LinearTriangleModel


def test_underconstrained_bar_raises_clear_singularity_error():
    model = BarModel("Singular bar")
    n1 = Node((0.0, 0.0))
    n2 = Node((1.0, 0.0))
    model.add_nodes([n1, n2])
    model.add_element(Bar((n1, n2), E=100.0, A=1.0))
    model.add_force(n2, (10.0,))

    with pytest.raises(
        np.linalg.LinAlgError,
        match="underconstrained|mechanism",
    ):
        model.solve()


def test_unconstrained_linear_triangle_raises_clear_singularity_error():
    model = LinearTriangleModel("Singular CST")
    n1 = Node((0.0, 0.0))
    n2 = Node((1.0, 0.0))
    n3 = Node((0.0, 1.0))
    model.add_nodes([n1, n2, n3])
    model.add_element(
        LinearTriangle((n1, n2, n3), E=1000.0, nu=0.25, t=0.5)
    )
    model.add_force(n2, (1.0, 0.0))

    with pytest.raises(
        np.linalg.LinAlgError,
        match="underconstrained|mechanism",
    ):
        model.solve()
