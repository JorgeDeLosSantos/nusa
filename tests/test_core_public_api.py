"""Regression tests for the 0.4 core public contract."""

import numpy as np
import pytest

from nusa import Element, Node
from nusa.element import Beam
from nusa.model import BeamModel, LinearTriangleModel, TrussModel


def test_node_and_element_are_solution_state_free():
    node = Node((0.0, 0.0))
    element = Element("mock")

    for name in (
        "ux", "uy", "ur", "fx", "fy", "m",
        "sx", "sy", "sxy", "seqv", "ex", "ey", "exy",
        "get_displacements", "set_displacements", "get_forces", "set_forces",
    ):
        assert not hasattr(node, name)

    for name in (
        "fx", "fy", "sx", "sy", "sxy",
        "set_element_forces", "get_element_forces",
    ):
        assert not hasattr(element, name)


def test_constraints_are_stored_only_on_model():
    model = TrussModel("constraints")
    node = Node((0.0, 0.0))
    model.add_node(node)

    model.add_constraint(node, ux=0.0)
    model.add_constraint(node, uy=0.25)

    assert model._prescribed_displacements[node] == {"ux": 0.0, "uy": 0.25}
    assert model.prescribed_displacement(node) == {"ux": 0.0, "uy": 0.25}
    assert not hasattr(node, "ux")


def test_linear_triangle_accepts_independent_constraints():
    model = LinearTriangleModel("constraints")
    node = Node((0.0, 0.0))
    model.add_node(node)

    model.add_constraint(node, ux=0.1)
    model.add_constraint(node, uy=-0.2)

    assert model.prescribed_displacement(node) == {"ux": 0.1, "uy": -0.2}


def test_beam_constraints_use_only_active_dofs():
    model = BeamModel("beam constraints")
    n1 = Node((0.0, 0.0))
    n2 = Node((1.0, 0.0))
    model.add_nodes([n1, n2])
    model.add_element(Beam((n1, n2), E=1.0, I=1.0))

    with pytest.raises(ValueError, match="Unsupported constraint"):
        model.add_constraint(n1, ux=0.0)

    model.add_constraint(n1, uy=0.0)
    model.add_constraint(n1, ur=0.0)
    assert model.prescribed_displacement(n1) == {"uy": 0.0, "ur": 0.0}

    free = model.prescribed_displacement(n2)
    assert np.isnan(free["uy"])
    assert np.isnan(free["ur"])
