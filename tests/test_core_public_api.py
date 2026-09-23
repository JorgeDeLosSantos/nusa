"""Regression tests for the cleaned core public API."""

import numpy as np

from nusa.core import Element, Node
from nusa.element import Beam, LinearTriangle, Spring
from nusa.model import BeamModel, LinearTriangleModel, TrussModel


def test_legacy_node_and_element_getters_setters_are_removed():
    node = Node((0.0, 0.0))
    element = Element("mock")

    for name in (
        "get_label",
        "set_label",
        "get_displacements",
        "set_displacements",
        "get_forces",
        "set_forces",
    ):
        assert not hasattr(node, name)

    for name in (
        "set_label",
        "set_element_forces",
        "get_element_forces",
        "get_nodes",
    ):
        assert not hasattr(element, name)


def test_truss_constraints_can_be_added_successively_without_resetting_components():
    model = TrussModel("Successive constraints")
    node = Node((0.0, 0.0))
    model.add_node(node)

    model.add_constraint(node, ux=0.0)
    assert np.isclose(node.ux, 0.0)
    assert np.isnan(node.uy)

    model.add_constraint(node, uy=0.25)
    assert np.isclose(node.ux, 0.0)
    assert np.isclose(node.uy, 0.25)
    assert model._prescribed_displacements[node] == {"ux": 0.0, "uy": 0.25}


def test_linear_triangle_accepts_independent_displacement_constraints():
    model = LinearTriangleModel("Independent CST constraints")
    node = Node((0.0, 0.0))
    model.add_node(node)

    model.add_constraint(node, ux=0.1)
    assert np.isclose(node.ux, 0.1)
    assert np.isnan(node.uy)

    model.add_constraint(node, uy=-0.2)
    assert np.isclose(node.ux, 0.1)
    assert np.isclose(node.uy, -0.2)
    assert model._prescribed_displacements[node] == {"ux": 0.1, "uy": -0.2}


def test_beam_constraints_use_only_active_displacement_dofs():
    model = BeamModel("Beam constraint contract")
    n1 = Node((0.0, 0.0))
    n2 = Node((1.0, 0.0))
    model.add_nodes([n1, n2])
    model.add_element(Beam((n1, n2), E=1.0, I=1.0))

    with np.testing.assert_raises_regex(ValueError, "Unsupported constraint"):
        model.add_constraint(n1, ux=0.0)

    model.add_constraint(n1, uy=0.0)
    model.add_constraint(n1, ur=0.0)

    n3 = Node((2.0, 0.0))
    model.add_node(n3)

    assert np.isclose(n1.uy, 0.0)
    assert np.isclose(n1.ur, 0.0)
    assert model._prescribed_displacements[n1] == {
        "uy": 0.0,
        "ur": 0.0,
    }
