"""Regression tests for the develop-branch API stabilization."""

import numpy as np

from nusa.core import Element, Node
from nusa.model import BeamModel, SpringModel
from nusa.version import __version__


class MockSpring(Element):
    def __init__(self, nodes, stiffness=1.0):
        super().__init__("spring")
        self.nodes = nodes
        self.stiffness = stiffness

    def get_element_stiffness(self):
        k = self.stiffness
        return np.array([[k, -k], [-k, k]], dtype=float)

    @property
    def fx(self):
        n1, n2 = self.nodes
        return self.get_element_stiffness() @ np.array([[n1.ux], [n2.ux]])


class MockBeam(Element):
    def __init__(self, nodes):
        super().__init__("beam")
        self.nodes = nodes

    def get_element_stiffness(self):
        # Euler-Bernoulli beam with E = I = L = 1.
        return np.array(
            [
                [12.0, 6.0, -12.0, 6.0],
                [6.0, 4.0, -6.0, 2.0],
                [-12.0, -6.0, 12.0, -6.0],
                [6.0, 2.0, -6.0, 4.0],
            ]
        )


def test_develop_version_remains_marked_as_development():
    assert __version__ == "0.3.0.dev0"


def test_spring_simple_report_uses_property_based_model_api():
    model = SpringModel("Regression spring")
    n1 = Node((0.0, 0.0))
    n2 = Node((1.0, 0.0))
    element = MockSpring((n1, n2), stiffness=10.0)

    model.add_nodes([n1, n2])
    model.add_element(element)
    model.add_constraint(n1, ux=0.0)
    model.add_force(n2, (10.0,))
    model.solve()

    report = model.simple_report(report_type="string")

    assert "Regression spring" in report
    assert "Number of nodes: 2" in report
    assert "Number of elements: 1" in report


def test_beam_solve_indexes_property_based_node_collection_with_integers():
    model = BeamModel("Regression beam")
    n1 = Node((0.0, 0.0))
    n2 = Node((1.0, 0.0))
    element = MockBeam((n1, n2))

    model.add_nodes([n1, n2])
    model.add_element(element)
    model.add_constraint(n1, ux=0.0, uy=0.0, ur=0.0)
    model.add_force(n2, (-1.0,))

    model.solve()

    assert np.isclose(n2.uy, -1.0 / 3.0)
    assert np.isclose(n2.ur, -0.5)
    assert np.isclose(n1.fy, 1.0)
    assert np.isclose(n1.m, 1.0)
