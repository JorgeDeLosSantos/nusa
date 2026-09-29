"""Regression tests retained from the 0.3 stabilization phase."""

import re
import numpy as np

from nusa import BeamModel, Node, Spring, SpringModel
from nusa.core import Element, Model
from nusa.version import __version__


class MockBeam(Element):
    def __init__(self, nodes):
        super().__init__("beam")
        self.nodes = tuple(nodes)

    def get_element_stiffness(self):
        return np.array([
            [12.0, 6.0, -12.0, 6.0],
            [6.0, 4.0, -6.0, 2.0],
            [-12.0, -6.0, 12.0, -6.0],
            [6.0, 2.0, -6.0, 4.0],
        ])

    def compute_results(self, u_e):
        actions = self.get_element_stiffness() @ np.asarray(u_e, dtype=float)
        return {
            "shear_force_i": float(actions[0]),
            "shear_force_j": float(actions[2]),
            "bending_moment_i": float(actions[1]),
            "bending_moment_j": float(actions[3]),
        }


def test_version_belongs_to_0_3_0_release_line_until_version_bump():
    assert re.fullmatch(
        r"0\.3\.0(?:\.dev\d+|(?:a|b|rc)\d+(?:\.dev\d+)?)?",
        __version__,
    )


def test_report_helpers_are_not_model_responsibilities():
    model = Model("Regression model", "mock")
    for name in (
        "_get_ndisplacements", "_get_nforces", "_get_nodes_info",
        "_get_elements_info", "_get_element_results",
    ):
        assert not hasattr(model, name)


def test_mock_beam_solves_through_result_contract():
    model = BeamModel("Regression beam")
    n1, n2 = Node((0, 0)), Node((1, 0))
    model.add_nodes([n1, n2])
    model.add_element(MockBeam((n1, n2)))
    model.add_constraint(n1, uy=0.0, ur=0.0)
    model.add_force(n2, (-1.0,))

    result = model.solve()

    np.testing.assert_allclose(result.displacements, [0.0, 0.0, -1.0 / 3.0, -0.5])
    np.testing.assert_allclose(result.reactions[:2], [1.0, 1.0])


def test_nonzero_prescribed_displacement_reduced_system_behavior():
    model = SpringModel("Reduced RHS")
    n1, n2, n3 = Node((0, 0)), Node((0, 0)), Node((0, 0))
    model.add_nodes([n1, n2, n3])
    model.add_elements([Spring((n1, n2), 100.0), Spring((n2, n3), 100.0)])
    model.add_constraint(n1, ux=0.0)
    model.add_constraint(n3, ux=0.03)

    result = model.solve()

    np.testing.assert_allclose(result.displacements, [0.0, 0.015, 0.03])
    assert not hasattr(model, "_K_reduced")
