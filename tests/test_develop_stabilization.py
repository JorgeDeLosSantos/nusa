"""Regression tests for the develop-branch API stabilization."""

import numpy as np

from nusa.core import Element, Model, Node
from nusa.element import Beam, Spring
from nusa.model import BeamModel, SpringModel
from nusa.version import __version__


class MockElement(Element):
    def __init__(self, nodes):
        super().__init__("mock")
        self.nodes = nodes
        self.f = 1.0
        self.s = 2.0


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
    assert __version__.startswith("0.3.0.dev")


def test_model_report_helpers_use_property_based_model_api():
    model = Model("Regression model", "mock")
    n1 = Node((0.0, 0.0))
    n2 = Node((1.0, 0.0))
    element = MockElement((n1, n2))

    model.add_nodes([n1, n2])
    model.add_element(element)

    options = {
        "headers": "firstrow",
        "tablefmt": "rst",
        "numalign": "right",
    }

    tables = (
        model._get_ndisplacements(options),
        model._get_nforces(options),
        model._get_nodes_info(options),
        model._get_elements_info(options),
    )

    assert all(isinstance(table, str) for table in tables)
    assert "Node" in tables[0]
    assert "Node" in tables[2]
    assert "Element" in tables[3]


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



def test_shared_solver_uses_vector_state():
    model = SpringModel("Solver state")
    n1 = Node((0.0, 0.0))
    n2 = Node((0.0, 0.0))
    element = Spring((n1, n2), 300.0)

    model.add_nodes([n1, n2])
    model.add_element(element)
    model.add_constraint(n1, ux=0.0)
    model.add_force(n2, (750.0,))
    model.solve()

    np.testing.assert_allclose(model._u, [0.0, 2.5])
    np.testing.assert_allclose(model._f, [0.0, 750.0])
    np.testing.assert_allclose(model._nodal_forces, [-750.0, 750.0])
    assert model._prescribed_dofs == [0]
    assert model._free_dofs == [1]
    np.testing.assert_allclose(model._K_reduced, [[300.0]])
    np.testing.assert_allclose(model._rhs_reduced, [750.0])

    for legacy_name in (
        "U",
        "F",
        "NF",
        "VU",
        "VF",
        "K2S",
        "F2S",
        "solved_u",
    ):
        assert not hasattr(model, legacy_name)



def test_reduced_rhs_includes_nonzero_prescribed_displacements():
    model = SpringModel("Reduced RHS")
    n1 = Node((0.0, 0.0))
    n2 = Node((0.0, 0.0))
    n3 = Node((0.0, 0.0))

    model.add_nodes([n1, n2, n3])
    model.add_elements([
        Spring((n1, n2), 100.0),
        Spring((n2, n3), 100.0),
    ])
    model.add_constraint(n1, ux=0.0)
    model.add_constraint(n3, ux=0.03)
    model.solve()

    assert model._prescribed_dofs == [0, 2]
    assert model._free_dofs == [1]
    np.testing.assert_allclose(model._K_reduced, [[200.0]])
    np.testing.assert_allclose(model._rhs_reduced, [3.0])
    assert np.isclose(n2.ux, 0.015)


def test_shared_assembly_accumulates_overlapping_beam_dofs():
    model = BeamModel("Assembly regression")
    n1 = Node((0.0, 0.0))
    n2 = Node((1.0, 0.0))
    n3 = Node((2.0, 0.0))
    e1 = Beam((n1, n2), E=1.0, I=1.0)
    e2 = Beam((n2, n3), E=1.0, I=1.0)

    model.add_nodes([n1, n2, n3])
    model.add_elements([e1, e2])
    model.assemble()

    expected = np.array(
        [
            [12.0, 6.0, -12.0, 6.0, 0.0, 0.0],
            [6.0, 4.0, -6.0, 2.0, 0.0, 0.0],
            [-12.0, -6.0, 24.0, 0.0, -12.0, 6.0],
            [6.0, 2.0, 0.0, 8.0, -6.0, 2.0],
            [0.0, 0.0, -12.0, -6.0, 12.0, -6.0],
            [0.0, 0.0, 6.0, 2.0, -6.0, 4.0],
        ]
    )

    np.testing.assert_allclose(model.stiffness_matrix, expected)
    assert model._is_assembled is True
