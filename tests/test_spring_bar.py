"""Numerical regression tests for spring and bar finite elements."""

import numpy as np
from nusa.core import Node
from nusa.element import Bar, Spring
from nusa.model import BarModel, SpringModel


class TestSpringElement:
    def test_element_stiffness_matrix(self):
        n1 = Node((0.0, 0.0))
        n2 = Node((0.0, 0.0))
        element = Spring((n1, n2), 300.0)

        expected = np.array([[300.0, -300.0], [-300.0, 300.0]])

        np.testing.assert_allclose(element.get_element_stiffness(), expected)

    def test_element_forces_from_nodal_displacements(self):
        n1 = Node((0.0, 0.0))
        n2 = Node((0.0, 0.0))
        n1.ux = 0.0
        n2.ux = 2.5
        element = Spring((n1, n2), 300.0)

        expected = np.array([[-750.0], [750.0]])

        np.testing.assert_allclose(element.fx, expected)


class TestSpringModel:
    def test_single_spring_reference_solution(self):
        """Single spring: P=750, k=300 -> u=P/k=2.5."""
        model = SpringModel("Single spring")
        n1 = Node((0.0, 0.0))
        n2 = Node((0.0, 0.0))
        element = Spring((n1, n2), 300.0)

        model.add_nodes([n1, n2])
        model.add_element(element)
        model.add_constraint(n1, ux=0.0)
        model.add_force(n2, (750.0,))
        model.solve()

        assert np.isclose(n1.ux, 0.0)
        assert np.isclose(n2.ux, 2.5)
        assert np.isclose(n1.fx, -750.0)
        assert np.isclose(n2.fx, 750.0)
        np.testing.assert_allclose(element.fx.ravel(), [-750.0, 750.0])

    def test_logan_example_2_1(self):
        """Regression of examples/spring/spring_01.py (Logan, Example 2.1)."""
        model = SpringModel("Logan 2.1")
        n1, n2, n3, n4 = [Node((0.0, 0.0)) for _ in range(4)]
        e1 = Spring((n1, n3), 1000.0)
        e2 = Spring((n3, n4), 2000.0)
        e3 = Spring((n4, n2), 3000.0)

        model.add_nodes([n1, n2, n3, n4])
        model.add_elements([e1, e2, e3])
        model.add_constraint(n1, ux=0.0)
        model.add_constraint(n2, ux=0.0)
        model.add_force(n4, (5000.0,))
        model.solve()

        np.testing.assert_allclose(
            [n3.ux, n4.ux],
            [10.0 / 11.0, 15.0 / 11.0],
        )
        np.testing.assert_allclose(
            [n1.fx, n2.fx],
            [-10000.0 / 11.0, -45000.0 / 11.0],
        )
        np.testing.assert_allclose(
            e1.fx.ravel(),
            [-10000.0 / 11.0, 10000.0 / 11.0],
        )
        np.testing.assert_allclose(
            e2.fx.ravel(),
            [-10000.0 / 11.0, 10000.0 / 11.0],
        )
        np.testing.assert_allclose(
            e3.fx.ravel(),
            [45000.0 / 11.0, -45000.0 / 11.0],
        )

    def test_nonzero_prescribed_displacement(self):
        """Regression target from examples/spring/spring_02.py (Logan, Example 2.2)."""
        model = SpringModel("Logan 2.2")
        nodes = [Node((0.0, 0.0)) for _ in range(5)]
        elements = [
            Spring((nodes[i], nodes[i + 1]), 200e3)
            for i in range(4)
        ]

        model.add_nodes(nodes)
        model.add_elements(elements)
        model.add_constraint(nodes[0], ux=0.0)
        model.add_constraint(nodes[4], ux=0.02)
        model.add_force(nodes[3], (4000.0,))
        model.solve()

        np.testing.assert_allclose(
            [node.ux for node in nodes],
            [0.0, 0.01, 0.02, 0.03, 0.02],
        )
        np.testing.assert_allclose(
            [node.fx for node in nodes],
            [-2000.0, 0.0, 0.0, 4000.0, -2000.0],
            atol=1e-10,
        )


class TestBarElement:
    def test_length_and_stiffness_matrix(self):
        n1 = Node((0.0, 0.0))
        n2 = Node((3.0, 4.0))
        element = Bar((n1, n2), E=200.0, A=2.0)

        assert np.isclose(element.L, 5.0)
        np.testing.assert_allclose(
            element.get_element_stiffness(),
            [[80.0, -80.0], [-80.0, 80.0]],
        )

    def test_element_force_and_stress(self):
        n1 = Node((0.0, 0.0))
        n2 = Node((2.0, 0.0))
        n1.ux = 0.0
        n2.ux = 0.01
        element = Bar((n1, n2), E=200e9, A=0.001)

        np.testing.assert_allclose(
            element.fx.ravel(),
            [-1.0e6, 1.0e6],
        )
        np.testing.assert_allclose(
            element.sx,
            [-1.0e9, 1.0e9],
        )


class TestBarModel:
    def test_logan_example_3_1(self):
        """Regression of examples/bar/bar_1.py (Logan, Example 3.1)."""
        model = BarModel("Logan 3.1")
        n1 = Node((0.0, 0.0))
        n2 = Node((30.0, 0.0))
        n3 = Node((60.0, 0.0))
        n4 = Node((90.0, 0.0))

        e1 = Bar((n1, n2), E=30e6, A=1.0)
        e2 = Bar((n2, n3), E=30e6, A=1.0)
        e3 = Bar((n3, n4), E=15e6, A=2.0)

        model.add_nodes([n1, n2, n3, n4])
        model.add_elements([e1, e2, e3])
        model.add_constraint(n1, ux=0.0)
        model.add_constraint(n4, ux=0.0)
        model.add_force(n2, (3000.0,))
        model.solve()

        np.testing.assert_allclose(
            [n1.ux, n2.ux, n3.ux, n4.ux],
            [0.0, 0.002, 0.001, 0.0],
        )
        np.testing.assert_allclose(
            [n1.fx, n2.fx, n3.fx, n4.fx],
            [-2000.0, 3000.0, 0.0, -1000.0],
            atol=1e-10,
        )
        np.testing.assert_allclose(e1.fx.ravel(), [-2000.0, 2000.0])
        np.testing.assert_allclose(e2.fx.ravel(), [1000.0, -1000.0])
        np.testing.assert_allclose(e3.fx.ravel(), [1000.0, -1000.0])
        np.testing.assert_allclose(e1.sx, [-2000.0, 2000.0])
        np.testing.assert_allclose(e2.sx, [1000.0, -1000.0])
        np.testing.assert_allclose(e3.sx, [500.0, -500.0])

    def test_nonzero_prescribed_displacement_with_single_unknown(self):
        """Regression of examples/bar/bar_2.py (Kattan, Example 3.1)."""
        E = 210e6
        A = 0.003

        model = BarModel("Kattan 3.1")
        n1 = Node((0.0, 0.0))
        n2 = Node((1.5, 0.0))
        n3 = Node((2.5, 0.0))
        e1 = Bar((n1, n2), E, A)
        e2 = Bar((n2, n3), E, A)

        model.add_nodes([n1, n2, n3])
        model.add_elements([e1, e2])
        model.add_constraint(n1, ux=0.0)
        model.add_constraint(n3, ux=0.002)
        model.add_force(n2, (-10.0,))
        model.solve()

        k1 = E * A / 1.5
        k2 = E * A / 1.0
        expected_u2 = (-10.0 + k2 * 0.002) / (k1 + k2)

        assert np.isclose(n2.ux, expected_u2)

        expected_u = np.array([0.0, expected_u2, 0.002])
        expected_k = np.array(
            [
                [k1, -k1, 0.0],
                [-k1, k1 + k2, -k2],
                [0.0, -k2, k2],
            ]
        )
        expected_nodal_forces = expected_k @ expected_u

        np.testing.assert_allclose(
            [n1.fx, n2.fx, n3.fx],
            expected_nodal_forces,
        )


def test_bar_nonzero_prescribed_displacement_with_multiple_unknowns():
    model = BarModel("Prescribed bar chain")
    nodes = [Node((float(i), 0.0)) for i in range(4)]
    elements = [Bar((nodes[i], nodes[i + 1]), E=100.0, A=1.0) for i in range(3)]

    model.add_nodes(nodes)
    model.add_elements(elements)
    model.add_constraint(nodes[0], ux=0.0)
    model.add_constraint(nodes[3], ux=0.03)
    model.solve()

    np.testing.assert_allclose(
        [node.ux for node in nodes],
        [0.0, 0.01, 0.02, 0.03],
        atol=1e-12,
    )
