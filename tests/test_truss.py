"""Numerical regression tests for the 2D truss element and model."""

import numpy as np

from nusa.core import Node
from nusa.element import Truss
from nusa.model import TrussModel


class TestTrussElement:
    def test_geometry_and_global_stiffness_matrix(self):
        n1 = Node((0.0, 0.0))
        n2 = Node((3.0, 4.0))
        element = Truss((n1, n2), E=200.0, A=2.0)

        assert np.isclose(element.L, 5.0)
        assert np.isclose(element.theta, np.arctan2(4.0, 3.0))

        expected = np.array(
            [
                [28.8, 38.4, -28.8, -38.4],
                [38.4, 51.2, -38.4, -51.2],
                [-28.8, -38.4, 28.8, 38.4],
                [-38.4, -51.2, 38.4, 51.2],
            ]
        )
        np.testing.assert_allclose(element.get_element_stiffness(), expected)

    def test_axial_force_and_stress_for_arbitrary_orientation(self):
        n1 = Node((0.0, 0.0))
        n2 = Node((3.0, 4.0))
        n1.ux = 0.0
        n1.uy = 0.0
        # Displacement of 0.01 along the element axis.
        n2.ux = 0.006
        n2.uy = 0.008

        element = Truss((n1, n2), E=200.0, A=2.0)

        assert np.isclose(element.f, 0.8)
        assert np.isclose(element.s, 0.4)

    def test_rigid_body_translation_produces_zero_axial_force(self):
        n1 = Node((0.0, 0.0))
        n2 = Node((3.0, 4.0))
        n1.ux, n1.uy = 0.25, -0.1
        n2.ux, n2.uy = 0.25, -0.1

        element = Truss((n1, n2), E=200.0, A=2.0)

        assert np.isclose(element.f, 0.0, atol=1e-12)
        assert np.isclose(element.s, 0.0, atol=1e-12)


class TestTrussModel:
    def test_three_member_reference_problem(self):
        """Regression of examples/truss/truss_01.py."""
        E = 30e6
        A = 2.0
        P = 10e3

        model = TrussModel("Three-member truss")
        n1 = Node((0.0, 0.0))
        n2 = Node((0.0, 120.0))
        n3 = Node((120.0, 120.0))
        n4 = Node((120.0, 0.0))

        e1 = Truss((n1, n2), E, A)
        e2 = Truss((n1, n3), E, A)
        e3 = Truss((n1, n4), E, A)

        model.add_nodes([n1, n2, n3, n4])
        model.add_elements([e1, e2, e3])

        model.add_force(n1, (0.0, -P))
        model.add_constraint(n2, ux=0.0, uy=0.0)
        model.add_constraint(n3, ux=0.0, uy=0.0)
        model.add_constraint(n4, ux=0.0, uy=0.0)
        model.solve()

        np.testing.assert_allclose(
            [n1.ux, n1.uy],
            [0.004142135623731, -0.015857864376269],
            rtol=1e-11,
            atol=1e-13,
        )

        np.testing.assert_allclose(
            [[n2.fx, n2.fy], [n3.fx, n3.fy], [n4.fx, n4.fy]],
            [
                [0.0, 7928.932188134525],
                [2071.067811865475, 2071.067811865475],
                [-2071.067811865475, 0.0],
            ],
            rtol=1e-11,
            atol=1e-8,
        )

        np.testing.assert_allclose(
            [e1.f, e2.f, e3.f],
            [7928.932188134525, 2928.9321881345245, -2071.067811865475],
            rtol=1e-11,
            atol=1e-8,
        )

        np.testing.assert_allclose(
            [e1.s, e2.s, e3.s],
            [3964.4660940672625, 1464.4660940672622, -1035.5339059327375],
            rtol=1e-11,
            atol=1e-8,
        )

        # Support reactions balance the applied load.
        support_reaction = np.array(
            [
                n2.fx + n3.fx + n4.fx,
                n2.fy + n3.fy + n4.fy,
            ]
        )
        np.testing.assert_allclose(support_reaction, [0.0, P], atol=1e-8)

    def test_kattan_problem_5_1(self):
        """Regression of examples/truss/truss_02.py."""
        E = 210e9
        A = 0.005

        model = TrussModel("Kattan 5.1")
        n1 = Node((0.0, 0.0))
        n2 = Node((5.0, 7.0))
        n3 = Node((5.0, 0.0))
        n4 = Node((10.0, 7.0))
        n5 = Node((10.0, 0.0))
        n6 = Node((15.0, 0.0))
        nodes = [n1, n2, n3, n4, n5, n6]

        elements = [
            Truss((n1, n2), E, A),
            Truss((n1, n3), E, A),
            Truss((n2, n3), E, A),
            Truss((n2, n4), E, A),
            Truss((n2, n5), E, A),
            Truss((n3, n5), E, A),
            Truss((n4, n5), E, A),
            Truss((n4, n6), E, A),
            Truss((n5, n6), E, A),
        ]

        model.add_nodes(nodes)
        model.add_elements(elements)
        model.add_constraint(n1, ux=0.0, uy=0.0)
        model.add_constraint(n6, ux=0.0, uy=0.0)
        model.add_force(n2, (20e3, 0.0))
        model.solve()

        expected_displacements = np.array(
            [
                [0.0, 0.0],
                [2.083428184225859e-4, -3.333837238599144e-5],
                [1.058201058201058e-5, -3.333837238599144e-5],
                [1.765967866765542e-4, 1.066263542454018e-5],
                [2.116402116402116e-5, -5.155958679768204e-5],
                [0.0, 0.0],
            ]
        )
        actual_displacements = np.array([[node.ux, node.uy] for node in nodes])
        np.testing.assert_allclose(
            actual_displacements,
            expected_displacements,
            rtol=1e-10,
            atol=1e-12,
        )

        expected_support_reactions = np.array(
            [
                [-8888.888888888889, -9333.333333333334],
                [-11111.111111111111, 9333.333333333334],
            ]
        )
        np.testing.assert_allclose(
            [[n1.fx, n1.fy], [n6.fx, n6.fy]],
            expected_support_reactions,
            rtol=1e-10,
            atol=1e-7,
        )

        expected_element_forces = np.array(
            [
                11469.7670227235,
                2222.222222222222,
                0.0,
                -6666.666666666662,
                -11469.767022723498,
                2222.222222222221,
                9333.333333333336,
                -11469.767022723501,
                -4444.444444444443,
            ]
        )
        np.testing.assert_allclose(
            [element.f for element in elements],
            expected_element_forces,
            rtol=1e-10,
            atol=1e-7,
        )

        # The vertical member n4 -> n5 is in tension with the current
        # positive-tension sign convention.
        assert elements[6].f > 0.0

        total_support_reaction = np.array(
            [
                n1.fx + n6.fx,
                n1.fy + n6.fy,
            ]
        )
        np.testing.assert_allclose(
            total_support_reaction,
            [-20e3, 0.0],
            atol=1e-7,
        )
