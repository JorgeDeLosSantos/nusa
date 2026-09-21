"""Numerical regression tests for Euler-Bernoulli beam elements and models."""

import numpy as np
from nusa.core import Node
from nusa.element import Beam
from nusa.model import BeamModel


class TestBeamElement:
    def test_length_and_stiffness_matrix(self):
        n1 = Node((0.0, 0.0))
        n2 = Node((2.0, 0.0))
        element = Beam((n1, n2), E=1.0, I=1.0)

        assert np.isclose(element.L, 2.0)

        expected = np.array(
            [
                [1.5, 1.5, -1.5, 1.5],
                [1.5, 2.0, -1.5, 1.0],
                [-1.5, -1.5, 1.5, -1.5],
                [1.5, 1.0, -1.5, 2.0],
            ]
        )
        np.testing.assert_allclose(element.get_element_stiffness(), expected)

    def test_element_end_forces_from_nodal_displacements(self):
        n1 = Node((0.0, 0.0))
        n2 = Node((2.0, 0.0))
        n1.uy = 0.0
        n1.ur = 0.0
        n2.uy = -8.0
        n2.ur = -6.0

        element = Beam((n1, n2), E=1.0, I=1.0)

        np.testing.assert_allclose(element.fy.ravel(), [3.0, -3.0])
        np.testing.assert_allclose(element.m.ravel(), [6.0, 0.0])


class TestBeamModel:
    def test_cantilever_tip_load_matches_closed_form_solution(self):
        """One-element cantilever with a tip load has the exact cubic solution."""
        E = 29e6
        I = 10.0
        L = 10.0
        P = 10e3

        model = BeamModel("Cantilever")
        n1 = Node((0.0, 0.0))
        n2 = Node((L, 0.0))
        element = Beam((n1, n2), E, I)

        model.add_nodes([n1, n2])
        model.add_element(element)
        model.add_constraint(n1, ux=0.0, uy=0.0, ur=0.0)
        model.add_force(n2, (-P,))
        model.solve()

        expected_tip_displacement = -P * L**3 / (3.0 * E * I)
        expected_tip_rotation = -P * L**2 / (2.0 * E * I)

        assert np.isclose(n2.uy, expected_tip_displacement)
        assert np.isclose(n2.ur, expected_tip_rotation)
        assert np.isclose(n1.fy, P)
        assert np.isclose(n1.m, P * L)
        assert np.isclose(n2.fy, -P)
        assert np.isclose(n2.m, 0.0, atol=1e-10)

        np.testing.assert_allclose(element.fy.ravel(), [P, -P])
        np.testing.assert_allclose(element.m.ravel(), [P * L, 0.0], atol=1e-10)

    def test_logan_example_4_4(self):
        """Regression of examples/beam/beam_2.py (Logan, Example 4.4)."""
        E = 210e9
        I = 4e-4
        P = 10e3
        M = 20e3
        L = 3.0

        model = BeamModel("Logan 4.4")
        n1 = Node((0.0, 0.0))
        n2 = Node((L, 0.0))
        n3 = Node((2.0 * L, 0.0))
        e1 = Beam((n1, n2), E, I)
        e2 = Beam((n2, n3), E, I)

        model.add_nodes([n1, n2, n3])
        model.add_elements([e1, e2])
        model.add_force(n2, (-P,))
        model.add_moment(n2, (M,))
        model.add_constraint(n1, ux=0.0, uy=0.0, ur=0.0)
        model.add_constraint(n3, ux=0.0, uy=0.0, ur=0.0)
        model.solve()

        np.testing.assert_allclose(
            [n2.uy, n2.ur],
            [-1.339285714285714e-4, 8.928571428571429e-5],
            rtol=1e-11,
            atol=1e-14,
        )

        np.testing.assert_allclose(
            [[n1.fy, n1.m], [n3.fy, n3.m]],
            [[1.0e4, 1.25e4], [0.0, -2.5e3]],
            rtol=1e-11,
            atol=1e-8,
        )

        np.testing.assert_allclose(e1.fy.ravel(), [1.0e4, -1.0e4], atol=1e-8)
        np.testing.assert_allclose(e1.m.ravel(), [1.25e4, 1.75e4], atol=1e-8)
        np.testing.assert_allclose(e2.fy.ravel(), [0.0, 0.0], atol=1e-8)
        np.testing.assert_allclose(e2.m.ravel(), [2.5e3, -2.5e3], atol=1e-8)

    def test_simply_supported_beam_with_eccentric_point_load(self):
        """Regression of examples/beam/beam_4.py against closed-form beam theory."""
        E = 29e6
        I = 291.0
        P = 35e3
        a = 5.0 * 12.0
        b = 10.0 * 12.0
        L = a + b

        model = BeamModel("Simply supported beam")
        n1 = Node((0.0, 0.0))
        n2 = Node((a, 0.0))
        n3 = Node((L, 0.0))
        e1 = Beam((n1, n2), E, I)
        e2 = Beam((n2, n3), E, I)

        model.add_nodes([n1, n2, n3])
        model.add_elements([e1, e2])
        model.add_force(n2, (-P,))
        model.add_constraint(n1, ux=0.0, uy=0.0)
        model.add_constraint(n3, uy=0.0)
        model.solve()

        expected_displacement = -(P * a**2 * b**2) / (3.0 * E * I * L)
        expected_left_reaction = P * b / L
        expected_right_reaction = P * a / L

        assert np.isclose(n2.uy, expected_displacement)
        assert np.isclose(n1.fy, expected_left_reaction)
        assert np.isclose(n3.fy, expected_right_reaction)
        assert np.isclose(n1.m, 0.0, atol=1e-7)
        assert np.isclose(n3.m, 0.0, atol=1e-7)

        assert np.isclose(n1.fy + n3.fy, P)

    def test_nonzero_support_settlement(self):
        """Known solver limitation: support settlements are not assembled correctly."""
        E = 1.0
        I = 1.0

        model = BeamModel("Support settlement")
        n1 = Node((0.0, 0.0))
        n2 = Node((1.0, 0.0))
        n3 = Node((2.0, 0.0))
        e1 = Beam((n1, n2), E, I)
        e2 = Beam((n2, n3), E, I)

        model.add_nodes([n1, n2, n3])
        model.add_elements([e1, e2])
        model.add_constraint(n1, ux=0.0, uy=0.0, ur=0.0)
        model.add_constraint(n3, ux=0.0, uy=0.1, ur=0.0)
        model.solve()

        # With no external load and equal EI spans, the exact middle-node
        # response follows from the partitioned system K_uu u_u = -K_uk u_k.
        np.testing.assert_allclose(
            [n2.uy, n2.ur],
            [0.05, 0.075],
            atol=1e-12,
        )
