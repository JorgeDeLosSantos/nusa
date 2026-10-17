"""Numerical regression tests for Euler-Bernoulli beam elements and models."""

import numpy as np

from nusa import Beam, BeamModel, Node

from nusa import Material, Section

def _make_beam(nodes, E, I):
    return Beam(nodes, material=Material(E=E), section=Section(I=I))



class TestBeamElement:
    def test_length_stiffness_and_explicit_actions(self):
        element = _make_beam((Node((0.0, 0.0)), Node((2.0, 0.0))), E=1.0, I=1.0)
        assert np.isclose(element.L, 2.0)
        np.testing.assert_allclose(
            element.get_element_stiffness(),
            [
                [1.5, 1.5, -1.5, 1.5],
                [1.5, 2.0, -1.5, 1.0],
                [-1.5, -1.5, 1.5, -1.5],
                [1.5, 1.0, -1.5, 2.0],
            ],
        )
        values = element.compute_results([0.0, 0.0, -8.0, -6.0])
        np.testing.assert_allclose(
            [
                values["shear_force_i"], values["shear_force_j"],
                values["bending_moment_i"], values["bending_moment_j"],
            ],
            [3.0, -3.0, 6.0, 0.0],
        )


class TestBeamModel:
    def test_cantilever_tip_load(self):
        E, I, L, P = 29e6, 10.0, 10.0, 10e3
        model = BeamModel("cantilever")
        n1, n2 = Node((0.0, 0.0)), Node((L, 0.0))
        element = _make_beam((n1, n2), E, I)
        model.add_nodes([n1, n2])
        model.add_element(element)
        model.add_constraint(n1, uy=0.0, ur=0.0)
        model.add_force(n2, (-P,))

        result = model.solve()

        assert np.isclose(
            result.displacement(n2)["uy"],
            -P * L**3 / (3.0 * E * I),
        )
        assert np.isclose(
            result.displacement(n2)["ur"],
            -P * L**2 / (2.0 * E * I),
        )
        assert np.isclose(result.reaction(n1)["fy"], P)

    def test_logan_example_4_4(self):
        E, I, P, M, L = 210e9, 4e-4, 10e3, 20e3, 3.0
        model = BeamModel("Logan 4.4")
        n1, n2, n3 = Node((0.0, 0.0)), Node((L, 0.0)), Node((2*L, 0.0))
        e1, e2 = _make_beam((n1, n2), E, I), _make_beam((n2, n3), E, I)
        model.add_nodes([n1, n2, n3])
        model.add_elements([e1, e2])
        model.add_force(n2, (-P,))
        model.add_moment(n2, (M,))
        model.add_constraint(n1, uy=0.0, ur=0.0)
        model.add_constraint(n3, uy=0.0, ur=0.0)

        result = model.solve()

        np.testing.assert_allclose(
            [result.displacement(n2)["uy"], result.displacement(n2)["ur"]],
            [-1.339285714285714e-4, 8.928571428571429e-5],
            rtol=1e-11,
            atol=1e-14,
        )
        np.testing.assert_allclose(
            [result.reaction(n1)["fy"], result.reaction(n1)["m"],
             result.reaction(n3)["fy"], result.reaction(n3)["m"]],
            [1.0e4, 1.25e4, 0.0, -2.5e3],
            atol=1e-8,
        )
        assert np.isclose(result.element_result(e1)["shear_force_i"], 1.0e4)

    def test_simply_supported_eccentric_load(self):
        E, I, P = 29e6, 291.0, 35e3
        a, b, L = 60.0, 120.0, 180.0
        model = BeamModel("simply supported")
        n1, n2, n3 = Node((0.0, 0.0)), Node((a, 0.0)), Node((L, 0.0))
        model.add_nodes([n1, n2, n3])
        model.add_elements([_make_beam((n1, n2), E, I), _make_beam((n2, n3), E, I)])
        model.add_force(n2, (-P,))
        model.add_constraint(n1, uy=0.0)
        model.add_constraint(n3, uy=0.0)

        result = model.solve()

        expected = -(P * a**2 * b**2) / (3.0 * E * I * L)
        assert np.isclose(result.displacement(n2)["uy"], expected)
        assert np.isclose(result.reaction(n1)["fy"], P * b / L)
        assert np.isclose(result.reaction(n3)["fy"], P * a / L)

    def test_nonzero_support_settlement(self):
        model = BeamModel("settlement")
        n1, n2, n3 = Node((0.0, 0.0)), Node((1.0, 0.0)), Node((2.0, 0.0))
        model.add_nodes([n1, n2, n3])
        model.add_elements([_make_beam((n1, n2), 1.0, 1.0), _make_beam((n2, n3), 1.0, 1.0)])
        model.add_constraint(n1, uy=0.0, ur=0.0)
        model.add_constraint(n3, uy=0.1, ur=0.0)

        result = model.solve()

        np.testing.assert_allclose(
            [result.displacement(n2)["uy"], result.displacement(n2)["ur"]],
            [0.05, 0.075],
            atol=1e-12,
        )
