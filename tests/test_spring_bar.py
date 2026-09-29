"""Numerical regression tests for spring and bar finite elements."""

import numpy as np

from nusa import Bar, BarModel, Node, Spring, SpringModel


class TestSpringElement:
    def test_element_stiffness_matrix(self):
        n1, n2 = Node((0.0, 0.0)), Node((0.0, 0.0))
        element = Spring((n1, n2), 300.0)
        np.testing.assert_allclose(
            element.get_element_stiffness(),
            [[300.0, -300.0], [-300.0, 300.0]],
        )

    def test_compute_results_requires_explicit_displacements(self):
        element = Spring((Node((0.0, 0.0)), Node((0.0, 0.0))), 300.0)
        assert element.compute_results([0.0, 2.5]) == {
            "force_i": -750.0,
            "force_j": 750.0,
        }


class TestSpringModel:
    def test_single_spring_reference_solution(self):
        model = SpringModel("Single spring")
        n1, n2 = Node((0.0, 0.0)), Node((0.0, 0.0))
        element = Spring((n1, n2), 300.0)
        model.add_nodes([n1, n2])
        model.add_element(element)
        model.add_constraint(n1, ux=0.0)
        model.add_force(n2, (750.0,))

        result = model.solve()

        np.testing.assert_allclose(result.displacements, [0.0, 2.5])
        np.testing.assert_allclose(result.nodal_forces, [-750.0, 750.0])
        assert result.element_result(element) == {
            "force_i": -750.0,
            "force_j": 750.0,
        }

    def test_logan_example_2_1(self):
        model = SpringModel("Logan 2.1")
        n1, n2, n3, n4 = [Node((0.0, 0.0)) for _ in range(4)]
        elements = [
            Spring((n1, n3), 1000.0),
            Spring((n3, n4), 2000.0),
            Spring((n4, n2), 3000.0),
        ]
        model.add_nodes([n1, n2, n3, n4])
        model.add_elements(elements)
        model.add_constraint(n1, ux=0.0)
        model.add_constraint(n2, ux=0.0)
        model.add_force(n4, (5000.0,))

        result = model.solve()

        np.testing.assert_allclose(
            [result.displacement(n3)["ux"], result.displacement(n4)["ux"]],
            [10.0 / 11.0, 15.0 / 11.0],
        )
        np.testing.assert_allclose(
            [result.nodal_force(n1)["fx"], result.nodal_force(n2)["fx"]],
            [-10000.0 / 11.0, -45000.0 / 11.0],
        )

    def test_nonzero_prescribed_displacement(self):
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

        result = model.solve()

        np.testing.assert_allclose(
            result.displacements,
            [0.0, 0.01, 0.02, 0.03, 0.02],
        )
        np.testing.assert_allclose(
            result.nodal_forces,
            [-2000.0, 0.0, 0.0, 4000.0, -2000.0],
            atol=1e-10,
        )


class TestBarElement:
    def test_length_stiffness_and_explicit_results(self):
        n1, n2 = Node((0.0, 0.0)), Node((2.0, 0.0))
        element = Bar((n1, n2), E=200e9, A=0.001)

        assert np.isclose(element.L, 2.0)
        values = element.compute_results([0.0, 0.01])
        np.testing.assert_allclose(
            [values["force_i"], values["force_j"]],
            [-1.0e6, 1.0e6],
        )
        assert values["axial_force"] == 1.0e6
        assert values["axial_stress"] == 1.0e9


class TestBarModel:
    def test_logan_example_3_1(self):
        model = BarModel("Logan 3.1")
        n1, n2, n3, n4 = (
            Node((0.0, 0.0)),
            Node((30.0, 0.0)),
            Node((60.0, 0.0)),
            Node((90.0, 0.0)),
        )
        elements = [
            Bar((n1, n2), E=30e6, A=1.0),
            Bar((n2, n3), E=30e6, A=1.0),
            Bar((n3, n4), E=15e6, A=2.0),
        ]
        model.add_nodes([n1, n2, n3, n4])
        model.add_elements(elements)
        model.add_constraint(n1, ux=0.0)
        model.add_constraint(n4, ux=0.0)
        model.add_force(n2, (3000.0,))

        result = model.solve()

        np.testing.assert_allclose(result.displacements, [0.0, 0.002, 0.001, 0.0])
        np.testing.assert_allclose(
            result.nodal_forces,
            [-2000.0, 3000.0, 0.0, -1000.0],
            atol=1e-10,
        )
        np.testing.assert_allclose(
            [result.element_result(e)["axial_force"] for e in elements],
            [2000.0, -1000.0, -1000.0],
        )

    def test_nonzero_prescribed_displacement_with_single_unknown(self):
        E, A = 210e6, 0.003
        model = BarModel("Kattan 3.1")
        n1, n2, n3 = Node((0.0, 0.0)), Node((1.5, 0.0)), Node((2.5, 0.0))
        model.add_nodes([n1, n2, n3])
        model.add_elements([Bar((n1, n2), E, A), Bar((n2, n3), E, A)])
        model.add_constraint(n1, ux=0.0)
        model.add_constraint(n3, ux=0.002)
        model.add_force(n2, (-10.0,))

        result = model.solve()
        k1, k2 = E * A / 1.5, E * A
        expected_u2 = (-10.0 + k2 * 0.002) / (k1 + k2)

        assert np.isclose(result.displacement(n2)["ux"], expected_u2)


def test_bar_chain_prescribed_displacement():
    model = BarModel("bar chain")
    nodes = [Node((float(i), 0.0)) for i in range(4)]
    model.add_nodes(nodes)
    model.add_elements([
        Bar((nodes[i], nodes[i + 1]), E=100.0, A=1.0)
        for i in range(3)
    ])
    model.add_constraint(nodes[0], ux=0.0)
    model.add_constraint(nodes[3], ux=0.03)

    result = model.solve()

    np.testing.assert_allclose(result.displacements, [0.0, 0.01, 0.02, 0.03])
