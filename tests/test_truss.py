"""Numerical regression tests for the 2D truss element and model."""

import numpy as np

from nusa import Node, Truss, TrussModel

from nusa import Material, Section

def _make_truss(nodes, E, A):
    return Truss(nodes, material=Material(E=E), section=Section(A=A))



class TestTrussElement:
    def test_geometry_and_stiffness(self):
        n1, n2 = Node((0.0, 0.0)), Node((3.0, 4.0))
        element = _make_truss((n1, n2), E=200.0, A=2.0)

        assert np.isclose(element.L, 5.0)
        assert np.isclose(element.theta, np.arctan2(4.0, 3.0))
        np.testing.assert_allclose(
            element.get_element_stiffness(),
            [
                [28.8, 38.4, -28.8, -38.4],
                [38.4, 51.2, -38.4, -51.2],
                [-28.8, -38.4, 28.8, 38.4],
                [-38.4, -51.2, 38.4, 51.2],
            ],
        )

    def test_explicit_axial_force_and_stress(self):
        element = _make_truss((Node((0.0, 0.0)), Node((3.0, 4.0))), E=200.0, A=2.0)

        values = element.compute_results([0.0, 0.0, 0.006, 0.008])
        assert np.isclose(values["axial_force"], 0.8)
        assert np.isclose(values["axial_stress"], 0.4)

        rigid = element.compute_results([0.25, -0.1, 0.25, -0.1])
        assert np.isclose(rigid["axial_force"], 0.0, atol=1e-12)


class TestTrussModel:
    def test_three_member_reference_problem(self):
        E, A, P = 30e6, 2.0, 10e3
        model = TrussModel("Three-member truss")
        n1, n2, n3, n4 = (
            Node((0.0, 0.0)),
            Node((0.0, 120.0)),
            Node((120.0, 120.0)),
            Node((120.0, 0.0)),
        )
        elements = [
            _make_truss((n1, n2), E, A),
            _make_truss((n1, n3), E, A),
            _make_truss((n1, n4), E, A),
        ]
        model.add_nodes([n1, n2, n3, n4])
        model.add_elements(elements)
        model.add_force(n1, (0.0, -P))
        for node in (n2, n3, n4):
            model.add_constraint(node, ux=0.0, uy=0.0)

        result = model.solve()

        np.testing.assert_allclose(
            [result.displacement(n1)["ux"], result.displacement(n1)["uy"]],
            [0.004142135623731, -0.015857864376269],
            rtol=1e-11,
            atol=1e-13,
        )
        np.testing.assert_allclose(
            [result.element_result(e)["axial_force"] for e in elements],
            [7928.932188134525, 2928.9321881345245, -2071.067811865475],
            rtol=1e-11,
            atol=1e-8,
        )

    def test_kattan_problem_5_1(self):
        E, A = 210e9, 0.005
        model = TrussModel("Kattan 5.1")
        nodes = [
            Node((0.0, 0.0)), Node((5.0, 7.0)), Node((5.0, 0.0)),
            Node((10.0, 7.0)), Node((10.0, 0.0)), Node((15.0, 0.0)),
        ]
        n1, n2, n3, n4, n5, n6 = nodes
        elements = [
            _make_truss((n1, n2), E, A), _make_truss((n1, n3), E, A),
            _make_truss((n2, n3), E, A), _make_truss((n2, n4), E, A),
            _make_truss((n2, n5), E, A), _make_truss((n3, n5), E, A),
            _make_truss((n4, n5), E, A), _make_truss((n4, n6), E, A),
            _make_truss((n5, n6), E, A),
        ]
        model.add_nodes(nodes)
        model.add_elements(elements)
        model.add_constraint(n1, ux=0.0, uy=0.0)
        model.add_constraint(n6, ux=0.0, uy=0.0)
        model.add_force(n2, (20e3, 0.0))

        result = model.solve()

        expected = np.array([
            [0.0, 0.0],
            [2.083428184225859e-4, -3.333837238599144e-5],
            [1.058201058201058e-5, -3.333837238599144e-5],
            [1.765967866765542e-4, 1.066263542454018e-5],
            [2.116402116402116e-5, -5.155958679768204e-5],
            [0.0, 0.0],
        ])
        np.testing.assert_allclose(
            result.displacements.reshape(-1, 2),
            expected,
            rtol=1e-10,
            atol=1e-12,
        )
        np.testing.assert_allclose(
            [result.element_result(e)["axial_force"] for e in elements],
            [
                11469.7670227235, 2222.222222222222, 0.0,
                -6666.666666666662, -11469.767022723498,
                2222.222222222221, 9333.333333333336,
                -11469.767022723501, -4444.444444444443,
            ],
            rtol=1e-10,
            atol=1e-7,
        )


def test_truss_nonzero_prescribed_displacement():
    model = TrussModel("prescribed truss")
    n1, n2, n3 = Node((0.0, 0.0)), Node((1.0, 0.0)), Node((2.0, 0.0))
    e1, e2 = _make_truss((n1, n2), 100.0, 1.0), _make_truss((n2, n3), 100.0, 1.0)
    model.add_nodes([n1, n2, n3])
    model.add_elements([e1, e2])
    model.add_constraint(n1, ux=0.0, uy=0.0)
    model.add_constraint(n2, uy=0.0)
    model.add_constraint(n3, ux=0.03, uy=0.0)

    result = model.solve()

    np.testing.assert_allclose(result.displacements[::2], [0.0, 0.015, 0.03])
    np.testing.assert_allclose(
        [result.element_result(e)["axial_force"] for e in (e1, e2)],
        [1.5, 1.5],
    )
