"""Numerical regression tests for the constant-strain linear triangle."""

import numpy as np

from nusa import LinearTriangle, LinearTriangleModel, Node

from nusa import Material

def _make_triangle(nodes, E, nu, t):
    return LinearTriangle(nodes, material=Material(E=E, nu=nu), thickness=t)



class TestLinearTriangleElement:
    def test_area_constitutive_B_and_stiffness(self):
        n1, n2, n3 = Node((0, 0)), Node((1, 0)), Node((0, 1))
        element = _make_triangle((n1, n2, n3), E=1000.0, nu=0.25, t=0.5)

        assert np.isclose(element.A, 0.5)
        np.testing.assert_allclose(
            element.B,
            [
                [-1.0, 0.0, 1.0, 0.0, 0.0, 0.0],
                [0.0, -1.0, 0.0, 0.0, 0.0, 1.0],
                [-1.0, -1.0, 0.0, 1.0, 1.0, 0.0],
            ],
        )
        np.testing.assert_allclose(
            element.D,
            [
                [1066.6666666666667, 266.6666666666667, 0.0],
                [266.6666666666667, 1066.6666666666667, 0.0],
                [0.0, 0.0, 400.0],
            ],
        )

    def test_affine_displacement_field_is_exact(self):
        n1, n2, n3 = Node((0, 0)), Node((1, 0)), Node((0, 1))
        element = _make_triangle((n1, n2, n3), E=200e9, nu=0.3, t=0.1)
        a, b, c = 1e-3, 2e-3, 0.1
        d, e, f = -0.5e-3, 3e-3, -0.2

        u_e = []
        for node in (n1, n2, n3):
            u_e.extend((
                a * node.x + b * node.y + c,
                d * node.x + e * node.y + f,
            ))

        expected_strain = np.array([a, e, b + d])
        expected_stress = element.D @ expected_strain

        np.testing.assert_allclose(element.compute_strain(u_e), expected_strain)
        np.testing.assert_allclose(element.compute_stress(u_e), expected_stress)

    def test_rigid_body_motion_produces_zero_strain_and_stress(self):
        n1, n2, n3 = Node((0, 0)), Node((1, 0)), Node((0, 1))
        element = _make_triangle((n1, n2, n3), E=200e9, nu=0.3, t=0.1)
        tx, ty, omega = 0.25, -0.4, 0.03
        u_e = []
        for node in (n1, n2, n3):
            u_e.extend((tx - omega * node.y, ty + omega * node.x))

        np.testing.assert_allclose(element.compute_strain(u_e), 0.0, atol=1e-14)
        np.testing.assert_allclose(element.compute_stress(u_e), 0.0, atol=1e-3)

    def test_clockwise_and_counterclockwise_are_equivalent(self):
        a, b, c = Node((0, 0)), Node((1, 0)), Node((0, 1))
        ccw = _make_triangle((a, b, c), E=1000.0, nu=0.25, t=0.5)
        cw = _make_triangle((a, c, b), E=1000.0, nu=0.25, t=0.5)

        dof_permutation = [0, 1, 4, 5, 2, 3]
        np.testing.assert_allclose(
            cw.get_element_stiffness()[np.ix_(dof_permutation, dof_permutation)],
            ccw.get_element_stiffness(),
            atol=1e-12,
        )

        u_ccw = [0.0, 0.0, 1e-3, -0.5e-3, 2e-3, 3e-3]
        u_cw = [0.0, 0.0, 2e-3, 3e-3, 1e-3, -0.5e-3]
        np.testing.assert_allclose(
            cw.compute_strain(u_cw),
            ccw.compute_strain(u_ccw),
            atol=1e-14,
        )

    def test_degenerate_triangle_rejected(self):
        with np.testing.assert_raises_regex(ValueError, "non-collinear"):
            _make_triangle(
                (Node((0, 0)), Node((1, 0)), Node((2, 0))),
                E=1000.0,
                nu=0.25,
                t=0.5,
            )


class TestLinearTriangleModel:
    def test_single_triangle_reference_problem(self):
        model = LinearTriangleModel("Single CST")
        n1, n2, n3 = Node((0, 0)), Node((1, 0.5)), Node((0, 1))
        element = _make_triangle((n1, n2, n3), E=200e9, nu=0.3, t=0.1)
        model.add_nodes([n1, n2, n3])
        model.add_element(element)
        model.add_constraint(n1, ux=0.0, uy=0.0)
        model.add_constraint(n3, ux=0.0, uy=0.0)
        model.add_force(n2, (1000.0, 0.0))

        result = model.solve()

        np.testing.assert_allclose(
            result.displacements,
            [0.0, 0.0, 9.1e-8, 0.0, 0.0, 0.0],
            atol=1e-14,
        )
        np.testing.assert_allclose(
            result.nodal_forces.reshape(-1, 2),
            [[-500.0, -300.0], [1000.0, 0.0], [-500.0, 300.0]],
            atol=1e-8,
        )
        values = result.element_result(element)
        np.testing.assert_allclose(
            [values["strain_xx"], values["strain_yy"], values["strain_xy"]],
            [9.1e-8, 0.0, 0.0],
            atol=1e-14,
        )
        np.testing.assert_allclose(
            [values["stress_xx"], values["stress_yy"], values["stress_xy"]],
            [20000.0, 6000.0, 0.0],
            atol=1e-7,
        )

    def test_three_element_plate_reference_problem(self):
        E, nu, t = 210e6, 0.3, 0.025
        model = LinearTriangleModel("Three-element plate")
        nodes = [
            Node((0.0, 0.0)), Node((0.5, 0.0)), Node((0.5, 0.25)),
            Node((0.0, 0.25)), Node((0.0, 0.5)),
        ]
        n1, n2, n3, n4, n5 = nodes
        elements = [
            _make_triangle((n1, n3, n4), E, nu, t),
            _make_triangle((n1, n2, n3), E, nu, t),
            _make_triangle((n4, n3, n5), E, nu, t),
        ]
        model.add_nodes(nodes)
        model.add_elements(elements)
        for node in (n1, n4, n5):
            model.add_constraint(node, ux=0.0, uy=0.0)
        model.add_force(n2, (9375.0, 0.0))
        model.add_force(n3, (9375.0, 0.0))

        result = model.solve()

        np.testing.assert_allclose(
            [
                [result.displacement(n2)["ux"], result.displacement(n2)["uy"]],
                [result.displacement(n3)["ux"], result.displacement(n3)["uy"]],
            ],
            [
                [0.005817197743558, 0.001967144381766],
                [0.003902081109925, 0.000931544442750],
            ],
            rtol=1e-11,
            atol=1e-14,
        )
        expected_stresses = np.array([
            [1800960.512273213, 540288.153681964, 150480.256136606],
            [2398078.975453576, -150480.256136607, -300960.512273212],
            [1800960.512273213, 540288.153681964, 150480.256136606],
        ])
        actual = np.array([
            [
                result.element_result(e)["stress_xx"],
                result.element_result(e)["stress_yy"],
                result.element_result(e)["stress_xy"],
            ]
            for e in elements
        ])
        np.testing.assert_allclose(actual, expected_stresses, rtol=1e-10, atol=1e-6)

    def test_nonzero_prescribed_displacement(self):
        model = LinearTriangleModel("Prescribed CST")
        n1, n2, n3 = Node((0, 0)), Node((1, 0)), Node((0, 1))
        model.add_nodes([n1, n2, n3])
        model.add_element(_make_triangle((n1, n2, n3), E=1000.0, nu=0.25, t=0.5))
        model.add_constraint(n1, ux=0.0, uy=0.0)
        model.add_constraint(n2, ux=0.01, uy=0.0)

        result = model.solve()

        np.testing.assert_allclose(
            [result.displacement(n3)["ux"], result.displacement(n3)["uy"]],
            [0.0, -0.0025],
            atol=1e-12,
        )
