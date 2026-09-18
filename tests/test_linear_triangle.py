"""Numerical regression tests for the constant-strain linear triangle element."""

import numpy as np
import pytest

from nusa.core import Node
from nusa.element import LinearTriangle
from nusa.model import LinearTriangleModel


class TestLinearTriangleElement:
    def test_area_constitutive_strain_displacement_and_stiffness(self):
        n1 = Node((0.0, 0.0))
        n2 = Node((1.0, 0.0))
        n3 = Node((0.0, 1.0))
        element = LinearTriangle((n1, n2, n3), E=1000.0, nu=0.25, t=0.5)

        assert np.isclose(element.A, 0.5)

        expected_B = np.array(
            [
                [-1.0, 0.0, 1.0, 0.0, 0.0, 0.0],
                [0.0, -1.0, 0.0, 0.0, 0.0, 1.0],
                [-1.0, -1.0, 0.0, 1.0, 1.0, 0.0],
            ]
        )
        expected_D = np.array(
            [
                [1066.6666666666667, 266.6666666666667, 0.0],
                [266.6666666666667, 1066.6666666666667, 0.0],
                [0.0, 0.0, 400.0],
            ]
        )
        expected_K = np.array(
            [
                [366.6666666666667, 166.66666666666669, -266.6666666666667, -100.0, -100.0, -66.66666666666667],
                [166.66666666666669, 366.6666666666667, -66.66666666666667, -100.0, -100.0, -266.6666666666667],
                [-266.6666666666667, -66.66666666666667, 266.6666666666667, 0.0, 0.0, 66.66666666666667],
                [-100.0, -100.0, 0.0, 100.0, 100.0, 0.0],
                [-100.0, -100.0, 0.0, 100.0, 100.0, 0.0],
                [-66.66666666666667, -266.6666666666667, 66.66666666666667, 0.0, 0.0, 266.6666666666667],
            ]
        )

        np.testing.assert_allclose(element.B, expected_B)
        np.testing.assert_allclose(element.D, expected_D)
        np.testing.assert_allclose(element.get_element_stiffness(), expected_K)

    def test_affine_displacement_field_gives_exact_constant_strain_and_stress(self):
        n1 = Node((0.0, 0.0))
        n2 = Node((1.0, 0.0))
        n3 = Node((0.0, 1.0))
        element = LinearTriangle((n1, n2, n3), E=200e9, nu=0.3, t=0.1)

        # u = a*x + b*y + c, v = d*x + e*y + f
        a, b, c = 1.0e-3, 2.0e-3, 0.1
        d, e, f = -0.5e-3, 3.0e-3, -0.2

        for node in (n1, n2, n3):
            node.ux = a * node.x + b * node.y + c
            node.uy = d * node.x + e * node.y + f

        expected_strain = np.array([a, e, b + d])
        expected_stress = element.D @ expected_strain

        np.testing.assert_allclose(element.get_element_strains(), expected_strain)
        np.testing.assert_allclose(element.get_element_stresses(), expected_stress)
        assert np.isclose(element.ex, expected_strain[0])
        assert np.isclose(element.ey, expected_strain[1])
        assert np.isclose(element.exy, expected_strain[2])

    def test_rigid_body_motion_produces_zero_strain_and_stress(self):
        n1 = Node((0.0, 0.0))
        n2 = Node((1.0, 0.0))
        n3 = Node((0.0, 1.0))
        element = LinearTriangle((n1, n2, n3), E=200e9, nu=0.3, t=0.1)

        tx, ty, omega = 0.25, -0.4, 0.03
        for node in (n1, n2, n3):
            node.ux = tx - omega * node.y
            node.uy = ty + omega * node.x

        np.testing.assert_allclose(element.get_element_strains(), 0.0, atol=1e-14)
        np.testing.assert_allclose(element.get_element_stresses(), 0.0, atol=1e-3)

    @pytest.mark.xfail(
        strict=True,
        reason=(
            "Clockwise node ordering gives a negative signed area and therefore "
            "a negative-semidefinite stiffness matrix."
        ),
    )
    def test_clockwise_connectivity_does_not_produce_negative_stiffness(self):
        n1 = Node((0.0, 0.0))
        n2 = Node((0.0, 1.0))
        n3 = Node((1.0, 0.0))
        element = LinearTriangle((n1, n2, n3), E=1000.0, nu=0.25, t=0.5)

        stiffness = element.get_element_stiffness()
        eigenvalues = np.linalg.eigvalsh(stiffness)

        assert eigenvalues.min() >= -1e-10


class TestLinearTriangleModel:
    def test_single_triangle_reference_problem(self):
        """Regression of examples/linear_triangle/simple_triangle/simple_triangle.py."""
        model = LinearTriangleModel("Single CST")
        n1 = Node((0.0, 0.0))
        n2 = Node((1.0, 0.5))
        n3 = Node((0.0, 1.0))
        element = LinearTriangle((n1, n2, n3), E=200e9, nu=0.3, t=0.1)

        model.add_nodes([n1, n2, n3])
        model.add_element(element)
        model.add_constraint(n1, ux=0.0, uy=0.0)
        model.add_constraint(n3, ux=0.0, uy=0.0)
        model.add_force(n2, (1000.0, 0.0))
        model.solve()

        np.testing.assert_allclose([n2.ux, n2.uy], [9.1e-8, 0.0], atol=1e-14)
        np.testing.assert_allclose(
            [[n1.fx, n1.fy], [n2.fx, n2.fy], [n3.fx, n3.fy]],
            [[-500.0, -300.0], [1000.0, 0.0], [-500.0, 300.0]],
            atol=1e-8,
        )

        np.testing.assert_allclose(
            element.get_element_strains(),
            [9.1e-8, 0.0, 0.0],
            atol=1e-14,
        )
        np.testing.assert_allclose(
            element.get_element_stresses(),
            [20000.0, 6000.0, 0.0],
            atol=1e-7,
        )

        np.testing.assert_allclose(
            [n1.fx + n2.fx + n3.fx, n1.fy + n2.fy + n3.fy],
            [0.0, 0.0],
            atol=1e-8,
        )

    def test_three_element_plate_reference_problem(self):
        """Regression of examples/linear_triangle/simple_plate/plate_1.py."""
        E = 210e6
        nu = 0.3
        t = 0.025

        model = LinearTriangleModel("Three-element plate")
        n1 = Node((0.0, 0.0))
        n2 = Node((0.5, 0.0))
        n3 = Node((0.5, 0.25))
        n4 = Node((0.0, 0.25))
        n5 = Node((0.0, 0.5))
        nodes = [n1, n2, n3, n4, n5]

        e1 = LinearTriangle((n1, n3, n4), E, nu, t)
        e2 = LinearTriangle((n1, n2, n3), E, nu, t)
        e3 = LinearTriangle((n4, n3, n5), E, nu, t)
        elements = [e1, e2, e3]

        model.add_nodes(nodes)
        model.add_elements(elements)
        model.add_constraint(n1, ux=0.0, uy=0.0)
        model.add_constraint(n4, ux=0.0, uy=0.0)
        model.add_constraint(n5, ux=0.0, uy=0.0)
        model.add_force(n2, (9375.0, 0.0))
        model.add_force(n3, (9375.0, 0.0))
        model.solve()

        np.testing.assert_allclose(
            [[n2.ux, n2.uy], [n3.ux, n3.uy]],
            [
                [0.005817197743558, 0.001967144381766],
                [0.003902081109925, 0.000931544442750],
            ],
            rtol=1e-11,
            atol=1e-14,
        )

        np.testing.assert_allclose(
            [[n1.fx, n1.fy], [n4.fx, n4.fy], [n5.fx, n5.fy]],
            [
                [-8434.498399146214, -2436.299359658486],
                [-11256.00320170758, -940.501600853790],
                [940.501600853790, 3376.800960512274],
            ],
            rtol=1e-10,
            atol=1e-7,
        )

        expected_stresses = np.array(
            [
                [1800960.512273213, 540288.153681964, 150480.256136606],
                [2398078.975453576, -150480.256136607, -300960.512273212],
                [1800960.512273213, 540288.153681964, 150480.256136606],
            ]
        )
        np.testing.assert_allclose(
            [element.get_element_stresses() for element in elements],
            expected_stresses,
            rtol=1e-10,
            atol=1e-6,
        )

        # After solve(), nodal forces contain K @ u. Applied nodal loads and
        # support reactions therefore balance over the complete model.
        total_nodal_force = np.array(
            [
                sum(node.fx for node in nodes),
                sum(node.fy for node in nodes),
            ]
        )
        np.testing.assert_allclose(total_nodal_force, [0.0, 0.0], atol=1e-7)
