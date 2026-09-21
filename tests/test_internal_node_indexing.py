"""Regression tests for public node labels and model-owned solver indices."""

import numpy as np

from nusa.core import Node
from nusa.element import Bar, Beam, LinearTriangle, Spring, Truss
from nusa.model import BarModel, BeamModel, LinearTriangleModel, SpringModel, TrussModel


def test_spring_solver_accepts_string_labels_and_label_mutation():
    model = SpringModel("Labeled spring")
    n1 = Node((0.0, 0.0))
    n2 = Node((0.0, 0.0))
    n1.label = "fixed"
    n2.label = "tip"
    element = Spring((n1, n2), 100.0)

    model.add_nodes([n1, n2])
    model.add_element(element)

    # Public labels can change without changing model-owned solver indices.
    n1.label = "support-A"
    n2.label = "load-point"

    model.add_constraint(n1, ux=0.0)
    model.add_force(n2, (50.0,))
    model.solve()

    assert np.isclose(n2.ux, 0.5)
    assert np.isclose(n1.fx, -50.0)
    assert np.isclose(n2.fx, 50.0)
    np.testing.assert_allclose(model._u, [0.0, 0.5])
    np.testing.assert_allclose(model._f, [0.0, 50.0])


def test_bar_solver_accepts_sparse_integer_labels():
    model = BarModel("Sparse labels")
    n1 = Node((0.0, 0.0))
    n2 = Node((1.0, 0.0))
    n3 = Node((2.0, 0.0))
    n1.label, n2.label, n3.label = 10, 30, 80

    e1 = Bar((n1, n2), E=1.0, A=1.0)
    e2 = Bar((n2, n3), E=1.0, A=1.0)

    model.add_nodes([n1, n2, n3])
    model.add_elements([e1, e2])
    model.add_constraint(n1, ux=0.0)
    model.add_force(n3, (1.0,))
    model.solve()

    np.testing.assert_allclose([n1.ux, n2.ux, n3.ux], [0.0, 1.0, 2.0])
    np.testing.assert_allclose([n1.fx, n2.fx, n3.fx], [-1.0, 0.0, 1.0])
    assert model.stiffness_matrix.shape == (3, 3)


def test_truss_solver_accepts_string_labels():
    model = TrussModel("Labeled truss")
    n1 = Node((0.0, 0.0))
    n2 = Node((2.0, 0.0))
    n1.label, n2.label = "A", "B"
    element = Truss((n1, n2), E=100.0, A=2.0)

    model.add_nodes([n1, n2])
    model.add_element(element)
    model.add_constraint(n1, uy=0.0)
    model.add_constraint(n2, uy=0.0)
    model.add_force(n2, (10.0, 0.0))
    model.solve()

    assert np.isclose(n2.ux, 0.1)
    assert np.isclose(element.f, 10.0)

    report = model.simple_report(report_type="string")
    assert "A" in report
    assert "B" in report


def test_beam_solver_accepts_string_labels():
    model = BeamModel("Labeled beam")
    n1 = Node((0.0, 0.0))
    n2 = Node((1.0, 0.0))
    n1.label, n2.label = "fixed", "tip"
    element = Beam((n1, n2), E=1.0, I=1.0)

    model.add_nodes([n1, n2])
    model.add_element(element)
    model.add_constraint(n1, uy=0.0, ur=0.0)
    model.add_force(n2, (-1.0,))
    model.solve()

    assert np.isclose(n2.uy, -1.0 / 3.0)
    assert np.isclose(n2.ur, -0.5)


def test_linear_triangle_solver_accepts_string_labels():
    model = LinearTriangleModel("Labeled CST")
    n1 = Node((0.0, 0.0))
    n2 = Node((1.0, 0.5))
    n3 = Node((0.0, 1.0))
    n1.label, n2.label, n3.label = "left-bottom", "loaded", "left-top"
    element = LinearTriangle((n1, n2, n3), E=200e9, nu=0.3, t=0.1)

    model.add_nodes([n1, n2, n3])
    model.add_element(element)
    model.add_constraint(n1, uy=0.0)
    model.add_constraint(n3, uy=0.0)
    model.add_force(n2, (1000.0, 0.0))
    model.solve()

    np.testing.assert_allclose([n2.ux, n2.uy], [9.1e-8, 0.0], atol=1e-14)
    np.testing.assert_allclose(
        [n1.fx + n2.fx + n3.fx, n1.fy + n2.fy + n3.fy],
        [0.0, 0.0],
        atol=1e-8,
    )

    triangulation = model._get_tri()
    np.testing.assert_array_equal(triangulation.triangles, [[0, 1, 2]])
