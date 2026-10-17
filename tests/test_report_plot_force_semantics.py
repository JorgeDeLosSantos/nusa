"""Regression tests for reporting and problem-plot force semantics."""

import matplotlib.pyplot as plt

import nusa.visualization as viz
from nusa import (
    LinearTriangle,
    LinearTriangleModel,
    Node,
    Spring,
    SpringModel,
    Truss,
    TrussModel,
    plot_model,
)

from nusa import Material, Section
from nusa import Material

def _make_triangle(nodes, E, nu, t):
    return LinearTriangle(nodes, material=Material(E=E, nu=nu), thickness=t)



def _make_truss(nodes, E, A):
    return Truss(nodes, material=Material(E=E), section=Section(A=A))



def test_spring_result_report_separates_force_quantities():
    model = SpringModel("Report semantics")
    n1, n2 = Node((0, 0)), Node((0, 0))
    model.add_nodes([n1, n2])
    model.add_element(Spring((n1, n2), 100))
    model.add_constraint(n1, ux=0)
    model.add_force(n2, (50,))

    report = model.solve().simple_report(report_type="string")

    assert "APPLIED LOADS" in report
    assert "NODAL FORCES (K @ U)" in report
    assert "REACTIONS" in report
    assert "50" in report and "-50" in report


def test_truss_problem_plot_shows_applied_load_direction(monkeypatch):
    model = TrussModel("Plot semantics")
    n1, n2 = Node((0, 0)), Node((1, 0))
    model.add_nodes([n1, n2])
    model.add_element(_make_truss((n1, n2), 100, 1))
    model.add_constraint(n1, ux=0, uy=0)
    model.add_constraint(n2, uy=0)
    model.add_force(n2, (10, 0))

    arrows = []
    monkeypatch.setattr(
        viz,
        "_draw_force_arrow",
        lambda ax, x, y, axis, direction, size:
            arrows.append((x, y, axis, direction)),
    )

    plot_model(model)

    assert arrows == [(1.0, 0.0, "x", 1)]
    plt.close("all")


def test_triangle_problem_plot_preserves_negative_load_direction(monkeypatch):
    model = LinearTriangleModel("Negative load plot")
    n1, n2, n3 = Node((0, 0)), Node((1, 0)), Node((0, 1))
    model.add_nodes([n1, n2, n3])
    model.add_element(_make_triangle((n1, n2, n3), 1000, 0.25, 0.5))
    model.add_force(n2, (-10, -5))

    arrows = []
    monkeypatch.setattr(
        viz,
        "_draw_force_arrow",
        lambda ax, x, y, axis, direction, size:
            arrows.append((axis, x, y, direction)),
    )

    plot_model(model)

    assert ("x", 1.0, 0.0, -1) in arrows
    assert ("y", 1.0, 0.0, -1) in arrows
    plt.close("all")
