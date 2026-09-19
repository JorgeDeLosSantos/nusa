"""Regression tests for report and plot force semantics."""

import matplotlib.pyplot as plt
import numpy as np

from nusa.core import Node
from nusa.element import LinearTriangle, Spring, Truss
from nusa.model import LinearTriangleModel, SpringModel, TrussModel


def test_spring_report_separates_applied_loads_nodal_forces_and_reactions():
    model = SpringModel("Report semantics")
    n1 = Node((0.0, 0.0))
    n2 = Node((0.0, 0.0))
    model.add_nodes([n1, n2])
    model.add_element(Spring((n1, n2), 100.0))
    model.add_constraint(n1, ux=0.0)
    model.add_force(n2, (50.0,))
    model.solve()

    report = model.simple_report(report_type="string")

    assert "APPLIED LOADS" in report
    assert "NODAL FORCES (K @ U)" in report
    assert "REACTIONS" in report
    assert "50" in report
    assert "-50" in report


def test_truss_plot_defaults_to_applied_loads_and_reactions_are_optional(monkeypatch):
    model = TrussModel("Plot semantics")
    n1 = Node((0.0, 0.0))
    n2 = Node((1.0, 0.0))
    model.add_nodes([n1, n2])
    model.add_element(Truss((n1, n2), E=100.0, A=1.0))
    model.add_constraint(n1, ux=0.0, uy=0.0)
    model.add_constraint(n2, uy=0.0)
    model.add_force(n2, (10.0, 0.0))
    model.solve()

    x_arrows = []
    monkeypatch.setattr(
        model,
        "_draw_xforce",
        lambda axes, x, y, ddir=1, reaction=False:
            x_arrows.append((x, y, ddir, reaction)),
    )
    monkeypatch.setattr(model, "_draw_yforce", lambda *args, **kwargs: None)
    monkeypatch.setattr(model, "_draw_xconstraint", lambda *args, **kwargs: None)
    monkeypatch.setattr(model, "_draw_yconstraint", lambda *args, **kwargs: None)

    model.plot_model()
    assert x_arrows == [(1.0, 0.0, 1, False)]
    plt.close("all")

    x_arrows.clear()
    model.plot_model(show_reactions=True)
    assert (1.0, 0.0, 1, False) in x_arrows
    assert (0.0, 0.0, -1, True) in x_arrows
    assert len(x_arrows) == 2
    plt.close("all")


def test_linear_triangle_plot_preserves_negative_applied_load_direction(monkeypatch):
    model = LinearTriangleModel("Negative load plot")
    n1 = Node((0.0, 0.0))
    n2 = Node((1.0, 0.0))
    n3 = Node((0.0, 1.0))
    model.add_nodes([n1, n2, n3])
    model.add_element(LinearTriangle((n1, n2, n3), E=1000.0, nu=0.25, t=0.5))
    model.add_force(n2, (-10.0, -5.0))

    arrows = []
    monkeypatch.setattr(
        model,
        "_draw_xforce",
        lambda axes, x, y, ddir=1, reaction=False:
            arrows.append(("x", x, y, ddir, reaction)),
    )
    monkeypatch.setattr(
        model,
        "_draw_yforce",
        lambda axes, x, y, ddir=1, reaction=False:
            arrows.append(("y", x, y, ddir, reaction)),
    )
    monkeypatch.setattr(model, "_draw_xyconstraint", lambda *args, **kwargs: None)

    model.plot_model()

    assert ("x", 1.0, 0.0, -1, False) in arrows
    assert ("y", 1.0, 0.0, -1, False) in arrows
    plt.close("all")
