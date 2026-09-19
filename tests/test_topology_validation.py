"""Regression tests for model topology validation before assembly."""

import pytest

from nusa.core import Node
from nusa.element import Bar, Beam, LinearTriangle, Spring, Truss
from nusa.model import (
    BarModel,
    BeamModel,
    LinearTriangleModel,
    SpringModel,
    TrussModel,
)


def _spring_case():
    model = SpringModel("orphan spring")
    n1 = Node((0.0, 0.0))
    n2 = Node((1.0, 0.0))
    orphan = Node((2.0, 0.0))
    model.add_nodes([n1, n2, orphan])
    model.add_element(Spring((n1, n2), 100.0))
    return model, orphan


def _bar_case():
    model = BarModel("orphan bar")
    n1 = Node((0.0, 0.0))
    n2 = Node((1.0, 0.0))
    orphan = Node((2.0, 0.0))
    model.add_nodes([n1, n2, orphan])
    model.add_element(Bar((n1, n2), E=100.0, A=1.0))
    return model, orphan


def _truss_case():
    model = TrussModel("orphan truss")
    n1 = Node((0.0, 0.0))
    n2 = Node((1.0, 0.0))
    orphan = Node((2.0, 0.0))
    model.add_nodes([n1, n2, orphan])
    model.add_element(Truss((n1, n2), E=100.0, A=1.0))
    return model, orphan


def _beam_case():
    model = BeamModel("orphan beam")
    n1 = Node((0.0, 0.0))
    n2 = Node((1.0, 0.0))
    orphan = Node((2.0, 0.0))
    model.add_nodes([n1, n2, orphan])
    model.add_element(Beam((n1, n2), E=1.0, I=1.0))
    return model, orphan


def _triangle_case():
    model = LinearTriangleModel("orphan triangle")
    n1 = Node((0.0, 0.0))
    n2 = Node((1.0, 0.0))
    n3 = Node((0.0, 1.0))
    orphan = Node((2.0, 2.0))
    model.add_nodes([n1, n2, n3, orphan])
    model.add_element(LinearTriangle((n1, n2, n3), E=1000.0, nu=0.25, t=0.5))
    return model, orphan


@pytest.mark.parametrize(
    "builder",
    [_spring_case, _bar_case, _truss_case, _beam_case, _triangle_case],
)
def test_assemble_rejects_orphan_nodes_for_all_public_models(builder):
    model, orphan = builder()

    with pytest.raises(ValueError, match=str(orphan.label)):
        model.assemble()

    assert model._is_assembled is False
    assert not hasattr(model, "_K")


@pytest.mark.parametrize(
    "model",
    [
        SpringModel("empty spring"),
        BarModel("empty bar"),
        TrussModel("empty truss"),
        BeamModel("empty beam"),
        LinearTriangleModel("empty triangle"),
    ],
)
def test_assemble_rejects_models_without_elements(model):
    with pytest.raises(ValueError, match="without elements"):
        model.assemble()

    assert model._is_assembled is False
    assert not hasattr(model, "_K")


def test_solve_reports_topology_error_before_solver_singularity():
    model, orphan = _spring_case()

    with pytest.raises(ValueError, match=str(orphan.label)):
        model.solve()

    assert model._is_assembled is False
    assert not hasattr(model, "_K")


def test_linear_triangle_does_not_auto_restrain_orphan_node():
    model, orphan = _triangle_case()

    assert pytest.approx(float("nan"), nan_ok=True) == orphan.ux
    assert pytest.approx(float("nan"), nan_ok=True) == orphan.uy

    with pytest.raises(ValueError, match=str(orphan.label)):
        model.solve()

    assert pytest.approx(float("nan"), nan_ok=True) == orphan.ux
    assert pytest.approx(float("nan"), nan_ok=True) == orphan.uy
    assert model._prescribed_displacements.get(orphan) is None
