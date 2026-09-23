"""Regression tests for 0.4 report integration with normalized element results."""

import pytest

from nusa import (
    Bar,
    BarModel,
    Beam,
    BeamModel,
    LinearTriangle,
    LinearTriangleModel,
    Model,
    Node,
    Spring,
    SpringModel,
    Truss,
    TrussModel,
)


def _spring_model():
    model = SpringModel("spring report")
    n1, n2 = Node((0.0, 0.0)), Node((0.0, 0.0))
    model.add_nodes([n1, n2])
    model.add_element(Spring((n1, n2), 100.0))
    model.add_constraint(n1, ux=0.0)
    model.add_force(n2, (10.0,))
    model.solve()
    return model


def _bar_model():
    model = BarModel("bar report")
    n1, n2 = Node((0.0, 0.0)), Node((2.0, 0.0))
    model.add_nodes([n1, n2])
    model.add_element(Bar((n1, n2), E=100.0, A=2.0))
    model.add_constraint(n1, ux=0.0)
    model.add_force(n2, (10.0,))
    model.solve()
    return model


def _truss_model():
    model = TrussModel("truss report")
    n1, n2 = Node((0.0, 0.0)), Node((2.0, 0.0))
    model.add_nodes([n1, n2])
    model.add_element(Truss((n1, n2), E=100.0, A=2.0))
    model.add_constraint(n1, ux=0.0, uy=0.0)
    model.add_constraint(n2, uy=0.0)
    model.add_force(n2, (10.0, 0.0))
    model.solve()
    return model


def _beam_model():
    model = BeamModel("beam report")
    n1, n2 = Node((0.0, 0.0)), Node((2.0, 0.0))
    model.add_nodes([n1, n2])
    model.add_element(Beam((n1, n2), E=100.0, I=1.0))
    model.add_constraint(n1, uy=0.0, ur=0.0)
    model.add_force(n2, (-10.0,))
    model.solve()
    return model


def _triangle_model():
    model = LinearTriangleModel("triangle report")
    n1 = Node((0.0, 0.0))
    n2 = Node((1.0, 0.5))
    n3 = Node((0.0, 1.0))
    model.add_nodes([n1, n2, n3])
    model.add_element(
        LinearTriangle((n1, n2, n3), E=200e9, nu=0.3, t=0.1)
    )
    model.add_constraint(n1, ux=0.0, uy=0.0)
    model.add_constraint(n3, ux=0.0, uy=0.0)
    model.add_force(n2, (1000.0, 0.0))
    model.solve()
    return model


@pytest.mark.parametrize(
    ("builder", "headers"),
    [
        (_spring_model, ("FORCE I", "FORCE J")),
        (
            _bar_model,
            ("FORCE I", "FORCE J", "AXIAL FORCE", "AXIAL STRESS"),
        ),
        (_truss_model, ("AXIAL FORCE", "AXIAL STRESS")),
        (
            _beam_model,
            (
                "SHEAR FORCE I",
                "SHEAR FORCE J",
                "BENDING MOMENT I",
                "BENDING MOMENT J",
            ),
        ),
        (
            _triangle_model,
            (
                "STRESS XX",
                "STRESS YY",
                "STRESS XY",
                "STRAIN XX",
                "STRAIN YY",
                "STRAIN XY",
            ),
        ),
    ],
)
def test_element_report_uses_canonical_result_names(builder, headers):
    model = builder()

    report = model.simple_report(report_type="string")

    for header in headers:
        assert header in report


def test_bar_report_exposes_physical_axial_force_and_stress():
    model = _bar_model()

    report = model.simple_report(report_type="string")

    assert "AXIAL FORCE" in report
    assert "AXIAL STRESS" in report
    assert "10" in report
    assert "5" in report


@pytest.mark.parametrize(
    "model_type",
    [SpringModel, BarModel, TrussModel, BeamModel, LinearTriangleModel],
)
def test_public_models_share_base_element_report_implementation(model_type):
    assert model_type._get_element_results is Model._get_element_results


def test_element_report_reads_normalized_model_result_api(monkeypatch):
    model = _truss_model()

    monkeypatch.setattr(
        model,
        "element_result",
        lambda element: {
            "axial_force": 1234.5,
            "axial_stress": 6789.0,
        },
    )

    report = model.simple_report(report_type="string")

    assert "1234.5" in report
    assert "6789" in report
