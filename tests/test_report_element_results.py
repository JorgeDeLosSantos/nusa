"""Regression tests for 0.4 report integration with normalized element results."""

import pytest

from nusa import (
    Bar,
    BarModel,
    Beam,
    BeamModel,
    LinearTriangle,
    LinearTriangleModel,
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


def test_model_report_delegates_to_latest_static_result():
    model = _truss_model()

    expected = model._last_result.simple_report(report_type="string")
    actual = model.simple_report(report_type="string")

    assert actual == expected


def test_report_uses_frozen_result_not_mutable_legacy_node_state():
    model = _truss_model()
    result = model._last_result
    element = model.elements[0]

    expected = result.element_result(element)
    model.nodes[-1].ux = 999.0
    model.nodes[-1].fx = 999.0

    report = result.simple_report(report_type="string")

    assert str(expected["axial_force"]) in report
    assert str(expected["axial_stress"]) in report
    assert "999" not in report
