"""Contract tests for result-based reporting."""

import numpy as np
import pytest

from nusa import Node, Spring, SpringModel, simple_report, solve


def _spring_result():
    model = SpringModel("result report")
    n1 = Node((0.0, 0.0))
    n2 = Node((1.0, 0.0))
    element = Spring((n1, n2), 100.0)
    model.add_nodes([n1, n2])
    model.add_element(element)
    model.add_constraint(n1, ux=0.0)
    model.add_force(n2, (25.0,))
    return model, n1, n2, element, solve(model)


def test_top_level_simple_report_consumes_static_result():
    model, _, _, _, result = _spring_result()

    report = simple_report(result, report_type="string")

    assert "NuSA Simple Report" in report
    assert "Model: result report" in report
    assert "NODAL DISPLACEMENTS" in report
    assert "APPLIED LOADS" in report
    assert "NODAL FORCES (K @ U)" in report
    assert "REACTIONS" in report
    assert "ELEMENT RESULTS" in report
    assert "FORCE I" in report
    assert "FORCE J" in report

    assert not hasattr(model, "displacements")


def test_static_result_report_matches_top_level_reporting_function():
    _, _, _, _, result = _spring_result()

    assert result.simple_report(report_type="string") == simple_report(
        result,
        report_type="string",
    )


def test_result_report_is_stable_after_model_and_node_mutation():
    model, n1, n2, element, result = _spring_result()
    baseline = result.simple_report(report_type="string")

    model.add_force(n2, (80.0,))
    n1.coordinates[:] = [99.0, 88.0]
    n2.coordinates[:] = [77.0, 66.0]
    element.label = "changed"

    assert result.simple_report(report_type="string") == baseline


def test_report_write_mode_uses_result_snapshot(tmp_path):
    _, _, _, _, result = _spring_result()
    path = tmp_path / "result-report.txt"

    returned = simple_report(result, report_type="write", fname=path)

    assert returned is None
    text = path.read_text(encoding="utf-8")
    assert "result report" in text
    assert "ELEMENT RESULTS" in text


def test_reporting_rejects_non_result_and_unknown_mode():
    with pytest.raises(TypeError, match="StaticResult"):
        simple_report(object(), report_type="string")

    _, _, _, _, result = _spring_result()
    with pytest.raises(ValueError, match="report_type"):
        simple_report(result, report_type="unknown")


def test_result_report_uses_frozen_numerical_values():
    _, _, _, _, result = _spring_result()

    report = result.simple_report(report_type="string")

    np.testing.assert_allclose(result.displacements, [0.0, 0.25])
    assert "-25" in report
    assert "25" in report
    assert "0.25" in report
