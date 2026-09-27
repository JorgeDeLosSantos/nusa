"""Regression tests for the transitional Model -> Analysis delegation."""

import numpy as np
import pytest

from nusa import (
    BarModel,
    BeamModel,
    LinearTriangleModel,
    Model,
    Node,
    Spring,
    SpringModel,
    StaticResult,
    TrussModel,
    solve,
)


@pytest.mark.parametrize(
    "model_type",
    [SpringModel, BarModel, TrussModel, BeamModel, LinearTriangleModel],
)
def test_public_models_share_base_assemble_and_solve(model_type):
    assert model_type.assemble is Model.assemble
    assert model_type.solve is Model.solve


def _spring_problem(load=10.0):
    model = SpringModel("delegation")
    n1 = Node((0.0, 0.0))
    n2 = Node((0.0, 0.0))
    element = Spring((n1, n2), 100.0)
    model.add_nodes([n1, n2])
    model.add_element(element)
    model.add_constraint(n1, ux=0.0)
    model.add_force(n2, (load,))
    return model, n1, n2, element


def test_model_solve_returns_static_result_and_keeps_legacy_surface():
    model, n1, n2, element = _spring_problem()

    result = model.solve()

    assert isinstance(result, StaticResult)
    assert model._last_result is result
    np.testing.assert_allclose(result.displacements, [0.0, 0.1])
    np.testing.assert_allclose(model.displacements, result.displacements)
    np.testing.assert_allclose(model.nodal_forces, result.nodal_forces)
    np.testing.assert_allclose(model.reactions, result.reactions)

    assert n1.ux == pytest.approx(0.0)
    assert n2.ux == pytest.approx(0.1)
    assert n1.fx == pytest.approx(-10.0)
    assert n2.fx == pytest.approx(10.0)
    assert model.element_result(element) == result.element_result(element)


def test_model_solve_matches_top_level_analysis_path():
    legacy_model, _, _, _ = _spring_problem()
    direct_model, _, _, _ = _spring_problem()

    delegated = legacy_model.solve()
    direct = solve(direct_model)

    np.testing.assert_allclose(delegated.displacements, direct.displacements)
    np.testing.assert_allclose(delegated.nodal_forces, direct.nodal_forces)
    np.testing.assert_allclose(delegated.reactions, direct.reactions)
    assert delegated.element_results == direct.element_results


def test_legacy_node_mutation_does_not_change_model_element_result_source():
    model, _, n2, element = _spring_problem()

    result = model.solve()
    frozen = result.element_result(element)

    # Historical element properties still read node state during the
    # transition, but the model-level normalized result now delegates to the
    # frozen StaticResult.
    n2.ux = 99.0

    assert model.element_result(element) == frozen
    assert model.element_result(element)["force_j"] == pytest.approx(10.0)


def test_input_change_invalidates_transitional_last_result():
    model, _, n2, _ = _spring_problem()
    old_result = model.solve()

    assert model._last_result is old_result

    model.add_force(n2, (20.0,))

    assert not hasattr(model, "_last_result")
    with pytest.raises(RuntimeError, match="after solve"):
        _ = model.element_results

    new_result = model.solve()
    assert new_result is model._last_result
    assert new_result is not old_result
    np.testing.assert_allclose(old_result.displacements, [0.0, 0.1])
    np.testing.assert_allclose(new_result.displacements, [0.0, 0.2])


def test_legacy_assemble_uses_analysis_assembly_and_returns_copy():
    model, _, _, _ = _spring_problem()

    model.assemble()

    expected = np.array([[100.0, -100.0], [-100.0, 100.0]])
    np.testing.assert_allclose(model.stiffness_matrix, expected)

    exposed = model.stiffness_matrix
    exposed[:] = 999.0
    np.testing.assert_allclose(model.stiffness_matrix, expected)
