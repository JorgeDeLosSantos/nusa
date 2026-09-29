"""Contract tests for stateless Model -> Analysis delegation."""

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
def test_public_models_share_base_solve(model_type):
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


def test_model_solve_returns_static_result_without_storing_it():
    model, _, _, element = _spring_problem()

    result = model.solve()

    assert isinstance(result, StaticResult)
    np.testing.assert_allclose(result.displacements, [0.0, 0.1])
    assert result.element_result(element)["force_j"] == pytest.approx(10.0)

    for name in (
        "_last_result", "_u", "_f", "_nodal_forces", "_reactions",
        "_K", "_K_reduced", "_rhs_reduced", "_free_dofs", "_prescribed_dofs",
        "displacements", "nodal_forces", "reactions", "element_results",
        "stiffness_matrix", "assemble", "simple_report",
    ):
        assert not hasattr(model, name)


def test_model_solve_matches_top_level_solve():
    model_a, _, _, _ = _spring_problem()
    model_b, _, _, _ = _spring_problem()

    delegated = model_a.solve()
    direct = solve(model_b)

    np.testing.assert_allclose(delegated.displacements, direct.displacements)
    np.testing.assert_allclose(delegated.nodal_forces, direct.nodal_forces)
    np.testing.assert_allclose(delegated.reactions, direct.reactions)
    assert delegated.element_results == direct.element_results


def test_each_model_solve_returns_independent_snapshot():
    model, _, n2, _ = _spring_problem()
    first = model.solve()

    model.add_force(n2, (20.0,))
    second = model.solve()

    assert first is not second
    np.testing.assert_allclose(first.displacements, [0.0, 0.1])
    np.testing.assert_allclose(second.displacements, [0.0, 0.2])
