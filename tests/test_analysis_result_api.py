"""Contract tests for the 0.4 analysis/result architecture."""

import numpy as np
import pytest

from nusa import (
    Bar,
    BarModel,
    Beam,
    BeamModel,
    LinearStaticAnalysis,
    LinearTriangle,
    LinearTriangleModel,
    Node,
    Spring,
    SpringModel,
    StaticResult,
    Truss,
    TrussModel,
    solve,
)


def _spring_problem(load=50.0):
    model = SpringModel("spring analysis")
    n1 = Node((0.0, 0.0))
    n2 = Node((0.0, 0.0))
    model.add_nodes([n1, n2])
    model.add_element(Spring((n1, n2), 100.0))
    model.add_constraint(n1, ux=0.0)
    model.add_force(n2, (load,))
    return model, n1, n2


def test_solve_returns_static_result_without_using_legacy_model_solve():
    model, n1, n2 = _spring_problem()

    result = solve(model)

    assert isinstance(result, StaticResult)
    np.testing.assert_allclose(result.displacements, [0.0, 0.5])
    np.testing.assert_allclose(result.nodal_forces, [-50.0, 50.0])
    np.testing.assert_allclose(result.reactions, [-50.0, 0.0])
    assert result.displacement(n2) == {"ux": pytest.approx(0.5)}
    assert result.reaction(n1) == {"fx": pytest.approx(-50.0)}

    # The new analysis path does not create legacy solved state on Model.
    with pytest.raises(RuntimeError, match="after solve"):
        _ = model.displacements


def test_explicit_analysis_matches_top_level_solve():
    model, _, _ = _spring_problem()

    direct = solve(model)
    explicit = LinearStaticAnalysis().solve(model)

    np.testing.assert_allclose(direct.displacements, explicit.displacements)
    np.testing.assert_allclose(direct.nodal_forces, explicit.nodal_forces)
    np.testing.assert_allclose(direct.reactions, explicit.reactions)


def test_new_analysis_does_not_write_solved_values_to_nodes():
    model, n1, n2 = _spring_problem()

    before = (
        n1.ux,
        n2.ux,
        n1.fx,
        n2.fx,
    )

    result = solve(model)

    after = (
        n1.ux,
        n2.ux,
        n1.fx,
        n2.fx,
    )

    assert np.isclose(before[0], after[0])
    assert np.isnan(before[1]) and np.isnan(after[1])
    assert np.isclose(before[2], after[2])
    assert np.isclose(before[3], after[3])
    assert result.displacement(n2)["ux"] == pytest.approx(0.5)


def test_result_is_stable_after_model_input_changes():
    model, _, n2 = _spring_problem(load=50.0)

    first = solve(model)
    model.add_force(n2, (80.0,))
    second = solve(model)

    np.testing.assert_allclose(first.applied_loads, [0.0, 50.0])
    np.testing.assert_allclose(first.displacements, [0.0, 0.5])
    np.testing.assert_allclose(first.reactions, [-50.0, 0.0])

    np.testing.assert_allclose(second.applied_loads, [0.0, 80.0])
    np.testing.assert_allclose(second.displacements, [0.0, 0.8])
    np.testing.assert_allclose(second.reactions, [-80.0, 0.0])


def test_result_arrays_are_defensive_copies():
    model, _, _ = _spring_problem()
    result = solve(model)

    arrays = (
        result.applied_loads,
        result.prescribed_displacements,
        result.displacements,
        result.nodal_forces,
        result.reactions,
        result.node_coordinates,
    )
    for values in arrays:
        values[...] = 999.0

    np.testing.assert_allclose(result.applied_loads, [0.0, 50.0])
    np.testing.assert_allclose(result.displacements, [0.0, 0.5])
    np.testing.assert_allclose(result.nodal_forces, [-50.0, 50.0])
    np.testing.assert_allclose(result.reactions, [-50.0, 0.0])
    np.testing.assert_allclose(result.node_coordinates, [[0.0, 0.0], [0.0, 0.0]])


def test_result_rejects_foreign_node_lookup():
    model, _, _ = _spring_problem()
    result = solve(model)

    with pytest.raises(ValueError, match="does not belong"):
        result.displacement(Node((0.0, 0.0)))


def test_result_freezes_geometry_connectivity_and_labels():
    model = TrussModel("snapshot")
    n1 = Node((0.0, 0.0))
    n2 = Node((2.0, 0.0))
    element = Truss((n1, n2), E=100.0, A=2.0)
    model.add_nodes([n1, n2])
    model.add_element(element)
    model.add_constraint(n1, ux=0.0, uy=0.0)
    model.add_constraint(n2, uy=0.0)
    model.add_force(n2, (10.0, 0.0))

    result = solve(model)

    assert result.model_name == "snapshot"
    assert result.node_labels == (0, 1)
    assert result.element_labels == (0,)
    assert result.element_types == ("truss",)
    assert result.connectivity == ((0, 1),)
    assert result.displacement_dofs == ("ux", "uy")
    assert result.force_dofs == ("fx", "fy")
    np.testing.assert_allclose(result.node_coordinates, [[0.0, 0.0], [2.0, 0.0]])


@pytest.mark.parametrize(
    "builder, expected",
    [
        (
            lambda: _bar_problem(),
            [0.0, 0.1],
        ),
        (
            lambda: _truss_problem(),
            [0.0, 0.0, 0.1, 0.0],
        ),
        (
            lambda: _beam_problem(),
            [0.0, 0.0, -0.26666666666666666, -0.2],
        ),
    ],
)
def test_new_analysis_preserves_reference_displacements(builder, expected):
    model = builder()

    result = solve(model)

    np.testing.assert_allclose(result.displacements, expected)


def _bar_problem():
    model = BarModel("bar")
    n1 = Node((0.0, 0.0))
    n2 = Node((2.0, 0.0))
    model.add_nodes([n1, n2])
    model.add_element(Bar((n1, n2), E=100.0, A=2.0))
    model.add_constraint(n1, ux=0.0)
    model.add_force(n2, (10.0,))
    return model


def _truss_problem():
    model = TrussModel("truss")
    n1 = Node((0.0, 0.0))
    n2 = Node((2.0, 0.0))
    model.add_nodes([n1, n2])
    model.add_element(Truss((n1, n2), E=100.0, A=2.0))
    model.add_constraint(n1, ux=0.0, uy=0.0)
    model.add_constraint(n2, uy=0.0)
    model.add_force(n2, (10.0, 0.0))
    return model


def _beam_problem():
    model = BeamModel("beam")
    n1 = Node((0.0, 0.0))
    n2 = Node((2.0, 0.0))
    model.add_nodes([n1, n2])
    model.add_element(Beam((n1, n2), E=100.0, I=1.0))
    model.add_constraint(n1, uy=0.0, ur=0.0)
    model.add_force(n2, (-10.0,))
    return model


def test_new_analysis_handles_nonzero_prescribed_displacement():
    model = SpringModel("nonzero prescribed")
    n1 = Node((0.0, 0.0))
    n2 = Node((0.0, 0.0))
    n3 = Node((0.0, 0.0))
    model.add_nodes([n1, n2, n3])
    model.add_elements([
        Spring((n1, n2), 100.0),
        Spring((n2, n3), 100.0),
    ])
    model.add_constraint(n1, ux=0.0)
    model.add_constraint(n3, ux=0.03)

    result = solve(model)

    np.testing.assert_allclose(result.displacements, [0.0, 0.015, 0.03])


def test_new_analysis_solves_linear_triangle_without_node_writeback():
    model = LinearTriangleModel("triangle")
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

    result = solve(model)

    np.testing.assert_allclose(
        result.displacements,
        [0.0, 0.0, 9.1e-8, 0.0, 0.0, 0.0],
        atol=1e-15,
    )
    assert np.isnan(n2.ux)
    assert np.isnan(n2.uy)


def test_singular_system_policy_matches_legacy_behavior():
    model = SpringModel("singular")
    n1 = Node((0.0, 0.0))
    n2 = Node((0.0, 0.0))
    model.add_nodes([n1, n2])
    model.add_element(Spring((n1, n2), 100.0))

    with pytest.raises(np.linalg.LinAlgError, match="Singular stiffness matrix"):
        solve(model)
