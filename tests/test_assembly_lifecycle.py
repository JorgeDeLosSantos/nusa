"""Regression tests for the stateless analysis lifecycle."""

import numpy as np

from nusa import Node, Spring, SpringModel, solve


def _spring_model():
    model = SpringModel("lifecycle")
    n1, n2 = Node((0.0, 0.0)), Node((0.0, 0.0))
    model.add_nodes([n1, n2])
    model.add_element(Spring((n1, n2), 100.0))
    model.add_constraint(n1, ux=0.0)
    model.add_force(n2, (50.0,))
    return model, n1, n2


def test_model_has_no_public_assembly_lifecycle():
    model, _, _ = _spring_model()

    assert not hasattr(model, "assemble")
    assert not hasattr(model, "stiffness_matrix")
    assert not hasattr(model, "_is_assembled")
    assert not hasattr(model, "_K")


def test_repeated_solve_assembles_fresh_without_model_cache():
    model, _, n2 = _spring_model()

    first = solve(model)
    model.add_force(n2, (80.0,))
    second = solve(model)

    np.testing.assert_allclose(first.displacements, [0.0, 0.5])
    np.testing.assert_allclose(second.displacements, [0.0, 0.8])
    assert not hasattr(model, "_K")


def test_constraint_change_affects_new_result_only():
    model, _, n2 = _spring_model()
    first = solve(model)

    model.add_constraint(n2, ux=0.1)
    second = solve(model)

    np.testing.assert_allclose(first.displacements, [0.0, 0.5])
    np.testing.assert_allclose(second.displacements, [0.0, 0.1])
