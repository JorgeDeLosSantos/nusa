"""Regression tests for model inputs and StaticResult force semantics."""

import numpy as np

from nusa import Beam, BeamModel, Node, Spring, SpringModel


def test_spring_force_semantics_live_on_result():
    model = SpringModel("force semantics")
    n1, n2 = Node((0.0, 0.0)), Node((0.0, 0.0))
    model.add_nodes([n1, n2])
    model.add_element(Spring((n1, n2), 100.0))
    model.add_constraint(n1, ux=0.0)
    model.add_force(n1, (10.0,))
    model.add_force(n2, (50.0,))

    np.testing.assert_allclose(model.applied_loads, [10.0, 50.0])
    assert model.applied_load(n1) == {"fx": 10.0}

    result = model.solve()

    np.testing.assert_allclose(result.nodal_forces, [-50.0, 50.0])
    np.testing.assert_allclose(result.reactions, [-60.0, 0.0])
    assert result.nodal_force(n1) == {"fx": -50.0}
    assert result.reaction(n1) == {"fx": -60.0}

    assert np.isclose(result.applied_loads.sum() + result.reactions.sum(), 0.0)


def test_beam_force_and_moment_components_are_result_owned():
    model = BeamModel("beam force semantics")
    n1, n2 = Node((0.0, 0.0)), Node((1.0, 0.0))
    model.add_nodes([n1, n2])
    model.add_element(Beam((n1, n2), E=1.0, I=1.0))
    model.add_constraint(n1, uy=0.0, ur=0.0)
    model.add_force(n2, (-1.0,))

    result = model.solve()

    np.testing.assert_allclose(result.applied_loads, [0.0, 0.0, -1.0, 0.0])
    np.testing.assert_allclose(result.nodal_forces, [1.0, 1.0, -1.0, 0.0])
    np.testing.assert_allclose(result.reactions, [1.0, 1.0, 0.0, 0.0])
    assert result.reaction(n1) == {"fy": 1.0, "m": 1.0}


def test_result_arrays_are_defensive_copies():
    model = SpringModel("copies")
    n1, n2 = Node((0.0, 0.0)), Node((0.0, 0.0))
    model.add_nodes([n1, n2])
    model.add_element(Spring((n1, n2), 10.0))
    model.add_constraint(n1, ux=0.0)
    model.add_force(n2, (5.0,))
    result = model.solve()

    for array in (
        result.applied_loads,
        result.displacements,
        result.nodal_forces,
        result.reactions,
    ):
        array[:] = 123.0

    np.testing.assert_allclose(result.applied_loads, [0.0, 5.0])
    np.testing.assert_allclose(result.displacements, [0.0, 0.5])
    np.testing.assert_allclose(result.nodal_forces, [-5.0, 5.0])
    np.testing.assert_allclose(result.reactions, [-5.0, 0.0])


def test_model_retains_only_problem_inputs_after_solve():
    model = BeamModel("inputs")
    n1, n2 = Node((0.0, 0.0)), Node((2.0, 0.0))
    model.add_nodes([n1, n2])
    model.add_element(Beam((n1, n2), E=1.0, I=1.0))
    model.add_constraint(n1, uy=0.0, ur=0.0)
    model.add_force(n2, (-3.0,))

    result = model.solve()

    np.testing.assert_allclose(model.applied_loads, [0.0, 0.0, -3.0, 0.0])
    np.testing.assert_allclose(model.prescribed_displacements[:2], [0.0, 0.0])
    assert np.isnan(model.prescribed_displacements[2:]).all()
    assert result.displacement(n2)["uy"] < 0.0
    assert not hasattr(model, "displacements")
