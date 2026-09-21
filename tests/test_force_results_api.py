"""Regression tests for applied loads, nodal forces, and reactions."""

import numpy as np
import pytest

from nusa.core import Node
from nusa.element import Beam, Spring
from nusa.model import BeamModel, SpringModel


def test_spring_force_api_separates_applied_load_nodal_force_and_reaction():
    model = SpringModel("Force semantics")
    n1 = Node((0.0, 0.0))
    n2 = Node((0.0, 0.0))
    model.add_nodes([n1, n2])
    model.add_element(Spring((n1, n2), 100.0))

    model.add_constraint(n1, ux=0.0)
    model.add_force(n1, (10.0,))
    model.add_force(n2, (50.0,))

    np.testing.assert_allclose(model.applied_loads, [10.0, 50.0])
    assert model.applied_load(n1) == {"fx": 10.0}
    assert model.applied_load(n2) == {"fx": 50.0}

    with pytest.raises(RuntimeError, match="after solve"):
        _ = model.nodal_forces
    with pytest.raises(RuntimeError, match="after solve"):
        _ = model.reactions

    model.solve()

    np.testing.assert_allclose(model.nodal_forces, [-50.0, 50.0])
    np.testing.assert_allclose(model.reactions, [-60.0, 0.0])

    assert model.nodal_force(n1) == {"fx": -50.0}
    assert model.nodal_force(n2) == {"fx": 50.0}
    assert model.reaction(n1) == {"fx": -60.0}
    assert model.reaction(n2) == {"fx": 0.0}

    # Node force attributes remain the solved generalized nodal forces (K @ u),
    # not the support reactions.
    assert np.isclose(n1.fx, -50.0)
    assert np.isclose(n2.fx, 50.0)

    # Applied loads plus reactions satisfy global equilibrium.
    assert np.isclose(model.applied_loads.sum() + model.reactions.sum(), 0.0)


def test_beam_force_api_uses_named_force_and_moment_components():
    model = BeamModel("Beam force semantics")
    n1 = Node((0.0, 0.0))
    n2 = Node((1.0, 0.0))
    model.add_nodes([n1, n2])
    model.add_element(Beam((n1, n2), E=1.0, I=1.0))

    model.add_constraint(n1, ux=0.0, uy=0.0, ur=0.0)
    model.add_force(n2, (-1.0,))
    model.solve()

    np.testing.assert_allclose(model.applied_loads, [0.0, 0.0, -1.0, 0.0])
    np.testing.assert_allclose(model.nodal_forces, [1.0, 1.0, -1.0, 0.0])
    np.testing.assert_allclose(model.reactions, [1.0, 1.0, 0.0, 0.0])

    assert model.applied_load(n1) == {"fy": 0.0, "m": 0.0}
    assert model.applied_load(n2) == {"fy": -1.0, "m": 0.0}
    assert model.reaction(n1) == {"fy": 1.0, "m": 1.0}
    assert model.reaction(n2) == {"fy": 0.0, "m": 0.0}


def test_force_result_arrays_are_returned_as_copies():
    model = SpringModel("Read-only by copy")
    n1 = Node((0.0, 0.0))
    n2 = Node((0.0, 0.0))
    model.add_nodes([n1, n2])
    model.add_element(Spring((n1, n2), 10.0))
    model.add_constraint(n1, ux=0.0)
    model.add_force(n2, (5.0,))
    model.solve()

    applied = model.applied_loads
    nodal = model.nodal_forces
    reactions = model.reactions

    applied[:] = 123.0
    nodal[:] = 123.0
    reactions[:] = 123.0

    np.testing.assert_allclose(model.applied_loads, [0.0, 5.0])
    np.testing.assert_allclose(model.nodal_forces, [-5.0, 5.0])
    np.testing.assert_allclose(model.reactions, [-5.0, 0.0])



def test_model_displacement_api_is_symmetric_with_force_results():
    model = BeamModel("Displacement semantics")
    n1 = Node((0.0, 0.0))
    n2 = Node((2.0, 0.0))
    model.add_nodes([n1, n2])
    model.add_element(Beam((n1, n2), E=1.0, I=1.0))

    model.add_constraint(n1, uy=0.0, ur=0.0)
    model.add_force(n2, (-3.0,))

    prescribed = model.prescribed_displacements
    assert prescribed.shape == (4,)
    np.testing.assert_allclose(prescribed[:2], [0.0, 0.0])
    assert np.isnan(prescribed[2])
    assert np.isnan(prescribed[3])
    assert model.prescribed_displacement(n1) == {"uy": 0.0, "ur": 0.0}
    assert np.isnan(model.prescribed_displacement(n2)["uy"])
    assert np.isnan(model.prescribed_displacement(n2)["ur"])

    with pytest.raises(RuntimeError, match="after solve"):
        _ = model.displacements
    with pytest.raises(RuntimeError, match="after solve"):
        model.displacement(n2)

    model.solve()

    np.testing.assert_allclose(
        model.displacements,
        [n1.uy, n1.ur, n2.uy, n2.ur],
    )
    assert model.displacement(n1) == {"uy": n1.uy, "ur": n1.ur}
    assert model.displacement(n2) == {"uy": n2.uy, "ur": n2.ur}

    copied = model.displacements
    copied[:] = 999.0
    np.testing.assert_allclose(
        model.displacements,
        [n1.uy, n1.ur, n2.uy, n2.ur],
    )


def test_result_vectors_are_invalidated_after_input_change():
    model = SpringModel("Result invalidation")
    n1 = Node((0.0, 0.0))
    n2 = Node((0.0, 0.0))
    model.add_nodes([n1, n2])
    model.add_element(Spring((n1, n2), 100.0))
    model.add_constraint(n1, ux=0.0)
    model.add_force(n2, (50.0,))
    model.solve()

    np.testing.assert_allclose(model.displacements, [0.0, 0.5])

    model.add_force(n2, (80.0,))

    with pytest.raises(RuntimeError, match="after solve"):
        _ = model.displacements
    with pytest.raises(RuntimeError, match="after solve"):
        _ = model.nodal_forces
    with pytest.raises(RuntimeError, match="after solve"):
        _ = model.reactions

    np.testing.assert_allclose(model.applied_loads, [0.0, 80.0])
    np.testing.assert_allclose(model.prescribed_displacements[:1], [0.0])
    assert np.isnan(model.prescribed_displacements[1])
