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
    assert model.get_applied_load(n1) == {"fx": 10.0}
    assert model.get_applied_load(n2) == {"fx": 50.0}

    with pytest.raises(RuntimeError, match="after solve"):
        _ = model.nodal_forces
    with pytest.raises(RuntimeError, match="after solve"):
        _ = model.reactions

    model.solve()

    np.testing.assert_allclose(model.nodal_forces, [-50.0, 50.0])
    np.testing.assert_allclose(model.reactions, [-60.0, 0.0])

    assert model.get_nodal_force(n1) == {"fx": -50.0}
    assert model.get_nodal_force(n2) == {"fx": 50.0}
    assert model.get_reaction(n1) == {"fx": -60.0}
    assert model.get_reaction(n2) == {"fx": 0.0}

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

    assert model.get_applied_load(n1) == {"fy": 0.0, "m": 0.0}
    assert model.get_applied_load(n2) == {"fy": -1.0, "m": 0.0}
    assert model.get_reaction(n1) == {"fy": 1.0, "m": 1.0}
    assert model.get_reaction(n2) == {"fy": 0.0, "m": 0.0}


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
