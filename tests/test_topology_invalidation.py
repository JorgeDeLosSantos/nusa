"""Regression tests for topology-driven analysis-state invalidation."""

import numpy as np

from nusa.core import Node
from nusa.element import Beam, Spring
from nusa.model import BeamModel, SpringModel


def test_topology_change_invalidates_and_rebuilds_spring_analysis_state():
    model = SpringModel("Topology invalidation")
    n1 = Node((0.0, 0.0))
    n2 = Node((0.0, 0.0))
    e1 = Spring((n1, n2), 100.0)

    model.add_nodes([n1, n2])
    model.add_element(e1)
    model.add_constraint(n1, ux=0.0)
    model.add_force(n2, (100.0,))
    model.solve()

    assert model._is_assembled is True
    assert np.isclose(n2.ux, 1.0)
    assert hasattr(model, "_K_reduced")
    assert hasattr(model, "_rhs_reduced")
    assert hasattr(model, "_free_dofs")
    assert hasattr(model, "_prescribed_dofs")

    n3 = Node((0.0, 0.0))
    model.add_node(n3)

    assert model._is_assembled is False
    assert not hasattr(model, "_K")
    assert not hasattr(model, "_K_reduced")
    assert not hasattr(model, "_rhs_reduced")
    assert not hasattr(model, "_free_dofs")
    assert not hasattr(model, "_prescribed_dofs")
    assert np.isclose(n1.ux, 0.0)
    assert np.isnan(n2.ux)
    assert np.isnan(n3.ux)

    e2 = Spring((n2, n3), 100.0)
    model.add_element(e2)

    assert model._is_assembled is False

    # solve() rebuilds automatically and restores only explicit inputs.
    model.solve()

    assert model._is_assembled is True
    assert model.stiffness_matrix.shape == (3, 3)
    np.testing.assert_allclose([n1.ux, n2.ux, n3.ux], [0.0, 1.0, 1.0])
    np.testing.assert_allclose([n1.fx, n2.fx, n3.fx], [-100.0, 100.0, 0.0])


def test_beam_rebuild_preserves_explicit_load_moment_and_constraints():
    model = BeamModel("Beam rebuild")
    n1 = Node((0.0, 0.0))
    n2 = Node((1.0, 0.0))
    e1 = Beam((n1, n2), E=1.0, I=1.0)

    model.add_nodes([n1, n2])
    model.add_element(e1)
    model.add_constraint(n1, ux=0.0, uy=0.0, ur=0.0)
    model.add_force(n2, (-1.0,))
    model.add_moment(n2, (0.5,))
    model.solve()

    old_tip_displacement = n2.uy
    assert not np.isnan(old_tip_displacement)

    n3 = Node((2.0, 0.0))
    model.add_node(n3)
    model.add_element(Beam((n2, n3), E=1.0, I=1.0))

    assert model._is_assembled is False
    assert np.isclose(n1.uy, 0.0)
    assert np.isclose(n1.ur, 0.0)
    assert np.isnan(n2.uy)
    assert np.isnan(n2.ur)

    model.solve()

    np.testing.assert_allclose(
        model._u[:4],
        [0.0, 0.0, n2.uy, n2.ur],
    )
    np.testing.assert_allclose(
        model._f[:4],
        [0.0, 0.0, -1.0, 0.5],
    )
    assert model.stiffness_matrix.shape == (6, 6)

    # The added segment is unloaded and free at n3, so it carries no
    # additional end forces and does not change the response at n2.
    assert np.isclose(n2.uy, old_tip_displacement)
