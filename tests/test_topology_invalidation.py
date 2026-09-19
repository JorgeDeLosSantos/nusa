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

    assert model.IS_KG_BUILDED is True
    assert np.isclose(n2.ux, 1.0)
    assert hasattr(model, "solved_u")

    n3 = Node((0.0, 0.0))
    model.add_node(n3)

    assert model.IS_KG_BUILDED is False
    assert not hasattr(model, "KG")
    assert not hasattr(model, "solved_u")
    assert np.isclose(n1.ux, 0.0)
    assert np.isnan(n2.ux)
    assert np.isnan(n3.ux)

    e2 = Spring((n2, n3), 100.0)
    model.add_element(e2)

    assert model.IS_KG_BUILDED is False

    # solve() rebuilds automatically and restores only explicit inputs.
    model.solve()

    assert model.IS_KG_BUILDED is True
    assert model.KG.shape == (3, 3)
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

    assert model.IS_KG_BUILDED is False
    assert np.isclose(n1.uy, 0.0)
    assert np.isclose(n1.ur, 0.0)
    assert np.isnan(n2.uy)
    assert np.isnan(n2.ur)

    model.solve()

    i1 = model._get_node_index(n1)
    i2 = model._get_node_index(n2)

    assert np.isclose(model.U[i1]["uy"], 0.0)
    assert np.isclose(model.U[i1]["ur"], 0.0)
    assert np.isclose(model.F[i2]["fy"], -1.0)
    assert np.isclose(model.F[i2]["m"], 0.5)
    assert model.KG.shape == (6, 6)
    assert not np.isclose(n2.uy, old_tip_displacement)
