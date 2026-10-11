"""Regression tests for snapshot independence under topology changes."""

import numpy as np

from nusa import Beam, BeamModel, Node, Spring, SpringModel

from nusa import Material, Section

def _make_beam(nodes, E, I):
    return Beam(nodes, material=Material(E=E), section=Section(I=I))



def test_topology_change_does_not_mutate_old_result():
    model = SpringModel("topology")
    n1, n2 = Node((0.0, 0.0)), Node((0.0, 0.0))
    model.add_nodes([n1, n2])
    model.add_element(Spring((n1, n2), 100.0))
    model.add_constraint(n1, ux=0.0)
    model.add_force(n2, (100.0,))

    old = model.solve()
    np.testing.assert_allclose(old.displacements, [0.0, 1.0])

    n3 = Node((0.0, 0.0))
    model.add_node(n3)
    model.add_element(Spring((n2, n3), 100.0))
    model.add_constraint(n3, ux=0.0)

    new = model.solve()

    np.testing.assert_allclose(old.displacements, [0.0, 1.0])
    assert len(old.node_labels) == 2
    assert len(new.node_labels) == 3


def test_beam_topology_change_preserves_problem_inputs_and_old_snapshot():
    model = BeamModel("beam topology")
    n1, n2 = Node((0.0, 0.0)), Node((1.0, 0.0))
    model.add_nodes([n1, n2])
    model.add_element(_make_beam((n1, n2), E=1.0, I=1.0))
    model.add_constraint(n1, uy=0.0, ur=0.0)
    model.add_force(n2, (-1.0,))
    model.add_moment(n2, (0.5,))

    old = model.solve()

    n3 = Node((2.0, 0.0))
    model.add_node(n3)
    model.add_element(_make_beam((n2, n3), E=1.0, I=1.0))
    model.add_constraint(n3, uy=0.0)

    new = model.solve()

    assert model.applied_load(n2) == {"fy": -1.0, "m": 0.5}
    assert old.node_labels == (0, 1)
    assert new.node_labels == (0, 1, 2)
