"""Regression tests for the model assembly/solution lifecycle."""

import numpy as np
import pytest

from nusa.core import Node
from nusa.element import Spring
from nusa.model import SpringModel


def _spring_model():
    model = SpringModel("Assembly lifecycle")
    n1 = Node((0.0, 0.0))
    n2 = Node((0.0, 0.0))
    model.add_nodes([n1, n2])
    model.add_element(Spring((n1, n2), 100.0))
    model.add_constraint(n1, ux=0.0)
    model.add_force(n2, (50.0,))
    return model, n1, n2


def test_stiffness_matrix_requires_assembly_and_is_returned_as_copy():
    model, _, _ = _spring_model()

    with pytest.raises(RuntimeError, match="assemble"):
        _ = model.stiffness_matrix

    model.assemble()
    np.testing.assert_allclose(
        model.stiffness_matrix,
        [[100.0, -100.0], [-100.0, 100.0]],
    )

    matrix = model.stiffness_matrix
    matrix[0, 0] = 999.0
    assert np.isclose(model.stiffness_matrix[0, 0], 100.0)


def test_solve_assembles_automatically():
    model, _, n2 = _spring_model()

    assert model._is_assembled is False
    model.solve()

    assert model._is_assembled is True
    assert np.isclose(n2.ux, 0.5)
    np.testing.assert_allclose(
        model.stiffness_matrix,
        [[100.0, -100.0], [-100.0, 100.0]],
    )


def test_load_change_invalidates_solution_but_preserves_assembly():
    model, n1, n2 = _spring_model()
    model.solve()

    assembled_matrix = model._K
    assert hasattr(model, "_nodal_forces")
    assert np.isclose(n2.ux, 0.5)

    model.add_force(n2, (80.0,))

    assert model._is_assembled is True
    assert model._K is assembled_matrix
    assert not hasattr(model, "_nodal_forces")
    assert not hasattr(model, "_reactions")
    assert np.isclose(n1.ux, 0.0)
    assert np.isnan(n2.ux)
    np.testing.assert_allclose(model._f, [0.0, 80.0])

    model.solve()
    assert model._K is assembled_matrix
    assert np.isclose(n2.ux, 0.8)


def test_topology_change_invalidates_assembly_and_legacy_names_are_absent():
    model, _, _ = _spring_model()
    model.assemble()

    model.add_node(Node((0.0, 0.0)))

    assert model._is_assembled is False
    assert not hasattr(model, "_K")
    with pytest.raises(RuntimeError, match="assemble"):
        _ = model.stiffness_matrix

    assert not hasattr(model, "KG")
    assert not hasattr(model, "IS_KG_BUILDED")
    assert not hasattr(model, "build_global_matrix")
