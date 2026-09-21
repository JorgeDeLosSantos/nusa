"""API-freeze contract tests for the 0.3.0 beta line."""

import numpy as np
import pytest

from nusa import (
    BarModel,
    BeamModel,
    LinearTriangleModel,
    Model,
    Node,
    SpringModel,
    TrussModel,
)


@pytest.mark.parametrize(
    "coordinates",
    [
        (0.0,),
        (0.0, 1.0, 2.0),
        (np.nan, 0.0),
        (np.inf, 0.0),
    ],
)
def test_node_requires_two_finite_coordinates(coordinates):
    with pytest.raises(ValueError):
        Node(coordinates)


def test_model_rejects_non_node_objects():
    model = SpringModel()

    with pytest.raises(TypeError, match="Node instances"):
        model.add_node((0.0, 0.0))


def test_node_element_registration_is_not_public():
    node = Node((0.0, 0.0))

    assert not hasattr(node, "add_element")
    assert hasattr(node, "_add_element")


@pytest.mark.parametrize(
    ("model", "force"),
    [
        (SpringModel(), (1.0, 2.0)),
        (BarModel(), (1.0, 2.0)),
        (TrussModel(), (1.0,)),
        (LinearTriangleModel(), (1.0,)),
        (BeamModel(), (1.0, 2.0)),
    ],
)
def test_force_component_count_is_strict(model, force):
    node = Node((0.0, 0.0))
    model.add_node(node)

    with pytest.raises(ValueError, match="requires exactly"):
        model.add_force(node, force)


@pytest.mark.parametrize(
    "model",
    [
        SpringModel(),
        BarModel(),
        TrussModel(),
        LinearTriangleModel(),
        BeamModel(),
    ],
)
def test_force_components_must_be_finite(model):
    node = Node((0.0, 0.0))
    model.add_node(node)

    values = (np.nan,) if isinstance(model, BeamModel) else [0.0] * len(model.force_dofs)
    if not isinstance(model, BeamModel):
        values[-1] = np.nan

    with pytest.raises(ValueError, match="finite"):
        model.add_force(node, values)


@pytest.mark.parametrize(
    ("model", "constraint"),
    [
        (SpringModel(), {"uy": 0.0}),
        (BarModel(), {"uy": 0.0}),
        (TrussModel(), {"ur": 0.0}),
        (LinearTriangleModel(), {"ur": 0.0}),
        (BeamModel(), {"ux": 0.0}),
    ],
)
def test_constraints_reject_inactive_dofs(model, constraint):
    node = Node((0.0, 0.0))
    model.add_node(node)

    with pytest.raises(ValueError, match="Unsupported constraint"):
        model.add_constraint(node, **constraint)


def test_beam_moment_accepts_exactly_one_finite_component():
    model = BeamModel()
    node = Node((0.0, 0.0))
    model.add_node(node)

    model.add_moment(node, (5.0,))
    assert model.applied_load(node) == {"fy": 0.0, "m": 5.0}

    with pytest.raises(ValueError, match="requires exactly"):
        model.add_moment(node, (1.0, 2.0))

    with pytest.raises(ValueError, match="finite"):
        model.add_moment(node, (np.nan,))


def test_base_model_rejects_non_element_objects():
    model = Model("base", "spring")
    node = Node((0.0, 0.0))
    model.add_node(node)

    with pytest.raises(TypeError, match="Element instances"):
        model.add_element(object())
