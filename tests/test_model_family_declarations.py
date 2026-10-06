"""Contract tests for declarative finite-element model families."""

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
    ("model_type", "element_type", "displacement_dofs", "force_dofs"),
    [
        (SpringModel, "spring", ("ux",), ("fx",)),
        (BarModel, "bar", ("ux",), ("fx",)),
        (TrussModel, "truss", ("ux", "uy"), ("fx", "fy")),
        (BeamModel, "beam", ("uy", "ur"), ("fy", "m")),
        (LinearTriangleModel, "triangle", ("ux", "uy"), ("fx", "fy")),
    ],
)
def test_model_family_contract_is_declared_on_class(
    model_type,
    element_type,
    displacement_dofs,
    force_dofs,
):
    model = model_type()

    assert model.element_type == element_type
    assert model.mtype == element_type
    assert model.displacement_dofs == displacement_dofs
    assert model.force_dofs == force_dofs
    assert model.dof == len(displacement_dofs)


@pytest.mark.parametrize(
    "model_type",
    [SpringModel, BarModel, TrussModel, BeamModel, LinearTriangleModel],
)
def test_model_families_share_generic_constraint_implementation(model_type):
    assert model_type.add_constraint is Model.add_constraint


@pytest.mark.parametrize(
    "model_type",
    [SpringModel, BarModel, TrussModel, LinearTriangleModel],
)
def test_standard_model_families_share_generic_force_implementation(model_type):
    assert model_type.add_force is Model.add_force


def test_beam_uses_generic_force_api_with_transverse_load_dofs_only():
    model = BeamModel()
    node = Node((0.0, 0.0))
    model.add_node(node)

    model.add_force(node, (-5.0,))
    model.add_moment(node, (2.0,))

    assert model.applied_load(node) == {"fy": -5.0, "m": 2.0}

    with pytest.raises(ValueError, match="exactly 1 component"):
        model.add_force(node, (-5.0, 2.0))


def test_base_model_remains_available_for_internal_generic_cases():
    model = Model("generic", "bar")

    assert model.mtype == "bar"
    assert model.dof == 0
