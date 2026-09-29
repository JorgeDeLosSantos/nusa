"""Contract tests for explicit element response and StaticResult ownership."""

import numpy as np
import pytest

from nusa import (
    Bar, BarModel, Beam, BeamModel, LinearTriangle, LinearTriangleModel,
    Node, Spring, SpringModel, Truss, TrussModel,
)


@pytest.mark.parametrize(
    ("model_builder", "expected_keys"),
    [
        ("spring", {"force_i", "force_j"}),
        ("bar", {"force_i", "force_j", "axial_force", "axial_stress"}),
        ("truss", {"axial_force", "axial_stress"}),
        ("beam", {"shear_force_i", "shear_force_j", "bending_moment_i", "bending_moment_j"}),
        ("triangle", {"stress_xx", "stress_yy", "stress_xy", "strain_xx", "strain_yy", "strain_xy"}),
    ],
)
def test_all_public_families_expose_canonical_result_keys(model_builder, expected_keys):
    if model_builder == "spring":
        model = SpringModel(); n1, n2 = Node((0,0)), Node((0,0)); e = Spring((n1,n2),100)
        model.add_nodes([n1,n2]); model.add_element(e); model.add_constraint(n1,ux=0); model.add_force(n2,(10,))
    elif model_builder == "bar":
        model = BarModel(); n1, n2 = Node((0,0)), Node((2,0)); e = Bar((n1,n2),100,2)
        model.add_nodes([n1,n2]); model.add_element(e); model.add_constraint(n1,ux=0); model.add_force(n2,(10,))
    elif model_builder == "truss":
        model = TrussModel(); n1, n2 = Node((0,0)), Node((2,0)); e = Truss((n1,n2),100,2)
        model.add_nodes([n1,n2]); model.add_element(e); model.add_constraint(n1,ux=0,uy=0); model.add_constraint(n2,uy=0); model.add_force(n2,(10,0))
    elif model_builder == "beam":
        model = BeamModel(); n1, n2 = Node((0,0)), Node((2,0)); e = Beam((n1,n2),100,1)
        model.add_nodes([n1,n2]); model.add_element(e); model.add_constraint(n1,uy=0,ur=0); model.add_force(n2,(-10,))
    else:
        model = LinearTriangleModel(); n1,n2,n3 = Node((0,0)),Node((1,.5)),Node((0,1)); e = LinearTriangle((n1,n2,n3),200e9,.3,.1)
        model.add_nodes([n1,n2,n3]); model.add_element(e); model.add_constraint(n1,ux=0,uy=0); model.add_constraint(n3,ux=0,uy=0); model.add_force(n2,(1000,0))

    result = model.solve()
    assert set(result.element_result(e)) == expected_keys


def test_bar_axial_sign_convention_is_positive_in_tension():
    e = Bar((Node((0,0)), Node((2,0))), E=100.0, A=2.0)
    values = e.compute_results([0.0, 0.1])
    assert values["axial_force"] > 0
    assert values["axial_stress"] > 0


def test_element_objects_do_not_expose_solution_dependent_properties():
    elements = [
        Spring((Node((0,0)), Node((0,0))), 1.0),
        Bar((Node((0,0)), Node((1,0))), 1.0, 1.0),
        Truss((Node((0,0)), Node((1,0))), 1.0, 1.0),
        Beam((Node((0,0)), Node((1,0))), 1.0, 1.0),
        LinearTriangle((Node((0,0)), Node((1,0)), Node((0,1))), 1.0, 0.25, 1.0),
    ]
    for element in elements:
        for name in ("f", "s", "fx", "fy", "m", "sx", "sy", "sxy", "ex", "ey", "exy"):
            assert not hasattr(element, name)
