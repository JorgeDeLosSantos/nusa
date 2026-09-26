"""Contract tests for frozen element results in StaticResult."""

import numpy as np
import pytest

from nusa import (
    Bar,
    BarModel,
    Beam,
    BeamModel,
    LinearTriangle,
    LinearTriangleModel,
    Node,
    Spring,
    SpringModel,
    Truss,
    TrussModel,
    solve,
)


def _spring():
    model = SpringModel("spring")
    n1, n2 = Node((0.0, 0.0)), Node((0.0, 0.0))
    element = Spring((n1, n2), 100.0)
    model.add_nodes([n1, n2])
    model.add_element(element)
    model.add_constraint(n1, ux=0.0)
    model.add_force(n2, (10.0,))
    return model, element


def _bar():
    model = BarModel("bar")
    n1, n2 = Node((0.0, 0.0)), Node((2.0, 0.0))
    element = Bar((n1, n2), E=100.0, A=2.0)
    model.add_nodes([n1, n2])
    model.add_element(element)
    model.add_constraint(n1, ux=0.0)
    model.add_force(n2, (10.0,))
    return model, element


def _truss():
    model = TrussModel("truss")
    n1, n2 = Node((0.0, 0.0)), Node((2.0, 0.0))
    element = Truss((n1, n2), E=100.0, A=2.0)
    model.add_nodes([n1, n2])
    model.add_element(element)
    model.add_constraint(n1, ux=0.0, uy=0.0)
    model.add_constraint(n2, uy=0.0)
    model.add_force(n2, (10.0, 0.0))
    return model, element


def _beam():
    model = BeamModel("beam")
    n1, n2 = Node((0.0, 0.0)), Node((2.0, 0.0))
    element = Beam((n1, n2), E=100.0, I=1.0)
    model.add_nodes([n1, n2])
    model.add_element(element)
    model.add_constraint(n1, uy=0.0, ur=0.0)
    model.add_force(n2, (-10.0,))
    return model, element


def _triangle():
    model = LinearTriangleModel("triangle")
    n1 = Node((0.0, 0.0))
    n2 = Node((1.0, 0.5))
    n3 = Node((0.0, 1.0))
    element = LinearTriangle((n1, n2, n3), E=200e9, nu=0.3, t=0.1)
    model.add_nodes([n1, n2, n3])
    model.add_element(element)
    model.add_constraint(n1, ux=0.0, uy=0.0)
    model.add_constraint(n3, ux=0.0, uy=0.0)
    model.add_force(n2, (1000.0, 0.0))
    return model, element


@pytest.mark.parametrize(
    "builder, expected",
    [
        (_spring, {"force_i": -10.0, "force_j": 10.0}),
        (
            _bar,
            {
                "force_i": -10.0,
                "force_j": 10.0,
                "axial_force": 10.0,
                "axial_stress": 5.0,
            },
        ),
        (_truss, {"axial_force": 10.0, "axial_stress": 5.0}),
        (
            _beam,
            {
                "shear_force_i": 10.0,
                "shear_force_j": -10.0,
                "bending_moment_i": 20.0,
                "bending_moment_j": 0.0,
            },
        ),
    ],
)
def test_static_result_exposes_canonical_element_results(builder, expected):
    model, element = builder()

    result = solve(model)
    values = result.element_result(element)

    assert set(values) == set(expected)
    for name, expected_value in expected.items():
        assert values[name] == pytest.approx(expected_value, abs=1e-12)


def test_static_result_exposes_triangle_stress_and_strain():
    model, element = _triangle()

    values = solve(model).element_result(element)

    assert values["stress_xx"] == pytest.approx(20000.0)
    assert values["stress_yy"] == pytest.approx(6000.0)
    assert values["stress_xy"] == pytest.approx(0.0, abs=1e-10)
    assert values["strain_xx"] == pytest.approx(9.1e-8)
    assert values["strain_yy"] == pytest.approx(0.0, abs=1e-15)
    assert values["strain_xy"] == pytest.approx(0.0, abs=1e-15)


def test_element_results_are_frozen_and_return_fresh_mappings():
    model, element = _truss()
    result = solve(model)

    first = result.element_result(element)
    first["axial_force"] = 999.0

    assert result.element_result(element)["axial_force"] == pytest.approx(10.0)

    bulk = result.element_results
    bulk[0]["axial_force"] = 777.0

    assert result.element_results[0]["axial_force"] == pytest.approx(10.0)


def test_element_result_rejects_foreign_element():
    model, _ = _spring()
    result = solve(model)

    foreign_nodes = (Node((0.0, 0.0)), Node((0.0, 0.0)))
    foreign = Spring(foreign_nodes, 100.0)

    with pytest.raises(ValueError, match="does not belong"):
        result.element_result(foreign)


def test_element_results_preserve_element_order():
    model = SpringModel("ordered")
    n1, n2, n3 = (
        Node((0.0, 0.0)),
        Node((0.0, 0.0)),
        Node((0.0, 0.0)),
    )
    e1 = Spring((n1, n2), 100.0)
    e2 = Spring((n2, n3), 200.0)
    model.add_nodes([n1, n2, n3])
    model.add_elements([e1, e2])
    model.add_constraint(n1, ux=0.0)
    model.add_force(n3, (20.0,))

    result = solve(model)

    assert len(result.element_results) == 2
    assert result.element_results[0] == result.element_result(e1)
    assert result.element_results[1] == result.element_result(e2)


def test_old_element_results_survive_model_changes_and_resolve():
    model, element = _spring()
    old = solve(model)

    model.add_force(model.nodes[-1], (20.0,))
    new = solve(model)

    assert old.element_result(element) == {
        "force_i": pytest.approx(-10.0),
        "force_j": pytest.approx(10.0),
    }
    assert new.element_result(element) == {
        "force_i": pytest.approx(-20.0),
        "force_j": pytest.approx(20.0),
    }


@pytest.mark.parametrize("builder", [_spring, _bar, _truss, _beam, _triangle])
def test_new_element_result_path_does_not_write_solved_node_state(builder):
    model, _ = builder()
    before = [
        (node.ux, node.uy, node.ur, node.fx, node.fy, node.m)
        for node in model.nodes
    ]

    solve(model)

    after = [
        (node.ux, node.uy, node.ur, node.fx, node.fy, node.m)
        for node in model.nodes
    ]

    for old, new in zip(before, after):
        for old_value, new_value in zip(old, new):
            if np.isnan(old_value):
                assert np.isnan(new_value)
            else:
                assert new_value == pytest.approx(old_value)


@pytest.mark.parametrize(
    "builder, size",
    [(_spring, 2), (_bar, 2), (_truss, 4), (_beam, 4), (_triangle, 6)],
)
def test_element_compute_results_requires_explicit_local_displacements(builder, size):
    _, element = builder()

    with pytest.raises(ValueError):
        element.compute_results(np.zeros(size - 1))
