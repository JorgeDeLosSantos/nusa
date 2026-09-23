"""Contract tests for the normalized 0.4 element-results API."""

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
)


def _solved_spring():
    model = SpringModel("spring results")
    n1 = Node((0.0, 0.0))
    n2 = Node((0.0, 0.0))
    element = Spring((n1, n2), 100.0)
    model.add_nodes([n1, n2])
    model.add_element(element)
    model.add_constraint(n1, ux=0.0)
    model.add_force(n2, (10.0,))
    model.solve()
    return model, element


def _solved_bar():
    model = BarModel("bar results")
    n1 = Node((0.0, 0.0))
    n2 = Node((2.0, 0.0))
    element = Bar((n1, n2), E=100.0, A=2.0)
    model.add_nodes([n1, n2])
    model.add_element(element)
    model.add_constraint(n1, ux=0.0)
    model.add_force(n2, (10.0,))
    model.solve()
    return model, element


def _solved_truss():
    model = TrussModel("truss results")
    n1 = Node((0.0, 0.0))
    n2 = Node((2.0, 0.0))
    element = Truss((n1, n2), E=100.0, A=2.0)
    model.add_nodes([n1, n2])
    model.add_element(element)
    model.add_constraint(n1, ux=0.0, uy=0.0)
    model.add_constraint(n2, uy=0.0)
    model.add_force(n2, (10.0, 0.0))
    model.solve()
    return model, element


def _solved_beam():
    model = BeamModel("beam results")
    n1 = Node((0.0, 0.0))
    n2 = Node((2.0, 0.0))
    element = Beam((n1, n2), E=100.0, I=1.0)
    model.add_nodes([n1, n2])
    model.add_element(element)
    model.add_constraint(n1, uy=0.0, ur=0.0)
    model.add_force(n2, (-10.0,))
    model.solve()
    return model, element


def _solved_triangle():
    model = LinearTriangleModel("triangle results")
    n1 = Node((0.0, 0.0))
    n2 = Node((1.0, 0.5))
    n3 = Node((0.0, 1.0))
    element = LinearTriangle((n1, n2, n3), E=200e9, nu=0.3, t=0.1)
    model.add_nodes([n1, n2, n3])
    model.add_element(element)
    model.add_constraint(n1, ux=0.0, uy=0.0)
    model.add_constraint(n3, ux=0.0, uy=0.0)
    model.add_force(n2, (1000.0, 0.0))
    model.solve()
    return model, element


def test_element_result_requires_solved_model():
    model = SpringModel()
    n1 = Node((0.0, 0.0))
    n2 = Node((0.0, 0.0))
    element = Spring((n1, n2), 100.0)
    model.add_nodes([n1, n2])
    model.add_element(element)

    with pytest.raises(RuntimeError, match="after solve"):
        model.element_result(element)

    with pytest.raises(RuntimeError, match="after solve"):
        _ = model.element_results


def test_element_result_rejects_element_outside_model():
    model, _ = _solved_spring()
    foreign_nodes = (Node((0.0, 0.0)), Node((0.0, 0.0)))
    foreign = Spring(foreign_nodes, 100.0)

    with pytest.raises(ValueError, match="does not belong"):
        model.element_result(foreign)


@pytest.mark.parametrize(
    ("builder", "expected_keys"),
    [
        (_solved_spring, {"force_i", "force_j"}),
        (
            _solved_bar,
            {"force_i", "force_j", "axial_force", "axial_stress"},
        ),
        (_solved_truss, {"axial_force", "axial_stress"}),
        (
            _solved_beam,
            {
                "shear_force_i",
                "shear_force_j",
                "bending_moment_i",
                "bending_moment_j",
            },
        ),
        (
            _solved_triangle,
            {
                "stress_xx",
                "stress_yy",
                "stress_xy",
                "strain_xx",
                "strain_yy",
                "strain_xy",
            },
        ),
    ],
)
def test_all_public_models_expose_canonical_scalar_element_results(
    builder, expected_keys
):
    model, element = builder()

    result = model.element_result(element)

    assert set(result) == expected_keys
    assert all(isinstance(value, float) for value in result.values())
    assert all(np.isfinite(value) for value in result.values())


def test_spring_results_preserve_element_end_actions():
    model, element = _solved_spring()

    result = model.element_result(element)

    assert result == {"force_i": -10.0, "force_j": 10.0}


def test_bar_results_distinguish_end_actions_from_physical_axial_results():
    model, element = _solved_bar()

    result = model.element_result(element)

    assert result["force_i"] == pytest.approx(-10.0)
    assert result["force_j"] == pytest.approx(10.0)
    assert result["axial_force"] == pytest.approx(10.0)
    assert result["axial_stress"] == pytest.approx(5.0)

    # Compatibility surface remains unchanged: Bar.sx contains signed
    # end-action/A values rather than one physical constant stress.
    np.testing.assert_allclose(element.sx, [-5.0, 5.0])


def test_bar_axial_result_sign_is_positive_in_tension_and_negative_in_compression():
    tension_model, tension_element = _solved_bar()
    assert tension_model.element_result(tension_element)["axial_force"] > 0.0

    model = BarModel("bar compression")
    n1 = Node((0.0, 0.0))
    n2 = Node((2.0, 0.0))
    element = Bar((n1, n2), E=100.0, A=2.0)
    model.add_nodes([n1, n2])
    model.add_element(element)
    model.add_constraint(n1, ux=0.0)
    model.add_force(n2, (-10.0,))
    model.solve()

    result = model.element_result(element)
    assert result["axial_force"] == pytest.approx(-10.0)
    assert result["axial_stress"] == pytest.approx(-5.0)


def test_truss_results_preserve_current_axial_force_and_stress():
    model, element = _solved_truss()

    result = model.element_result(element)

    assert result["axial_force"] == pytest.approx(element.f)
    assert result["axial_stress"] == pytest.approx(element.s)
    assert result["axial_force"] == pytest.approx(10.0)
    assert result["axial_stress"] == pytest.approx(5.0)


def test_beam_results_preserve_end_action_signs():
    model, element = _solved_beam()

    result = model.element_result(element)

    assert result["shear_force_i"] == pytest.approx(10.0)
    assert result["shear_force_j"] == pytest.approx(-10.0)
    assert result["bending_moment_i"] == pytest.approx(20.0)
    assert result["bending_moment_j"] == pytest.approx(0.0, abs=1e-12)


def test_linear_triangle_results_preserve_reference_stress_and_strain():
    model, element = _solved_triangle()

    result = model.element_result(element)

    assert result["stress_xx"] == pytest.approx(20000.0)
    assert result["stress_yy"] == pytest.approx(6000.0)
    assert result["stress_xy"] == pytest.approx(0.0, abs=1e-10)
    assert result["strain_xx"] == pytest.approx(9.1e-8)
    assert result["strain_yy"] == pytest.approx(0.0, abs=1e-15)
    assert result["strain_xy"] == pytest.approx(0.0, abs=1e-15)


def test_element_results_preserve_model_insertion_order():
    model = SpringModel("ordered results")
    n1 = Node((0.0, 0.0))
    n2 = Node((0.0, 0.0))
    n3 = Node((0.0, 0.0))
    e1 = Spring((n1, n2), 100.0)
    e2 = Spring((n2, n3), 200.0)

    model.add_nodes([n1, n2, n3])
    model.add_elements([e1, e2])
    model.add_constraint(n1, ux=0.0)
    model.add_force(n3, (20.0,))
    model.solve()

    results = model.element_results

    assert isinstance(results, tuple)
    assert len(results) == 2
    assert results[0] == model.element_result(e1)
    assert results[1] == model.element_result(e2)


def test_returned_element_result_mapping_is_fresh():
    model, element = _solved_truss()

    result = model.element_result(element)
    result["axial_force"] = 999.0

    assert model.element_result(element)["axial_force"] == pytest.approx(10.0)


@pytest.mark.parametrize("mutation", ["load", "constraint"])
def test_input_change_invalidates_element_results(mutation):
    model, element = _solved_spring()
    node = model.nodes[-1]

    if mutation == "load":
        model.add_force(node, (20.0,))
    else:
        model.add_constraint(node, ux=0.1)

    with pytest.raises(RuntimeError, match="after solve"):
        model.element_result(element)

    with pytest.raises(RuntimeError, match="after solve"):
        _ = model.element_results


def test_topology_change_invalidates_element_results():
    model, element = _solved_spring()
    extra_node = Node((0.0, 0.0))

    model.add_node(extra_node)

    with pytest.raises(RuntimeError, match="after solve"):
        model.element_result(element)


def test_resolve_restores_element_results_after_input_change():
    model, element = _solved_spring()
    model.add_force(model.nodes[-1], (20.0,))

    with pytest.raises(RuntimeError, match="after solve"):
        model.element_result(element)

    model.solve()

    assert model.element_result(element) == {
        "force_i": pytest.approx(-20.0),
        "force_j": pytest.approx(20.0),
    }
