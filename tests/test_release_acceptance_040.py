"""Release-acceptance contract for the NuSA 0.4.0 architecture."""

from pathlib import Path
import re

import numpy as np
import pytest

from nusa import (
    Bar,
    BarModel,
    Beam,
    BeamModel,
    LinearStaticAnalysis,
    LinearTriangle,
    LinearTriangleModel,
    Node,
    Spring,
    SpringModel,
    StaticResult,
    Truss,
    TrussModel,
    solve,
)


ROOT = Path(__file__).resolve().parents[1]
EXAMPLES = ROOT / "examples"

from nusa import Material, Section

def _make_bar(nodes, E, A):
    return Bar(nodes, material=Material(E=E), section=Section(A=A))

def _make_truss(nodes, E, A):
    return Truss(nodes, material=Material(E=E), section=Section(A=A))



def _spring_case():
    model = SpringModel("spring acceptance")
    n1, n2 = Node((0.0, 0.0)), Node((0.0, 0.0))
    element = Spring((n1, n2), 100.0)
    model.add_nodes([n1, n2])
    model.add_element(element)
    model.add_constraint(n1, ux=0.0)
    model.add_force(n2, (50.0,))
    return model, n2, element


def _bar_case():
    model = BarModel("bar acceptance")
    n1, n2 = Node((0.0, 0.0)), Node((2.0, 0.0))
    element = _make_bar((n1, n2), E=100.0, A=2.0)
    model.add_nodes([n1, n2])
    model.add_element(element)
    model.add_constraint(n1, ux=0.0)
    model.add_force(n2, (10.0,))
    return model, n2, element


def _truss_case():
    model = TrussModel("truss acceptance")
    n1, n2 = Node((0.0, 0.0)), Node((2.0, 0.0))
    element = _make_truss((n1, n2), E=100.0, A=2.0)
    model.add_nodes([n1, n2])
    model.add_element(element)
    model.add_constraint(n1, ux=0.0, uy=0.0)
    model.add_constraint(n2, uy=0.0)
    model.add_force(n2, (10.0, 0.0))
    return model, n2, element


def _beam_case():
    model = BeamModel("beam acceptance")
    n1, n2 = Node((0.0, 0.0)), Node((1.0, 0.0))
    element = Beam((n1, n2), E=1.0, I=1.0)
    model.add_nodes([n1, n2])
    model.add_element(element)
    model.add_constraint(n1, uy=0.0, ur=0.0)
    model.add_force(n2, (-1.0,))
    return model, n2, element


def _triangle_case():
    model = LinearTriangleModel("triangle acceptance")
    n1 = Node((0.0, 0.0))
    n2 = Node((1.0, 0.5))
    n3 = Node((0.0, 1.0))
    element = LinearTriangle((n1, n2, n3), E=200e9, nu=0.3, t=0.1)
    model.add_nodes([n1, n2, n3])
    model.add_element(element)
    model.add_constraint(n1, ux=0.0, uy=0.0)
    model.add_constraint(n3, ux=0.0, uy=0.0)
    model.add_force(n2, (1000.0, 0.0))
    return model, n2, element


@pytest.mark.parametrize(
    ("builder", "expected_displacement", "result_key"),
    [
        (_spring_case, 0.5, "force_j"),
        (_bar_case, 0.1, "axial_force"),
        (_truss_case, 0.1, "axial_force"),
        (_beam_case, -1.0 / 3.0, "shear_force_j"),
        (_triangle_case, 9.1e-8, "stress_xx"),
    ],
)
def test_all_public_families_complete_model_analysis_result_flow(
    builder,
    expected_displacement,
    result_key,
):
    model, node, element = builder()

    result = solve(model)

    assert isinstance(result, StaticResult)
    assert result.model_name == model.name
    assert result_key in result.element_result(element)

    first_component = model.displacement_dofs[0]
    assert result.displacement(node)[first_component] == pytest.approx(
        expected_displacement
    )

    for state_name in (
        "displacements",
        "nodal_forces",
        "reactions",
        "element_results",
        "stiffness_matrix",
        "_last_result",
    ):
        assert not hasattr(model, state_name)

    for state_name in ("ux", "uy", "ur", "fx", "fy", "m"):
        assert not hasattr(node, state_name)


def test_model_solve_and_explicit_analysis_are_equivalent():
    model_a, _, _ = _spring_case()
    model_b, _, _ = _spring_case()

    delegated = model_a.solve()
    explicit = LinearStaticAnalysis().solve(model_b)

    np.testing.assert_allclose(delegated.displacements, explicit.displacements)
    np.testing.assert_allclose(delegated.reactions, explicit.reactions)
    assert delegated.element_results == explicit.element_results


def test_result_snapshot_survives_later_problem_mutation():
    model, node, element = _spring_case()
    result = model.solve()

    original_displacements = result.displacements
    original_element = result.element_result(element)

    model.add_force(node, (80.0,))
    node.coordinates[:] = (5.0, 6.0)

    np.testing.assert_allclose(result.displacements, original_displacements)
    assert result.element_result(element) == original_element
    np.testing.assert_allclose(result.node_coordinates[1], [0.0, 0.0])


def test_removed_0_3_solved_state_api_is_absent():
    model, node, element = _spring_case()

    for name in (
        "assemble",
        "stiffness_matrix",
        "simple_report",
        "plot_deformed_shape",
        "element_result",
        "reaction",
    ):
        assert not hasattr(model, name)

    for name in ("ux", "fx", "get_displacements", "get_forces"):
        assert not hasattr(node, name)

    for name in ("fx", "sx", "_result_values"):
        assert not hasattr(element, name)


def test_examples_do_not_reference_removed_0_3_solved_state_api():
    forbidden = (
        re.compile(r"\.stiffness_matrix\b"),
        re.compile(r"\b(?:model|m\d*|ms)\.simple_report\s*\("),
        re.compile(r"\b(?:model|m\d*|ms)\.plot_(?:deformed_shape|nodal_result|element_result|moment_diagram|shear_diagram)\s*\("),
        re.compile(r"\b(?:n\d+|node|nodos\[[^\]]+\])\.(?:ux|uy|ur|fx|fy|m)\b"),
        re.compile(r"\b(?:e\d+|element)\.(?:fx|fy|sx|sy|sxy|m)\b"),
    )

    offenders = []
    for path in EXAMPLES.rglob("*.py"):
        text = path.read_text(encoding="utf-8")
        for pattern in forbidden:
            if pattern.search(text):
                offenders.append((path.relative_to(ROOT), pattern.pattern))

    assert offenders == []
