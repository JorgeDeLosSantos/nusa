"""Contract tests for result-based post-processing and visualization."""

import matplotlib.pyplot as plt
import numpy as np
import pytest

from nusa import (
    Beam,
    BeamModel,
    LinearTriangle,
    LinearTriangleModel,
    Node,
    element_field,
    nodal_field,
    plot_element_field,
    plot_moment_diagram,
    plot_nodal_field,
    plot_shear_diagram,
    solve,
)


def _triangle_result():
    model = LinearTriangleModel("triangle fields")
    n1 = Node((0.0, 0.0))
    n2 = Node((1.0, 0.5))
    n3 = Node((0.0, 1.0))
    element = LinearTriangle((n1, n2, n3), E=200e9, nu=0.3, t=0.1)
    model.add_nodes([n1, n2, n3])
    model.add_element(element)
    model.add_constraint(n1, ux=0.0, uy=0.0)
    model.add_constraint(n3, ux=0.0, uy=0.0)
    model.add_force(n2, (1000.0, 0.0))
    return model, element, solve(model)


def _beam_result():
    model = BeamModel("beam diagrams")
    n1 = Node((0.0, 0.0))
    n2 = Node((2.0, 0.0))
    element = Beam((n1, n2), E=100.0, I=1.0)
    model.add_nodes([n1, n2])
    model.add_element(element)
    model.add_constraint(n1, uy=0.0, ur=0.0)
    model.add_force(n2, (-10.0,))
    return model, element, solve(model)


def test_element_field_uses_canonical_result_values_and_legacy_aliases():
    _, _, result = _triangle_result()

    np.testing.assert_allclose(element_field(result, "stress_xx"), [20000.0])
    np.testing.assert_allclose(element_field(result, "sxx"), [20000.0])
    np.testing.assert_allclose(result.element_field("strain_xx"), [9.1e-8])


def test_nodal_field_recovers_cst_values_by_arithmetic_average():
    _, _, result = _triangle_result()

    np.testing.assert_allclose(
        nodal_field(result, "stress_xx"),
        [20000.0, 20000.0, 20000.0],
    )
    np.testing.assert_allclose(
        result.nodal_field("sxx"),
        [20000.0, 20000.0, 20000.0],
    )
    np.testing.assert_allclose(
        result.nodal_field("ux"),
        [0.0, 9.1e-8, 0.0],
        atol=1e-15,
    )


def test_nodal_derived_fields_are_computed_from_result_snapshot():
    _, _, result = _triangle_result()

    magnitude = result.nodal_field("usum")
    np.testing.assert_allclose(magnitude, [0.0, 9.1e-8, 0.0], atol=1e-15)

    expected_vm = np.sqrt(20000.0**2 - 20000.0 * 6000.0 + 6000.0**2)
    np.testing.assert_allclose(
        result.nodal_field("seqv"),
        [expected_vm, expected_vm, expected_vm],
    )


def test_unknown_fields_and_recovery_methods_raise_clear_errors():
    _, _, result = _triangle_result()

    with pytest.raises(ValueError, match="not available"):
        element_field(result, "does_not_exist")
    with pytest.raises(ValueError, match="recovery"):
        nodal_field(result, "stress_xx", recovery="area_weighted")


def test_triangle_field_plots_consume_static_result():
    _, _, result = _triangle_result()

    ax_element = plot_element_field(result, "sxx")
    assert len(ax_element.collections) >= 1
    plt.close("all")

    ax_nodal = plot_nodal_field(result, "ux")
    assert len(ax_nodal.collections) >= 1
    plt.close("all")


def test_beam_diagrams_use_frozen_element_actions():
    model, element, result = _beam_result()

    expected = result.element_result(element)

    # Mutate the problem after the result exists. The diagrams must continue
    # to use the frozen element actions from the original result.
    model.add_force(model.nodes[-1], (-20.0,))

    ax_m = plot_moment_diagram(result)
    _, moment_values = ax_m.lines[0].get_data()
    np.testing.assert_allclose(
        moment_values,
        [-expected["bending_moment_i"], expected["bending_moment_j"]],
        atol=1e-12,
    )
    plt.close("all")

    ax_v = plot_shear_diagram(result)
    _, shear_values = ax_v.lines[0].get_data()
    np.testing.assert_allclose(
        shear_values,
        [expected["shear_force_i"], -expected["shear_force_j"]],
        atol=1e-12,
    )
    plt.close("all")


def test_result_visualization_remains_stable_after_geometry_mutation():
    model, _, result = _triangle_result()
    frozen = result.node_coordinates

    for node in model.nodes:
        node.coordinates[:] += 100.0

    ax = result.plot_deformed_shape(scale=1.0)

    first_original_x, first_original_y = ax.lines[0].get_data()
    first_connectivity = list(result.connectivity[0]) + [result.connectivity[0][0]]
    np.testing.assert_allclose(
        np.column_stack((first_original_x, first_original_y)),
        frozen[first_connectivity],
    )
    plt.close("all")
