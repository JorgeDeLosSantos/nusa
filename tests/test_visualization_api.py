"""Regression tests for result-based visualization API."""

import matplotlib.pyplot as plt
import numpy as np
import pytest

from nusa import (
    Beam,
    BeamModel,
    LinearTriangleModel,
    Node,
    StaticResult,
    TrussModel,
    plot_deformed_shape,
)


def test_legacy_model_visualization_names_remain_as_compatibility_wrappers():
    truss = TrussModel()
    beam = BeamModel()
    triangle = LinearTriangleModel()

    assert hasattr(truss, "plot_deformed_shape")
    assert hasattr(beam, "plot_deformed_shape")
    assert hasattr(beam, "plot_moment_diagram")
    assert hasattr(beam, "plot_shear_diagram")
    assert hasattr(triangle, "plot_nodal_result")
    assert hasattr(triangle, "plot_element_result")

    assert not hasattr(beam, "plot_disp")
    assert not hasattr(triangle, "plot_nsol")
    assert not hasattr(triangle, "plot_esol")


def test_solution_visualization_requires_solved_model_for_legacy_wrapper():
    model = BeamModel()
    with pytest.raises(RuntimeError, match="after solve"):
        model.plot_deformed_shape()


def test_beam_model_deformed_shape_delegates_to_static_result():
    model = BeamModel("Scaled deformation")
    n1 = Node((0.0, 0.0))
    n2 = Node((1.0, 0.0))
    model.add_nodes([n1, n2])
    model.add_element(Beam((n1, n2), E=1.0, I=1.0))
    model.add_constraint(n1, uy=0.0, ur=0.0)
    model.add_force(n2, (-0.6,))

    result = model.solve()
    assert isinstance(result, StaticResult)

    ax = model.plot_deformed_shape(scale=10.0)
    assert ax is not None
    plt.close("all")


def test_top_level_deformed_plot_consumes_static_result():
    model = BeamModel("Top-level deformation")
    n1 = Node((0.0, 0.0))
    n2 = Node((1.0, 0.0))
    model.add_nodes([n1, n2])
    model.add_element(Beam((n1, n2), E=1.0, I=1.0))
    model.add_constraint(n1, uy=0.0, ur=0.0)
    model.add_force(n2, (-0.6,))

    result = model.solve()
    ax = plot_deformed_shape(result, scale=2.0)

    assert len(ax.lines) == 2
    undeformed_x, undeformed_y = ax.lines[0].get_data()
    deformed_x, deformed_y = ax.lines[1].get_data()
    np.testing.assert_allclose(undeformed_x, [0.0, 1.0])
    np.testing.assert_allclose(undeformed_y, [0.0, 0.0])
    np.testing.assert_allclose(deformed_x, [0.0, 1.0])
    np.testing.assert_allclose(
        deformed_y,
        [0.0, 2.0 * result.displacements[2]],
    )
    plt.close("all")
