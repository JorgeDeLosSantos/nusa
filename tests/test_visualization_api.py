"""Regression tests for result-owned visualization API."""

import matplotlib.pyplot as plt
import numpy as np

from nusa import (
    Beam,
    BeamModel,
    LinearTriangleModel,
    Node,
    StaticResult,
    TrussModel,
    plot_deformed_shape,
    plot_model,
)


def test_all_visualization_is_not_model_owned():
    truss = TrussModel()
    beam = BeamModel()
    triangle = LinearTriangleModel()

    for model, names in (
        (truss, ("plot_model", "plot_deformed_shape")),
        (beam, ("plot_model", "plot_deformed_shape", "plot_moment_diagram", "plot_shear_diagram")),
        (triangle, ("plot_model", "plot_nodal_result", "plot_element_result")),
    ):
        for name in names:
            assert not hasattr(model, name)


def _beam_result():
    model = BeamModel("deformation")
    n1 = Node((0.0, 0.0))
    n2 = Node((1.0, 0.0))
    model.add_nodes([n1, n2])
    model.add_element(Beam((n1, n2), E=1.0, I=1.0))
    model.add_constraint(n1, uy=0.0, ur=0.0)
    model.add_force(n2, (-0.6,))
    return model.solve()


def test_static_result_exposes_visualization_convenience():
    result = _beam_result()

    assert isinstance(result, StaticResult)
    ax = result.plot_deformed_shape(scale=10.0)
    assert ax is not None
    plt.close("all")


def test_top_level_deformed_plot_consumes_static_result():
    result = _beam_result()
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


def test_top_level_plot_model_consumes_problem_definition():
    model = BeamModel("problem plot")
    n1 = Node((0.0, 0.0))
    n2 = Node((1.0, 0.0))
    model.add_nodes([n1, n2])
    model.add_element(Beam((n1, n2), E=1.0, I=1.0))
    model.add_constraint(n1, uy=0.0, ur=0.0)
    model.add_force(n2, (-1.0,))

    ax = plot_model(model)

    assert ax is not None
    assert len(ax.lines) >= 1
    assert len(ax.patches) >= 1
    plt.close("all")
