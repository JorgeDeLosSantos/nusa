"""Regression tests for normalized visualization API names."""

import matplotlib.pyplot as plt
import numpy as np

from nusa.core import Node
from nusa.element import Beam
from nusa.model import BeamModel, LinearTriangleModel, TrussModel


def test_visualization_api_uses_descriptive_public_names():
    truss = TrussModel()
    beam = BeamModel()
    triangle = LinearTriangleModel()

    assert hasattr(truss, "plot_deformed_shape")
    assert hasattr(beam, "plot_deformed_shape")
    assert hasattr(triangle, "plot_nodal_result")
    assert hasattr(triangle, "plot_element_result")

    assert not hasattr(beam, "plot_disp")
    assert not hasattr(triangle, "plot_nsol")
    assert not hasattr(triangle, "plot_esol")


def test_geometry_and_scaling_helpers_are_private():
    truss = TrussModel()
    beam = BeamModel()
    triangle = LinearTriangleModel()

    for model in (truss, beam, triangle):
        assert hasattr(model, "_rect_region")
        assert not hasattr(model, "rect_region")

    assert hasattr(truss, "_calculate_deformed_factor")
    assert hasattr(triangle, "_calculate_deformed_factor")
    assert not hasattr(triangle, "calculate_deformed_factor")


def test_beam_plot_deformed_shape_uses_scale_keyword():
    model = BeamModel("Scaled deformation")
    n1 = Node((0.0, 0.0))
    n2 = Node((1.0, 0.0))
    model.add_nodes([n1, n2])
    model.add_element(Beam((n1, n2), E=1.0, I=1.0))

    n1.uy = 0.0
    n2.uy = 0.2

    model.plot_deformed_shape(scale=10.0)
    ax = plt.gcf().axes[0]
    xdata, ydata = ax.lines[0].get_data()

    np.testing.assert_allclose(xdata, [0.0, 1.0])
    np.testing.assert_allclose(ydata, [0.0, 2.0])
    plt.close("all")
