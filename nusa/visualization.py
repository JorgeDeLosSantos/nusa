"""Visualization helpers for NuSA analysis results."""

from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from matplotlib.collections import PatchCollection
from matplotlib.patches import Polygon

from .core import Model
from .post import element_field, nodal_field
from .result import StaticResult


def _require_result(result):
    if not isinstance(result, StaticResult):
        raise TypeError("visualization requires a StaticResult")


def _axes(ax=None):
    if ax is not None:
        return ax
    _, ax = plt.subplots()
    return ax


def _component_matrix(result):
    width = len(result.displacement_dofs)
    return result.displacements.reshape(-1, width)


def _component(result, name):
    try:
        index = result.displacement_dofs.index(name)
    except ValueError:
        return np.zeros(len(result.node_labels), dtype=float)
    return _component_matrix(result)[:, index]


def _set_equal_limits(ax, coordinates):
    if coordinates.size == 0:
        return
    xmin, ymin = coordinates.min(axis=0)
    xmax, ymax = coordinates.max(axis=0)
    span_x = xmax - xmin
    span_y = ymax - ymin
    margin_x = span_x / 7.0 if span_x else 1.0 / 7.0
    margin_y = span_y / 7.0 if span_y else 1.0 / 7.0
    ax.set_xlim(xmin - margin_x, xmax + margin_x)
    ax.set_ylim(ymin - margin_y, ymax + margin_y)
    ax.set_aspect("equal")




def _require_model(model):
    if not isinstance(model, Model):
        raise TypeError("plot_model() requires a Model")


def _problem_coordinates(model):
    return np.asarray([[node.x, node.y] for node in model.nodes], dtype=float)


def _arrow_size(coordinates):
    if coordinates.size == 0:
        return 1.0
    xmin, ymin = coordinates.min(axis=0)
    xmax, ymax = coordinates.max(axis=0)
    span = max(xmax - xmin, ymax - ymin)
    return 0.08 * span if span else 0.08


def _draw_force_arrow(ax, x, y, axis, direction, size):
    if axis == "x":
        ax.arrow(
            x, y, direction * size, 0.0,
            head_width=size / 4.0,
            head_length=size / 3.0,
            length_includes_head=True,
        )
    elif axis == "y":
        ax.arrow(
            x, y, 0.0, direction * size,
            head_width=size / 4.0,
            head_length=size / 3.0,
            length_includes_head=True,
        )


def _draw_constraint_marker(ax, x, y, dof):
    marker = {"ux": "<", "uy": "v", "ur": "s"}.get(dof, "x")
    ax.plot(x, y, marker=marker, linestyle="None", markersize=9, alpha=0.7)


def plot_model(model, ax=None):
    """Plot model geometry, applied translational loads, and constraints."""
    _require_model(model)
    ax = _axes(ax)
    coordinates = _problem_coordinates(model)

    if model.mtype == "triangle":
        patches = [
            Polygon(
                [[node.x, node.y] for node in element.nodes],
                closed=True,
            )
            for element in model.elements
        ]
        collection = PatchCollection(patches, alpha=0.35, edgecolor="k")
        ax.add_collection(collection)
    else:
        for element in model.elements:
            xy = np.asarray([[node.x, node.y] for node in element.nodes], dtype=float)
            ax.plot(xy[:, 0], xy[:, 1], marker="o")

    size = _arrow_size(coordinates)
    for node in model.nodes:
        loads = model.applied_load(node)
        for component, axis in (("fx", "x"), ("fy", "y")):
            value = loads.get(component, 0.0)
            if value != 0.0:
                _draw_force_arrow(
                    ax,
                    node.x,
                    node.y,
                    axis,
                    1 if value > 0.0 else -1,
                    size,
                )

        prescribed = model.prescribed_displacement(node)
        for dof, value in prescribed.items():
            if np.isfinite(value):
                _draw_constraint_marker(ax, node.x, node.y, dof)

    _set_equal_limits(ax, coordinates)
    ax.set_title(model.name)
    return ax


def plot_deformed_shape(result, scale=1.0, ax=None, **kwargs):
    """Plot frozen undeformed and deformed geometry."""
    _require_result(result)
    ax = _axes(ax)
    coordinates = result.node_coordinates
    ux = _component(result, "ux")
    uy = _component(result, "uy")
    deformed = coordinates + scale * np.column_stack((ux, uy))

    undeformed_kwargs = {"marker": "o", "linestyle": "-", "alpha": 0.5}
    deformed_kwargs = {"marker": "o", "linestyle": "--"}
    deformed_kwargs.update(kwargs)

    for connectivity in result.connectivity:
        indices = list(connectivity)
        if len(indices) >= 3:
            indices = indices + [indices[0]]
        original = coordinates[indices]
        current = deformed[indices]
        ax.plot(original[:, 0], original[:, 1], **undeformed_kwargs)
        ax.plot(current[:, 0], current[:, 1], **deformed_kwargs)

    all_coordinates = np.vstack((coordinates, deformed))
    _set_equal_limits(ax, all_coordinates)
    return ax


def _triangle_connectivity(result):
    triangles = np.asarray(result.connectivity, dtype=int)
    if triangles.ndim != 2 or triangles.shape[1] != 3:
        raise ValueError("This visualization requires triangular connectivity")
    return triangles


def plot_nodal_field(result, field, ax=None, recovery="average"):
    """Plot a scalar nodal field over a triangular mesh."""
    _require_result(result)
    ax = _axes(ax)
    coordinates = result.node_coordinates
    triangles = _triangle_connectivity(result)
    values = nodal_field(result, field, recovery=recovery)

    triangulation = mtri.Triangulation(
        coordinates[:, 0],
        coordinates[:, 1],
        triangles=triangles,
    )
    contour = ax.tricontourf(triangulation, values, cmap="jet")
    ax.figure.colorbar(contour, ax=ax)
    _set_equal_limits(ax, coordinates)
    ax.set_title(
        f"{field} (Max:{values.max():0.3e}, Min:{values.min():0.3e})",
        fontsize=8,
    )
    return ax


def plot_element_field(result, field, ax=None):
    """Plot a scalar element field over a triangular mesh."""
    _require_result(result)
    ax = _axes(ax)
    coordinates = result.node_coordinates
    triangles = _triangle_connectivity(result)
    values = element_field(result, field)

    patches = [
        Polygon(coordinates[connectivity], closed=True)
        for connectivity in triangles
    ]
    collection = PatchCollection(patches, cmap="jet", alpha=1.0)
    collection.set_array(values)
    ax.add_collection(collection)
    ax.figure.colorbar(collection, ax=ax)
    _set_equal_limits(ax, coordinates)
    ax.set_title(
        f"{field} (Max:{values.max():0.3e}, Min:{values.min():0.3e})",
        fontsize=8,
    )
    return ax


def _beam_diagram_data(result, first_key, second_key, flip_first=False, flip_second=False):
    _require_result(result)
    if any(element_type != "beam" for element_type in result.element_types):
        raise ValueError("Beam diagrams require beam element results")

    coordinates = result.node_coordinates
    x_values = []
    field_values = []
    cursor = 0.0

    for connectivity, record in zip(result.connectivity, result.element_results):
        if len(connectivity) != 2:
            raise ValueError("Beam diagrams require two-node elements")
        i, j = connectivity
        length = float(np.linalg.norm(coordinates[j] - coordinates[i]))
        first = float(record[first_key])
        second = float(record[second_key])
        if flip_first:
            first = -first
        if flip_second:
            second = -second
        x_values.extend((cursor, cursor + length))
        field_values.extend((first, second))
        cursor += length

    return np.asarray(x_values), np.asarray(field_values)


def plot_moment_diagram(result, ax=None):
    """Plot the legacy beam end-moment diagram from frozen element actions."""
    ax = _axes(ax)
    x_values, moments = _beam_diagram_data(
        result,
        "bending_moment_i",
        "bending_moment_j",
        flip_first=True,
    )
    ax.plot(x_values, moments)
    ax.fill_between(x_values, moments)
    return ax


def plot_shear_diagram(result, ax=None):
    """Plot the legacy beam end-shear diagram from frozen element actions."""
    ax = _axes(ax)
    x_values, shear = _beam_diagram_data(
        result,
        "shear_force_i",
        "shear_force_j",
        flip_second=True,
    )
    ax.plot(x_values, shear)
    ax.fill_between(x_values, shear)
    return ax
