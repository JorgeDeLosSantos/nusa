"""Post-processing helpers for NuSA analysis results."""

from __future__ import annotations

import numpy as np

from .result import StaticResult


def _require_result(result):
    if not isinstance(result, StaticResult):
        raise TypeError("post-processing requires a StaticResult")


def element_field(result, name):
    """Return one canonical scalar element field in result order."""
    _require_result(result)
    records = result.element_results
    if not records:
        return np.empty(0, dtype=float)
    if any(name not in record for record in records):
        raise ValueError(f"Element field {name!r} is not available")
    return np.array([record[name] for record in records], dtype=float)


def _displacement_component(result, name):
    try:
        component = result.displacement_dofs.index(name)
    except ValueError:
        return None
    width = len(result.displacement_dofs)
    values = result.displacements.reshape(-1, width)
    return values[:, component].copy()


def nodal_field(result, name, recovery="average"):
    """Return a nodal field from primary or recovered result quantities.

    Continuum element fields are recovered to nodes using an arithmetic
    average of adjacent element values. This is intentionally the first,
    explicit recovery policy; richer recovery methods can be added later.
    """
    _require_result(result)

    direct = _displacement_component(result, name)
    if direct is not None:
        return direct

    if name == "displacement_magnitude":
        ux = _displacement_component(result, "ux")
        uy = _displacement_component(result, "uy")
        if ux is None and uy is None:
            raise ValueError("Displacement magnitude is not available")
        if ux is None:
            ux = np.zeros_like(uy)
        if uy is None:
            uy = np.zeros_like(ux)
        return np.sqrt(ux**2 + uy**2)

    if name == "von_mises_stress":
        sxx = nodal_field(result, "stress_xx", recovery=recovery)
        syy = nodal_field(result, "stress_yy", recovery=recovery)
        sxy = nodal_field(result, "stress_xy", recovery=recovery)
        return np.sqrt(sxx**2 - sxx * syy + syy**2 + 3.0 * sxy**2)

    if recovery != "average":
        raise ValueError(
            f"Unknown recovery method {recovery!r}; expected 'average'"
        )

    element_values = element_field(result, name)
    totals = np.zeros(len(result.node_labels), dtype=float)
    counts = np.zeros(len(result.node_labels), dtype=int)

    for value, connectivity in zip(element_values, result.connectivity):
        for node_index in connectivity:
            totals[node_index] += value
            counts[node_index] += 1

    if np.any(counts == 0):
        raise ValueError("Cannot recover a nodal field for unconnected nodes")

    return totals / counts
