"""Text reporting for NuSA analysis results."""

from __future__ import annotations

import numpy as np
from tabulate import tabulate

from .result import StaticResult


_TABLE_OPTIONS = {
    "headers": "firstrow",
    "tablefmt": "rst",
    "numalign": "right",
}


def _node_vector_table(labels, vector, component_names):
    values = np.asarray(vector, dtype=float).reshape(len(labels), len(component_names))
    rows = [["Node"] + [name.upper() for name in component_names]]
    for label, components in zip(labels, values):
        rows.append([label] + list(components))
    return tabulate(rows, **_TABLE_OPTIONS)


def _element_results_table(result):
    records = result.element_results
    keys = tuple(records[0].keys()) if records else ()
    headers = ["Element"] + [name.replace("_", " ").upper() for name in keys]
    rows = [headers]
    for label, record in zip(result.element_labels, records):
        rows.append([label] + [record[key] for key in keys])
    return tabulate(rows, **_TABLE_OPTIONS)


def _nodes_table(result):
    rows = [["Node", "X", "Y"]]
    for label, coordinates in zip(result.node_labels, result.node_coordinates):
        rows.append([label, coordinates[0], coordinates[1]])
    return tabulate(rows, **_TABLE_OPTIONS)


def _elements_table(result):
    max_nodes = max((len(indices) for indices in result.connectivity), default=0)
    rows = [["Element"] + [f"N{k + 1}" for k in range(max_nodes)]]
    node_labels = result.node_labels

    for label, indices in zip(result.element_labels, result.connectivity):
        connected_labels = [node_labels[index] for index in indices]
        rows.append(
            [label]
            + connected_labels
            + [""] * (max_nodes - len(connected_labels))
        )

    return tabulate(rows, **_TABLE_OPTIONS)


def simple_report(result, report_type="print", fname="nusa_rpt.txt"):
    """Generate a compact text report from a :class:`StaticResult`."""
    if not isinstance(result, StaticResult):
        raise TypeError("simple_report() requires a StaticResult")

    valid_report_types = {"print", "string", "write"}
    if report_type not in valid_report_types:
        raise ValueError(
            f"Unknown report_type {report_type!r}; "
            f"expected one of {sorted(valid_report_types)}"
        )

    sections = [
        "==========================",
        "    NuSA Simple Report",
        "==========================",
        "",
        f"Model: {result.model_name}",
        f"Number of nodes: {len(result.node_labels)}",
        f"Number of elements: {len(result.element_labels)}",
        "",
        "RESULTS",
        "",
        "NODAL DISPLACEMENTS",
        _node_vector_table(
            result.node_labels,
            result.displacements,
            result.displacement_dofs,
        ),
        "",
        "APPLIED LOADS",
        _node_vector_table(
            result.node_labels,
            result.applied_loads,
            result.force_dofs,
        ),
        "",
        "NODAL FORCES (K @ U)",
        _node_vector_table(
            result.node_labels,
            result.nodal_forces,
            result.force_dofs,
        ),
        "",
        "REACTIONS",
        _node_vector_table(
            result.node_labels,
            result.reactions,
            result.force_dofs,
        ),
        "",
        "ELEMENT RESULTS",
        _element_results_table(result),
        "",
        "FINITE ELEMENT MODEL INFO",
        "",
        "NODES",
        _nodes_table(result),
        "",
        "ELEMENTS",
        _elements_table(result),
    ]
    report = "\n".join(sections) + "\n"

    if report_type == "print":
        print(report)
        return None
    if report_type == "write":
        with open(fname, "w", encoding="utf-8") as report_file:
            report_file.write(report)
        return None
    return report
