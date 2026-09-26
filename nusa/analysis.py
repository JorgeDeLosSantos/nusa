"""Linear-static finite-element analysis orchestration."""

from __future__ import annotations

import numpy as np

from .core import Model
from .result import StaticResult


def _element_dof_indices(model, element):
    indices = []
    for node in element.nodes:
        node_index = model._get_node_index(node)
        base = model.dof * node_index
        indices.extend(base + component for component in range(model.dof))
    return indices


def _assemble_stiffness(model):
    model._validate_topology()
    size = model.dof * model.n_nodes
    stiffness = np.zeros((size, size), dtype=float)

    for element in model.elements:
        element_stiffness = np.asarray(
            element.get_element_stiffness(), dtype=float
        )
        indices = _element_dof_indices(model, element)
        stiffness[np.ix_(indices, indices)] += element_stiffness

    return stiffness


def _partition_system(stiffness, loads, prescribed):
    prescribed_dofs = np.flatnonzero(~np.isnan(prescribed))
    free_dofs = np.flatnonzero(np.isnan(prescribed))

    reduced = stiffness[np.ix_(free_dofs, free_dofs)]
    rhs = loads[free_dofs].copy()

    if prescribed_dofs.size:
        coupling = stiffness[np.ix_(free_dofs, prescribed_dofs)]
        rhs -= coupling @ prescribed[prescribed_dofs]

    return prescribed_dofs, free_dofs, reduced, rhs


class LinearStaticAnalysis:
    """Perform a small-displacement linear-static finite-element analysis."""

    def solve(self, model):
        """Solve *model* and return an independent :class:`StaticResult`."""
        if not isinstance(model, Model):
            raise TypeError("LinearStaticAnalysis.solve() requires a Model")

        nodes = tuple(model.nodes)
        elements = tuple(model.elements)
        stiffness = _assemble_stiffness(model)
        loads = np.asarray(model.applied_loads, dtype=float).copy()
        prescribed = np.asarray(
            model.prescribed_displacements, dtype=float
        ).copy()

        (
            prescribed_dofs,
            free_dofs,
            reduced,
            rhs,
        ) = _partition_system(stiffness, loads, prescribed)

        if reduced.size:
            if np.linalg.matrix_rank(reduced) < reduced.shape[0]:
                raise np.linalg.LinAlgError(
                    "Singular stiffness matrix: the model may be "
                    "underconstrained or contain a mechanism."
                )
            free_displacements = np.linalg.solve(reduced, rhs)
        else:
            free_displacements = np.empty(0, dtype=float)

        displacements = prescribed.copy()
        displacements[free_dofs] = free_displacements

        nodal_forces = stiffness @ displacements
        reactions = np.zeros_like(nodal_forces)
        reactions[prescribed_dofs] = (
            nodal_forces[prescribed_dofs] - loads[prescribed_dofs]
        )

        node_indices = {node: index for index, node in enumerate(nodes)}
        connectivity = tuple(
            tuple(node_indices[node] for node in element.nodes)
            for element in elements
        )

        return StaticResult(
            model_name=model.name,
            node_objects=nodes,
            node_labels=tuple(node.label for node in nodes),
            node_coordinates=np.array(
                [[node.x, node.y] for node in nodes], dtype=float
            ),
            element_objects=elements,
            element_labels=tuple(element.label for element in elements),
            element_types=tuple(element.etype for element in elements),
            connectivity=connectivity,
            displacement_dofs=model.displacement_dofs,
            force_dofs=model.force_dofs,
            applied_loads=loads,
            prescribed_displacements=prescribed,
            displacements=displacements,
            nodal_forces=nodal_forces,
            reactions=reactions,
        )


def solve(model):
    """Solve *model* using the default linear-static analysis."""
    return LinearStaticAnalysis().solve(model)
