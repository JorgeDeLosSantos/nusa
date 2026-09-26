"""Result objects for NuSA finite-element analyses."""

from __future__ import annotations

import numpy as np


class StaticResult:
    """Stable snapshot of one completed linear-static analysis.

    Instances are produced by :func:\`nusa.solve\` or
    :class:\`nusa.LinearStaticAnalysis\`. Public array properties return
    copies so callers cannot mutate the stored result.
    """

    def __init__(
        self,
        *,
        model_name,
        node_objects,
        node_labels,
        node_coordinates,
        element_objects,
        element_labels,
        element_types,
        connectivity,
        displacement_dofs,
        force_dofs,
        applied_loads,
        prescribed_displacements,
        displacements,
        nodal_forces,
        reactions,
    ):
        self._model_name = str(model_name)
        self._node_labels = tuple(node_labels)
        self._node_coordinates = np.asarray(node_coordinates, dtype=float).copy()
        self._element_labels = tuple(element_labels)
        self._element_types = tuple(element_types)
        self._connectivity = tuple(tuple(indices) for indices in connectivity)
        self._displacement_dofs = tuple(displacement_dofs)
        self._force_dofs = tuple(force_dofs)

        self._applied_loads = np.asarray(applied_loads, dtype=float).copy()
        self._prescribed_displacements = np.asarray(
            prescribed_displacements, dtype=float
        ).copy()
        self._displacements = np.asarray(displacements, dtype=float).copy()
        self._nodal_forces = np.asarray(nodal_forces, dtype=float).copy()
        self._reactions = np.asarray(reactions, dtype=float).copy()

        self._node_index = {
            node: index for index, node in enumerate(tuple(node_objects))
        }
        self._element_index = {
            element: index for index, element in enumerate(tuple(element_objects))
        }

    @property
    def model_name(self):
        """Model name captured at solve time."""
        return self._model_name

    @property
    def node_labels(self):
        """Node labels in frozen result order."""
        return self._node_labels

    @property
    def node_coordinates(self):
        """Node coordinates in frozen result order."""
        return self._node_coordinates.copy()

    @property
    def element_labels(self):
        """Element labels in frozen result order."""
        return self._element_labels

    @property
    def element_types(self):
        """Element type names in frozen result order."""
        return self._element_types

    @property
    def connectivity(self):
        """Element connectivity as frozen result-node indices."""
        return self._connectivity

    @property
    def displacement_dofs(self):
        """Displacement component names in per-node order."""
        return self._displacement_dofs

    @property
    def force_dofs(self):
        """Generalized force component names in per-node order."""
        return self._force_dofs

    @property
    def applied_loads(self):
        """Applied nodal-load vector captured for this analysis."""
        return self._applied_loads.copy()

    @property
    def prescribed_displacements(self):
        """Prescribed-displacement vector captured for this analysis."""
        return self._prescribed_displacements.copy()

    @property
    def displacements(self):
        """Complete solved global displacement vector."""
        return self._displacements.copy()

    @property
    def nodal_forces(self):
        """Solved generalized nodal-force vector K @ u."""
        return self._nodal_forces.copy()

    @property
    def reactions(self):
        """Support-reaction vector for prescribed degrees of freedom."""
        return self._reactions.copy()

    def _node_position(self, node):
        try:
            return self._node_index[node]
        except KeyError as exc:
            raise ValueError("Node does not belong to this result") from exc

    def _node_components(self, node, vector, names):
        index = self._node_position(node)
        width = len(names)
        start = width * index
        return {
            name: float(vector[start + component])
            for component, name in enumerate(names)
        }

    def applied_load(self, node):
        """Return applied load components for one result node."""
        return self._node_components(node, self._applied_loads, self._force_dofs)

    def prescribed_displacement(self, node):
        """Return prescribed displacement components for one result node."""
        return self._node_components(
            node,
            self._prescribed_displacements,
            self._displacement_dofs,
        )

    def displacement(self, node):
        """Return solved displacement components for one result node."""
        return self._node_components(
            node,
            self._displacements,
            self._displacement_dofs,
        )

    def nodal_force(self, node):
        """Return solved generalized nodal-force components for one node."""
        return self._node_components(node, self._nodal_forces, self._force_dofs)

    def reaction(self, node):
        """Return support-reaction components for one node."""
        return self._node_components(node, self._reactions, self._force_dofs)

    def __repr__(self):
        return (
            f"StaticResult(model_name={self.model_name!r}, "
            f"nodes={len(self._node_labels)}, "
            f"elements={len(self._element_labels)})"
        )
