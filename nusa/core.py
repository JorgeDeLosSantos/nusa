# ***********************************
#  Author: Pedro Jorge De Los Santos    
#  E-mail: delossantosmfq@gmail.com 
#  License: MIT License
# ***********************************
import numpy as np

#~ ===========================  MODEL  ===========================
class Model:
    """
    Base class for all Finite Element Analysis (FEA) models.
    This class provides a base container for nodes and elements, enabling derived models to construct and manipulate FEA structures.
    """
    def __init__(self,name,mtype):
        """
        Initialize a new FEA model.

        Parameters
        ----------
        name : str
            Name of the model.
        mtype : str
            Type of model (e.g., 'bar', 'truss', 'beam').
        """
        self.mtype = mtype # Model type
        self.name = name # Name 
        self._nodes = [] # Nodes in model insertion order
        self._node_index = {} # Node object -> contiguous internal solver index
        self._elements = {} # Dictionary for elements {number: ElementObject}
        self._applied_forces = {} # Node -> explicitly applied nodal loads
        self._prescribed_displacements = {} # Node -> explicitly prescribed DOFs
        
    def add_node(self,node):
        """
        Add a node to the model.

        Parameters
        ----------
        node : :class:`~nusa.core.Node`
            Instance of a Node to be added.

        Returns
        -------
        None
        """
        if not isinstance(node, Node):
            raise TypeError("Model nodes must be Node instances")

        labels = [current.label for current in self._nodes]
        if node.label is None:
            label = 0
            while label in labels:
                label += 1
            node.label = label
        elif node.label in labels:
            raise ValueError(
                f"Node label {node.label!r} already exists in this model"
            )

        self._node_index[node] = len(self._nodes)
        self._nodes.append(node)

    def add_nodes(self, nodes):
        """
        Add multiple nodes to the model.

        Parameters
        ----------  

        nodes : list
            List of Node instances to be added.
        """
        for node in nodes:
            self.add_node(node)
        
    def add_element(self,element):
        """
        Add an element to the model.

        Parameters
        ----------
        element : :class:`~nusa.core.Element`
            Instance of an Element to be added.

        Raises
        ------
        ValueError
            If the element type does not match the model type.

        Example
        -------
        >>> m1 = BarModel()
        >>> E, A = 200e9, 0.001
        >>> n1 = Node((0,0))
        >>> n2 = Node((1,0))
        >>> e1 = Bar((n1,n2), E, A)
        >>> m1.add_element(e1)
        """

        if not isinstance(element, Element):
            raise TypeError("Model elements must be Element instances")

        if element.etype != self.mtype:
            raise ValueError(
                f"Element type '{element.etype}' incompatible with model '{self.mtype}'"
            )

        if element in self._elements.values():
            raise ValueError("Element already belongs to this model")

        missing_nodes = [node for node in element.nodes if node not in self._node_index]
        if missing_nodes:
            raise ValueError(
                "Element references nodes that do not belong to this model"
            )

        labels = set(self._elements)
        if element.label is None:
            label = 0
            while label in labels:
                label += 1
            element.label = label
        elif element.label in labels:
            raise ValueError(
                f"Element label {element.label!r} already exists in this model"
            )

        self._elements[element.label] = element

    def add_elements(self, elements):
        """
        Add multiple elements to the model.

        Parameters
        ----------
        elements : list
            List of Element instances to be added.
        """
        for element in elements:
            self.add_element(element)

    def _validate_topology(self):
        """Validate structural connectivity before global assembly."""
        if not self._elements:
            raise ValueError("Cannot assemble a model without elements")

        connected_nodes = {
            node
            for element in self.elements
            for node in element.nodes
        }
        orphan_nodes = [
            node for node in self.nodes
            if node not in connected_nodes
        ]
        if orphan_nodes:
            labels = [node.label for node in orphan_nodes]
            raise ValueError(
                "Model contains nodes not connected to any element: "
                f"{labels}"
            )

    @property
    def nodes(self):
        """
        Return a list of node objects.

        Returns
        -------
        list
            List of Node instances.
        """
        return list(self._nodes)

    @property
    def n_nodes(self):
        """
        Return the number of nodes in the model.

        Returns
        -------
        int
            Total number of nodes.
        """
        return len(self._nodes)

    def _get_node_index(self, node):
        """Return the model-owned contiguous index for a node."""
        try:
            return self._node_index[node]
        except KeyError:
            raise ValueError("Node does not belong to this model")

    def _global_dof_index(self, node, variable, dof_names):
        """Return the global vector index for one nodal degree of freedom."""
        try:
            component = dof_names.index(variable)
        except ValueError:
            raise ValueError(
                f"Unknown degree of freedom {variable!r}; expected one of {dof_names}"
            )
        return self.dof * self._get_node_index(node) + component


    def _validated_component_vector(self, values, names, quantity):
        """Return finite numeric components in the declared model order."""
        try:
            array = np.asarray(values, dtype=float)
        except (TypeError, ValueError) as exc:
            raise ValueError(
                f"{quantity} must contain {len(names)} finite numeric component(s)"
            ) from exc

        if array.ndim == 0:
            array = array.reshape(1)
        else:
            array = array.reshape(-1)

        if array.size != len(names):
            raise ValueError(
                f"{quantity} requires exactly {len(names)} component(s) "
                f"{names}; got {array.size}"
            )
        if not np.isfinite(array).all():
            raise ValueError(f"{quantity} components must be finite")

        return {
            name: float(array[index])
            for index, name in enumerate(names)
        }

    def _validated_named_components(self, values, names, quantity):
        """Validate finite named components against the model's active DOFs."""
        unknown = set(values) - set(names)
        if unknown:
            unknown_names = ", ".join(sorted(unknown))
            expected = ", ".join(names)
            raise ValueError(
                f"Unsupported {quantity} component(s): {unknown_names}; "
                f"expected only: {expected}"
            )

        validated = {}
        for name, value in values.items():
            try:
                scalar = float(value)
            except (TypeError, ValueError) as exc:
                raise ValueError(
                    f"{quantity} component {name!r} must be a finite scalar"
                ) from exc
            if not np.isfinite(scalar):
                raise ValueError(
                    f"{quantity} component {name!r} must be a finite scalar"
                )
            validated[name] = scalar
        return validated

    def _record_applied_forces(self, node, **values):
        """Persist explicitly applied nodal loads."""
        self._get_node_index(node)
        self._applied_forces.setdefault(node, {}).update(values)

    def _record_prescribed_displacements(self, node, **values):
        """Persist explicitly prescribed nodal degrees of freedom."""
        self._get_node_index(node)
        self._prescribed_displacements.setdefault(node, {}).update(values)

    @property
    def applied_loads(self):
        """Return the global vector of explicitly applied nodal loads.

        Components follow node insertion order and ``force_dofs``.
        """
        vector = np.zeros(self.dof * self.n_nodes, dtype=float)
        for node, values in self._applied_forces.items():
            if node not in self._node_index:
                continue
            for variable, value in values.items():
                if variable in self.force_dofs:
                    index = self._global_dof_index(node, variable, self.force_dofs)
                    vector[index] = value
        return vector

    @property
    def prescribed_displacements(self):
        """Return the global prescribed-displacement vector.

        Free degrees of freedom are represented by ``numpy.nan``. Components
        follow node insertion order and ``displacement_dofs``.
        """
        vector = np.full(self.dof * self.n_nodes, np.nan, dtype=float)
        for node, values in self._prescribed_displacements.items():
            if node not in self._node_index:
                continue
            for variable, value in values.items():
                if variable in self.displacement_dofs:
                    index = self._global_dof_index(
                        node, variable, self.displacement_dofs
                    )
                    vector[index] = value
        return vector

    def solve(self):
        """Run a linear-static analysis and return a new StaticResult.

        The model remains a problem definition: solving does not write
        displacements, forces, reactions, matrices, or result caches back to
        the model or its nodes.
        """
        from .analysis import LinearStaticAnalysis

        return LinearStaticAnalysis().solve(self)

    def applied_load(self, node):
        """Return explicitly applied load components for one node."""
        self._get_node_index(node)
        values = self._applied_forces.get(node, {})
        return {name: values.get(name, 0.0) for name in self.force_dofs}

    def prescribed_displacement(self, node):
        """Return prescribed displacement components for one node.

        Unprescribed degrees of freedom are returned as numpy.nan.
        """
        self._get_node_index(node)
        values = self._prescribed_displacements.get(node, {})
        return {
            name: values.get(name, np.nan)
            for name in self.displacement_dofs
        }

    @property
    def elements(self):
        """
        Return a list of element objects.

        Returns
        -------
        list
            List of Element instances.
        """
        return list(self._elements.values())

    @property
    def n_elements(self):
        """
        Return the number of elements in the model.

        Returns
        -------
        int
            Total number of elements.
        """
        return len(self._elements)

    
    def __str__(self):
        """
        Return a string representation of the model.

        Returns
        -------
        str
            Model name and number of nodes/elements.
        """
        return (
            f"Model: {self.name}\n"
            f"Nodes: {self.n_nodes}\n"
            f"Elements: {self.n_elements}"
        )
    
    def __repr__(self):
        """
        Return a string representation of the model.

        Returns
        -------
        str
            Model name and number of nodes/elements.
        """
        return (
            f"Model: {self.name}\n"
            f"Nodes: {self.n_nodes}\n"
            f"Elements: {self.n_elements}"
        )




#~ =========================== ELEMENT ===========================

class Element:
    """Base class for finite elements.

    Elements own formulation, connectivity, and physical properties. Solved
    response belongs to analysis results, not to the element instance.
    """

    def __init__(self, etype):
        self.etype = etype
        self.label = None

    def __str__(self):
        return str(self.__class__)


#~ =========================== NODE ===========================

class Node:
    """Geometry/topology point used by finite-element models."""

    __slots__ = ("coordinates", "_label")

    def __init__(self, coordinates):
        try:
            coordinates = np.asarray(coordinates, dtype=float)
        except (TypeError, ValueError) as exc:
            raise ValueError(
                "Node coordinates must contain two finite numbers"
            ) from exc
        if coordinates.shape != (2,):
            raise ValueError("Node coordinates must contain exactly two values")
        if not np.isfinite(coordinates).all():
            raise ValueError("Node coordinates must be finite")

        self.coordinates = coordinates.copy()
        self._label = None

    @property
    def x(self):
        return self.coordinates[0]

    @property
    def y(self):
        return self.coordinates[1]

    @property
    def label(self):
        return self._label

    @label.setter
    def label(self, value):
        self._label = value

    def __str__(self):
        return f"Node {self.label}: ({self.x},{self.y})"

    def __repr__(self):
        return f"<Node {self.label}: ({self.x},{self.y})>"


if __name__=='__main__':
    pass
