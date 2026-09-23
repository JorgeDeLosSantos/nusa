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
        self._is_assembled = False
        
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
        self._invalidate_assembly()

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

        for node in element.nodes:
            node._add_element(element)
        self._invalidate_assembly()

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
        """Persist explicitly applied nodal loads and invalidate solved state."""
        self._get_node_index(node)
        self._applied_forces.setdefault(node, {}).update(values)
        self._invalidate_solution()

    def _record_prescribed_displacements(self, node, **values):
        """Persist explicitly prescribed nodal DOFs and invalidate solved state."""
        self._get_node_index(node)
        self._prescribed_displacements.setdefault(node, {}).update(values)
        self._invalidate_solution()

    def _restore_input_state(self):
        """Restore explicit loads and prescribed DOFs into vectors and Node state."""
        for node, values in self._applied_forces.items():
            if node not in self._node_index:
                continue
            for variable, value in values.items():
                setattr(node, variable, value)
                if hasattr(self, "_f") and variable in self.force_dofs:
                    index = self._global_dof_index(node, variable, self.force_dofs)
                    self._f[index] = value

        for node, values in self._prescribed_displacements.items():
            if node not in self._node_index:
                continue
            for variable, value in values.items():
                setattr(node, variable, value)
                if hasattr(self, "_u") and variable in self.displacement_dofs:
                    index = self._global_dof_index(
                        node, variable, self.displacement_dofs
                    )
                    self._u[index] = value

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

    @property
    def displacements(self):
        """Return the solved global displacement vector.

        Results are available only after ``solve()``.
        """
        if not hasattr(self, "_nodal_forces"):
            raise RuntimeError("Displacements are available only after solve()")
        return self._u.copy()

    @property
    def nodal_forces(self):
        """Return the solved generalized nodal-force vector ``K @ u``."""
        if not hasattr(self, "_nodal_forces"):
            raise RuntimeError("Nodal forces are available only after solve()")
        return self._nodal_forces.copy()

    @property
    def reactions(self):
        """Return the solved global support-reaction vector.

        Entries are nonzero only at prescribed solver degrees of freedom and
        are computed as ``K @ u - applied_loads`` at those DOFs.
        """
        if not hasattr(self, "_reactions"):
            raise RuntimeError("Reactions are available only after solve()")
        return self._reactions.copy()

    def _get_node_vector_components(self, node, vector, names):
        """Return named components from a model-ordered global vector."""
        node_index = self._get_node_index(node)
        start = self.dof * node_index
        return {
            name: vector[start + component]
            for component, name in enumerate(names)
        }

    def applied_load(self, node):
        """Return explicitly applied load components for one node."""
        self._get_node_index(node)
        values = self._applied_forces.get(node, {})
        return {name: values.get(name, 0.0) for name in self.force_dofs}

    def prescribed_displacement(self, node):
        """Return prescribed displacement components for one node.

        Unprescribed degrees of freedom are returned as ``numpy.nan``.
        """
        self._get_node_index(node)
        values = self._prescribed_displacements.get(node, {})
        return {
            name: values.get(name, np.nan)
            for name in self.displacement_dofs
        }

    def displacement(self, node):
        """Return solved displacement components for one node."""
        return self._get_node_vector_components(
            node, self.displacements, self.displacement_dofs
        )

    def nodal_force(self, node):
        """Return solved generalized nodal-force components ``K @ u``."""
        return self._get_node_vector_components(
            node, self.nodal_forces, self.force_dofs
        )

    def reaction(self, node):
        """Return solved support-reaction components for one node."""
        return self._get_node_vector_components(
            node, self.reactions, self.force_dofs
        )

    def element_result(self, element):
        """Return normalized solved results for one model element."""
        if not hasattr(self, "_nodal_forces"):
            raise RuntimeError("Element results are available only after solve()")
        if element not in self._elements.values():
            raise ValueError("Element does not belong to this model")
        return {
            name: float(value)
            for name, value in element._result_values().items()
        }

    @property
    def element_results(self):
        """Return normalized solved results in element insertion order."""
        if not hasattr(self, "_nodal_forces"):
            raise RuntimeError("Element results are available only after solve()")
        return tuple(self.element_result(element) for element in self.elements)

    @property
    def stiffness_matrix(self):
        """Return a copy of the assembled global stiffness matrix."""
        if not self._is_assembled or not hasattr(self, "_K"):
            raise RuntimeError(
                "Stiffness matrix is available only after assemble() or solve()"
            )
        return self._K.copy()

    def _reset_input_vectors(self):
        """Rebuild numeric input vectors from persistent loads and constraints."""
        if not self._is_assembled:
            return
        matrix_size = self.dof * self.n_nodes
        self._f = np.zeros(matrix_size, dtype=float)
        self._u = np.full(matrix_size, np.nan, dtype=float)

    def _reset_node_state(self):
        """Clear solved nodal state and restore explicit model inputs."""
        for node in self._nodes:
            node.ux = np.nan
            node.uy = np.nan
            node.ur = np.nan
            node.fx = 0.0
            node.fy = 0.0
            node.m = 0.0

        self._restore_input_state()

    def _invalidate_solution(self):
        """Invalidate solved state while preserving a valid assembly."""
        for attribute in (
            "_K_reduced",
            "_rhs_reduced",
            "_free_dofs",
            "_prescribed_dofs",
            "_nodal_forces",
            "_reactions",
        ):
            if hasattr(self, attribute):
                delattr(self, attribute)

        if self._is_assembled:
            self._reset_input_vectors()
        else:
            for attribute in ("_u", "_f"):
                if hasattr(self, attribute):
                    delattr(self, attribute)

        self._reset_node_state()

    def _invalidate_assembly(self):
        """Invalidate global assembly and every dependent solved result."""
        self._is_assembled = False
        if hasattr(self, "_K"):
            del self._K
        self._invalidate_solution()

    
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

    def simple_report(self, report_type="print", fname="nusa_rpt.txt"):
        """Generate a compact text report for a solved finite-element model."""
        if not hasattr(self, "_nodal_forces"):
            raise RuntimeError("simple_report() is available only after solve()")

        valid_report_types = {"print", "string", "write"}
        if report_type not in valid_report_types:
            raise ValueError(
                f"Unknown report_type {report_type!r}; "
                f"expected one of {sorted(valid_report_types)}"
            )

        options = {
            "headers": "firstrow",
            "tablefmt": "rst",
            "numalign": "right",
        }
        sections = [
            "==========================",
            "    NuSA Simple Report",
            "==========================",
            "",
            f"Model: {self.name}",
            f"Number of nodes: {self.n_nodes}",
            f"Number of elements: {self.n_elements}",
            "",
            "RESULTS",
            "",
            "NODAL DISPLACEMENTS",
            self._get_ndisplacements(options),
            "",
            "APPLIED LOADS",
            self._get_applied_loads(options),
            "",
            "NODAL FORCES (K @ U)",
            self._get_nforces(options),
            "",
            "REACTIONS",
            self._get_reactions(options),
            "",
            "ELEMENT RESULTS",
            self._get_element_results(options),
            "",
            "FINITE ELEMENT MODEL INFO",
            "",
            "NODES",
            self._get_nodes_info(options),
            "",
            "ELEMENTS",
            self._get_elements_info(options),
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

    def _get_ndisplacements(self, options):
        """Generate a table of solved nodal displacement components."""
        from tabulate import tabulate

        dof_names = getattr(self, "displacement_dofs", ("ux", "uy"))
        headers = ["Node"] + [name.upper() for name in dof_names]
        rows = [headers]
        for node in self.nodes:
            rows.append([node.label] + [getattr(node, name) for name in dof_names])
        return tabulate(rows, **options)

    def _get_force_table(self, options, getter):
        """Generate a named-component nodal force table."""
        from tabulate import tabulate

        headers = ["Node"] + [name.upper() for name in self.force_dofs]
        rows = [headers]
        for node in self.nodes:
            values = getter(node)
            rows.append([node.label] + [values[name] for name in self.force_dofs])
        return tabulate(rows, **options)

    def _get_applied_loads(self, options):
        """Generate a table of explicitly applied nodal loads."""
        return self._get_force_table(options, self.applied_load)

    def _get_nforces(self, options):
        """Generate a table of solved generalized nodal forces (K @ u)."""
        if hasattr(self, "force_dofs") and hasattr(self, "_nodal_forces"):
            return self._get_force_table(options, self.nodal_force)

        from tabulate import tabulate

        rows = [["Node", "FX", "FY"]]
        for node in self.nodes:
            rows.append([node.label, node.fx, node.fy])
        return tabulate(rows, **options)

    def _get_reactions(self, options):
        """Generate a table of support reactions."""
        return self._get_force_table(options, self.reaction)

    def _get_element_results(self, options):
        """Generate the model-specific element-results table."""
        raise NotImplementedError(
            f"{self.__class__.__name__} must implement _get_element_results()"
        )

    def _get_nodes_info(self, options):
        """Generate a table of node coordinates."""
        from tabulate import tabulate

        rows = [["Node", "X", "Y"]]
        for node in self.nodes:
            rows.append([node.label, node.x, node.y])
        return tabulate(rows, **options)

    def _get_elements_info(self, options):
        """Generate a table of element connectivity."""
        from tabulate import tabulate

        max_nodes = max((len(element.nodes) for element in self.elements), default=0)
        headers = ["Element"] + [f"N{k + 1}" for k in range(max_nodes)]
        rows = [headers]
        for element in self.elements:
            labels = [node.label for node in element.nodes]
            rows.append(
                [element.label]
                + labels
                + [""] * (max_nodes - len(labels))
            )
        return tabulate(rows, **options)



#~ =========================== ELEMENT ===========================

class Element:
    """
    Superclass for all Elements
    """
    def __init__(self,etype):
        self.etype = etype # element type
        self.label = None
        self._fx = 0.0
        self._fy = 0.0
        self._sx = 0.0
        self._sy = 0.0
        self._sxy = 0.0
        
    @property
    def fx(self):
        return self._fx
        
    @fx.setter
    def fx(self,val):
        self._fx = val
        
    @property
    def fy(self):
        return self._fy
        
    @fy.setter
    def fy(self,val):
        self._fy = val
        
    def __str__(self):
        _str = str(self.__class__)
        return _str


#~ =========================== NODE ===========================

class Node:
    """
    Class for node object.
    """
    def __init__(self,coordinates):
        """
        Initialize a node with given coordinates.

        Parameters
        ----------
        coordinates : tuple
            A tuple containing the (x, y) coordinates of the node.
        """
        try:
            coordinates = np.asarray(coordinates, dtype=float)
        except (TypeError, ValueError) as exc:
            raise ValueError("Node coordinates must contain two finite numbers") from exc
        if coordinates.shape != (2,):
            raise ValueError("Node coordinates must contain exactly two values")
        if not np.isfinite(coordinates).all():
            raise ValueError("Node coordinates must be finite")
        self.coordinates = coordinates.copy()
        self._label = None

        # DOF
        self._ux = np.nan
        self._uy = np.nan
        self._ur = np.nan
        # Nodal forces
        self._fx = 0.0
        self._fy = 0.0
        self._m = 0.0
        # Nodal stresses
        self._sx = 0.0
        self._sy = 0.0
        self._sxy = 0.0
        self._seqv = 0.0 
        # strain
        self._ex = 0.0
        self._ey = 0.0
        self._exy = 0.0
        # Elements ¿what?
        self._elements = []
        
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
    def label(self,val):
        self._label = val

    def _add_element(self, element):
        """Register an attached element for internal nodal post-processing."""
        self._elements.append(element)
        
    @property
    def ux(self):
        """
        Return the x-displacement of the node.
        """
        return self._ux
    
    @ux.setter
    def ux(self,val):
        self._ux = val
    
    @property
    def uy(self):
        """
        Return the y-displacement of the node.
        """
        return self._uy
    
    @uy.setter
    def uy(self,val):
        self._uy = val
    
    @property
    def ur(self):
        """
        Return the rotational displacement of the node.
        """
        return self._ur
    
    @ur.setter
    def ur(self,val):
        self._ur = val
        
    @property
    def fx(self):
        """Solved generalized nodal x-force (``K @ u``), not a reaction."""
        return self._fx
    
    @fx.setter
    def fx(self,val):
        self._fx = val
    
    @property
    def fy(self):
        """Solved generalized nodal y-force (``K @ u``), not a reaction."""
        return self._fy
    
    @fy.setter
    def fy(self,val):
        self._fy = val
        
    @property
    def m(self):
        """Solved generalized nodal moment (``K @ u``), not a reaction."""
        return self._m
    
    @m.setter
    def m(self,val):
        self._m = val
        
    @property
    def sx(self):
        elements = self._elements
        if elements == []:
            self._sx = 0.0
        else:
            self._sx = sum([el.sx for el in elements])/len(elements)
        return self._sx
    
    @sx.setter
    def sx(self,val):
        self._sx = val
        
    @property
    def sy(self):
        elements = self._elements
        if elements == []:
            self._sy = 0
        else:
            self._sy = sum([el.sy for el in elements])/len(elements)
        return self._sy
    
    @sy.setter
    def sy(self,val):
        self._sy = val
        
    @property
    def sxy(self):
        elements = self._elements
        if elements == []:
            self._sxy = 0
        else:
            self._sxy = sum([el.sxy for el in elements])/len(elements)
        return self._sxy
    
    @sxy.setter
    def sxy(self,val):
        self._sxy = val
        
    @property
    def seqv(self):
        sxx, syy, sxy = self.sx, self.sy, self.sxy
        seqv = np.sqrt(sxx**2 - sxx*syy + syy**2 + 3*sxy**2)
        return seqv

    @property
    def ex(self):
        elements = self._elements
        if elements == []:
            self._ex = 0
        else:
            self._ex = sum([el.ex for el in elements])/len(elements)
        return self._ex
    
    @ex.setter
    def ex(self,val):
        self._ex = val

    @property
    def ey(self):
        elements = self._elements
        if elements == []:
            self._ey = 0
        else:
            self._ey = sum([el.ey for el in elements])/len(elements)
        return self._ey
    
    @ey.setter
    def ey(self,val):
        self._ey = val

    @property
    def exy(self):
        elements = self._elements
        if elements == []:
            self._exy = 0
        else:
            self._exy = sum([el.exy for el in elements])/len(elements)
        return self._exy
    
    @exy.setter
    def exy(self,val):
        self._exy = val

    def __str__(self):
        _str = self.__class__
        _str = "%s\nU:(%g,%g)\n"%(_str,self.ux, self.uy)
        _str = "%sF:(%g,%g)"%(_str,self.fx,self.fy)
        return _str
    
    def __repr__(self):
        return f"<Node {self.label}: ({self.x},{self.y})>"
        

if __name__=='__main__':
    pass
