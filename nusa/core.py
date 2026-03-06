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
        self._nodes = {} # Dictionary for nodes {number: NodeObject}
        self._elements = {} # Dictionary for elements {number: ElementObject}
        
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
        if node.label is None:
            node.label = self.n_nodes

        self._nodes[node.label] = node

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

        if element.etype != self.mtype:
            raise ValueError(
                f"Element type '{element.etype}' incompatible with model '{self.mtype}'"
            )

        if element.label is None:
            element.label = self.n_elements

        self._elements[element.label] = element

        for node in element.nodes:
            node.add_element(element)

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
    
    @property
    def nodes(self):
        """
        Return a list of node objects.

        Returns
        -------
        list
            List of Node instances.
        """
        return list(self._nodes.values())

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

    def simple_report(self,report_type="print",fname="nusa_rpt.txt"):
        """
        Placeholder for a future implementation of a simple report.

        Parameters
        ----------
        report_type : str, optional
            Type of report to generate ('print', 'file', etc.).
        fname : str, optional
            Output filename for file-based reports.
        """
        pass
        
    def _get_ndisplacements(self,options):
        """
        Generate a table of node displacements.

        Parameters
        ----------
        options : dict
            Tabulate formatting options.

        Returns
        -------
        str
            Tabulated string of displacements.
        """
        from tabulate import tabulate
        D = [["Node","UX","UY"]]
        for n in self.get_nodes():
            D.append([n.label+1,n.ux,n.uy])
        return tabulate(D, **options)
        
    def _get_nforces(self,options):
        """
        Generate a table of nodal forces.

        Parameters
        ----------
        options : dict
            Tabulate formatting options.

        Returns
        -------
        str
            Tabulated string of nodal forces.
        """
        from tabulate import tabulate
        F = [["Node","FX","FY"]]
        for n in self.get_nodes():
            F.append([n.label+1,n.fx,n.fy])
        return tabulate(F, **options)
        
    def _get_eforces(self,options):
        """
        Generate a table of element internal forces.

        Parameters
        ----------
        options : dict
            Tabulate formatting options.

        Returns
        -------
        str
            Tabulated string of element forces.
        """
        from tabulate import tabulate
        F = [["Element","F"]]
        for elm in self.get_elements():
            F.append([elm.label+1, elm.f])
        return tabulate(F, **options)
        
    def _get_estresses(self,options):
        """
        Generate a table of element stresses.

        Parameters
        ----------
        options : dict
            Tabulate formatting options.

        Returns
        -------
        str
            Tabulated string of element stresses.
        """
        from tabulate import tabulate
        S = [["Element","S"]]
        for elm in self.get_elements():
            S.append([elm.label+1, elm.s])
        return tabulate(S, **options)
    
    def _get_nodes_info(self,options):
        """
        Generate a table of node coordinates.

        Parameters
        ----------
        options : dict
            Tabulate formatting options.

        Returns
        -------
        str
            Tabulated string of node positions.
        """
        from tabulate import tabulate
        F = [["Node","X","Y"]]
        for n in self.get_nodes():
            F.append([n.label+1, n.x, n.y])
        return tabulate(F, **options)
    
    def _get_elements_info(self,options):
        """
        Generate a table of element connectivity.

        Parameters
        ----------
        options : dict
            Tabulate formatting options.

        Returns
        -------
        str
            Tabulated string of element-node relationships.
        """
        from tabulate import tabulate
        S = [["Element","NI","NJ"]]
        for elm in self.get_elements():
            ni, nj = elm.get_nodes()
            S.append([elm.label+1, ni.label+1, nj.label+1])
        return tabulate(S, **options)
            



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
        
    def set_label(self,label):
        """
        Set the label property
        
        *label* : int
            Label, must be an integer
        """
        self.label = label
        
    def set_element_forces(self,fx=0.0,fy=0.0):
        """
        Set element forces
        
        *fx* : float
            Force in x-dir
        *fy* : float
            Force in y-dir
        
        Normally this method is used by the `solve` method to 
        update computed element-forces.
        """
        self._fx = fx
        self._fy = fy
        
    def get_element_forces(self):
        """
        Returns a tuple with element forces:  (fx, fy)
        """
        return self._fx, self._fy
        
    def get_nodes(self):
        return self.nodes
        
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
        self.coordinates = np.asanyarray(coordinates, dtype=float)
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

    def add_element(self,element):
        self._elements.append(element)
        
    @property
    def ux(self):
        return self._ux
    
    @ux.setter
    def ux(self,val):
        self._ux = val
    
    @property
    def uy(self):
        return self._uy
    
    @uy.setter
    def uy(self,val):
        self._uy = val
    
    @property
    def ur(self):
        return self._ur
    
    @ur.setter
    def ur(self,val):
        if True:#type(val) in [int,float]:
            self._ur = val
        else:
            raise ValueError("Value must be float or int")
        
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
        
    @property
    def m(self):
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

    def get_label(self):
        return self._label
    
    def set_label(self,label):
        self._label = label
    
    def get_displacements(self):
        return self._ux,self._uy,self._ur
        
    def set_displacements(self,ux=np.nan, uy=np.nan, ur=np.nan):
        self._ux = ux
        self._uy = uy
        self._ur = ur
    
    def get_forces(self):
        return (self._fx,self._fy)
    
    def set_forces(self,fx=np.nan,fy=np.nan):
        self._fx = fx
        self._fy = fy
        
    def __str__(self):
        _str = self.__class__
        _str = "%s\nU:(%g,%g)\n"%(_str,self.ux, self.uy)
        _str = "%sF:(%g,%g)"%(_str,self.fx,self.fy)
        return _str
    
    def __repr__(self):
        return f"<Node {self.label}: ({self.x},{self.y})>"
        

if __name__=='__main__':
    pass
