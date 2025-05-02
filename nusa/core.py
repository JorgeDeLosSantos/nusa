# ***********************************
#  Author: Pedro Jorge De Los Santos    
#  E-mail: delossantosmfq@gmail.com 
#  Blog: numython.github.io
#  License: MIT License
# ***********************************
import numpy as np

#~ ===========================  MODEL  ===========================
class Model(object):
    """
    Superclass for all Finite Element Analysis (FEA) models.

    This class serves as a base container to manage nodes and elements,
    allowing derived models to build and manipulate FEA structures.
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
        self.nodes = {} # Dictionary for nodes {number: NodeObject}
        self.elements = {} # Dictionary for elements {number: ElementObject}
        
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
        current_label = self.get_number_of_nodes()
        if node.label is "":
            node.set_label(current_label)
        self.nodes[node.label] = node
        
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
        if self.mtype != element.etype:
            raise ValueError("Element type must be "+self.mtype)
        current_label = self.get_number_of_elements()
        if element.label is "":
            element.set_label(current_label)
        self.elements[element.label] = element
        # Assign this element to "xxxx" 
        for node in element.get_nodes():
            node._elements.append(element)

    def get_number_of_nodes(self):
        """
        Return the number of nodes in the model.

        Returns
        -------
        int
            Total number of nodes.
        """
        return len(self.nodes)
        
    def get_number_of_elements(self):
        """
        Return the number of elements in the model.

        Returns
        -------
        int
            Total number of elements.
        """
        return len(self.elements)
        
    def get_nodes(self):
        """
        Return a list of node objects.

        Returns
        -------
        list
            List of Node instances.
        """
        return self.nodes.values()
        
    def get_elements(self):
        """
        Return a list of element objects.

        Returns
        -------
        list
            List of Element instances.
        """
        return self.elements.values()
    
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
            f"Nodes: {self.get_number_of_nodes()}\n"
            f"Elements: {self.get_number_of_elements()}"
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
            f"Nodes: {self.get_number_of_nodes()}\n"
            f"Elements: {self.get_number_of_elements()}"
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

class Element(object):
    """
    Superclass for all Elements
    """
    def __init__(self,etype):
        self.etype = etype # element type
        self.label = "" # label (reassignment -> Model.addElement)
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

class Node(object):
    """
    Class for node object.
    
    *coordinates* : `tuple`, `list`
        Coordinates of node
    
    *label* : int
        Label of node
        
    ::
    
        n1 = Node((0,0))
        n2 = Node((0,0))
    
    """
    def __init__(self,coordinates):
        self.coordinates = coordinates
        self.x = coordinates[0] # usable prop
        self.y = coordinates[1] # usable prop
        self._label = ""
        self._ux = np.nan
        self._uy = np.nan
        self._ur = np.nan
        self._fx = 0.0
        self._fy = 0.0
        self._m = 0.0
        # Nodal stresses
        self._sx = 0.0
        self._sy = 0.0
        self._sxy = 0.0
        self._seqv = 0.0 
        # Elements ¿what?
        self._elements = []
        
        
    @property
    def label(self):
        return self._label
        
    @label.setter
    def label(self,val):
        """
        Experimental setter for adjust range of labels: TO DO
        """
        self._label = val
        
    @property
    def ux(self):
        return self._ux
    
    @ux.setter
    def ux(self,val):
        if True:#type(val) in [int,float]:
            self._ux = val
        else:
            raise ValueError("Value must be float or int")
    
    @property
    def uy(self):
        return self._uy
    
    @uy.setter
    def uy(self,val):
        if True:#type(val) in [int,float]:
            self._uy = val
        else:
            raise ValueError("Value must be float or int")
    
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
        

if __name__=='__main__':
    pass
