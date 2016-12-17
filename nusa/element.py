# -*- coding: utf-8 -*-
# ***********************************
#  Author: Pedro Jorge De Los Santos     
#  E-mail: delossantosmfq@gmail.com 
#  Blog: numython.github.io
#  License: MIT License
# ***********************************
from __future__ import division
import numpy as np
import numpy.linalg as la
from core import Element
import templates as tmp

class Spring(Element):
    """
    Spring element for finite element analysis
    
    *nodes* : tuple
        Connectivity for element given as tuple of 
        :class:`~core.base.Node` objects
    
    *ke* : float
        Spring stiffness
    
    Example ::
    
        n1 = Node((0,0))
        n2 = Node((0,0))
        e1 = Spring((n1,n2), 1000)
    
    """
    def __init__(self,nodes,ke):
        Element.__init__(self, etype="spring")
        self.nodes = nodes
        self.k = k
    
    @property
    def fx(self):
        ke = self.getElementStiffness() # Element stiffness
        n1, n2 = self.getNodes()
        un = np.array([[n1.ux],[n2.ux]]) # Nodal displacements
        return np.dot(ke, un) # Return  {fxe} = [Ke]{uxe}
        
    @fx.setter
    def fx(self,val):
        self._fx = val
        
    def getElementStiffness(self):
        r"""
        Get stiffness matrix for this element.
        
        The stiffness matrix for a spring element is defined by:
        
        .. math::
        
            [k]_e = \begin{bmatrix}
            k  & -k \\
            -k &  k \\
            \end{bmatrix}
        
        where *k* is the spring stiffness.
        
        Return a numpy array.
        """
        self._KE = np.array([[self.k,-self.k],[-self.k,self.k]])
        return self._KE
    
    def getNodes(self):
        """
        Returns a tuple of nodes
        """
        return self.nodes



class Bar(Element):
    """
    Bar element for finite element analysis
    
    *nodes* : :class:`~core.base.Node`
        Connectivity for element
    
    *E* : float
        Young's modulus
        
    *A* : float
        Area of element
        
    *L* : float
        Length of element
    
    
    Example::
    
        n1 = Node((0,0))
        n2 = Node((1,0))
        e1 = Bar((n1,n2),200e9,0.02,1)
    
    """
    def __init__(self,nodes,E,A,L):
        Element.__init__(self,etype="bar")
        self.nodes = nodes
        self.E = E # Elastic modulus
        self.A = A # Cross-section
        self.L = L # Length
        
    @property
    def fx(self):
        """
        Compute force in x-dir (axial-dir)
        """
        ke = self.getElementStiffness() # Element stiffness
        n1, n2 = self.getNodes()
        un = np.array([[n1.ux],[n2.ux]]) # Nodal displacements
        return np.dot(ke, un) # Return  {fxe} = [Ke]{uxe}
        
    @fx.setter
    def fx(self,val):
        self._fx = val
        
    @property
    def sx(self):
        r"""
        Compute normal stress in x-dir
        
        Given by:
        
        .. math::
            
            S_x = K \left( \frac{u}{A} \right)
        
        where
        
        * K - Element stiffness matrix
        * u - Nodal displacements
        * A - Cross-section of element
        """
        ke = self.getElementStiffness() # Element stiffness
        na, nb = self.getNodes()
        u = np.array([na.ux, nb.ux]) # Nodes displacements
        sx = np.dot(ke, u/self.A) # matrix multiplication
        return sx
        
    @sx.setter
    def sx(self,val):
        self._sx = val
        
        
    def getElementStiffness(self):
        r"""
        Get stiffness matrix for this element
        
        The stiffness matrix for bar element is given by:
        
        .. math::
        
            [k]_e = \frac{AE}{L} \begin{bmatrix} 1 & -1 \\ -1 & 1 \end{bmatrix}
        
        where
        
        * A - Cross-section of element
        * E - Young's Modulus
        * L - Length of element
        """
        self._KE = (self.A*self.E/self.L)*np.array([[1,-1],[-1,1]])
        return self._KE
        
    def getNodes(self):
        """
        Returns a tuple of nodes
        """
        return self.nodes




class Beam(Element):
    """
    Beam element for finite element analysis
    
    *nodes* : :class:`~core.base.Node`
        Connectivity for element
    
    *E* : float
        Young's modulus
        
    *I* : float
        Moment of inertia
        
    *L* : float
        Length of element
    
    Example::
    
        E = 30e6
        I = 500.0
        L = 100
        n1 = Node((0,0))
        n2 = Node((1,0))
        e1 = Bar((n1,n2),E,I,L)
    
    """
    def __init__(self,nodes,E,I,L):
        Element.__init__(self,etype="beam")
        self.nodes = nodes
        self.E = E
        self.I = I
        self.L = L
        
    def getElementStiffness(self):
        """
        Get stiffness matrix for this element
        
        """
        multiplier = (self.I*self.E/self.L**3)
        a = 6*self.L
        b = 4*self.L**2
        c = 2*self.L**2
        self._K = multiplier*np.array([[ 12, a, -12, a],
                                       [  a, b,  -a, c],
                                       [-12,-a,  12,-a],
                                       [  a, c,  -a, b]])
        return self._K

    def _compute_element_forces(self):
        """
        Just that 
        
        Set fy and m properties.
        """
        ke = self.getElementStiffness() # Element stiffness
        n1, n2 = self.getNodes()
        un = np.array([[n1.uy, n1.ur, n2.uy, n2.ur]]).transpose() # Nodal displacements
        EF = np.dot(ke, un) # Return  {fxe} = [Ke]{uxe}
        self.fy = EF[::2] # Set fy
        self.m = EF[1::2] # Set m
        
    @property
    def fy(self):
        """
        Compute y-force 
        """
        self._compute_element_forces()
        return self._fy
        
    @fy.setter
    def fy(self,val):
        self._fy = val
        
    @property
    def m(self):
        """
        Compute moment 
        """
        self._compute_element_forces()
        return self._m
        
    @m.setter
    def m(self,val):
        self._m = val
        
    def getNodes(self):
        return self.nodes



class Truss(Element):
    """
    Truss element for finite element analysis
    
    *nodes* : :class:`~core.base.Node`
        Connectivity for element
    
    *E* : float
        Young modulus
        
    *A* : float
        Area of element
        
    *L* : float
        Length of element
    
    
    Example::
    
        n1 = Node((0,0),1)
        n2 = Node((1,0),1)
        e1 = Bar((n1,n2),200e9,0.02,1)
    
    """
    def __init__(self,nodes,E,A,L,theta,label=""):
        Element.__init__(self,etype="truss",label=label)
        self.nodes = nodes
        self.E = E
        self.A = A
        self.L = L
        self.theta = theta
        
    def getElementStiffness(self):
        """
        Get stiffness matrix for this element
        """
        multiplier = (self.A*self.E/self.L)
        C = np.cos(self.theta)
        S = np.sin(self.theta)
        CS = C*S
        self._K = multiplier*np.array([[C**2 , CS   , -C**2, -CS  ],
                                       [CS   , S**2 , -CS  , -S**2],
                                       [-C**2, -CS  , C**2 , CS   ],
                                       [-CS  , -S**2,  CS  , S**2 ]])
        return self._K
        
    def getNodes(self):
        return self.nodes





if __name__=='__main__':
    pass
