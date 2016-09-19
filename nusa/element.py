# -*- coding: utf-8 -*-
# ***********************************
#  Author: Pedro Jorge De Los Santos     
#  E-mail: delossantosmfq@gmail.com 
#  Blog: numython.blogspot.mx
#  License: MIT License
# ***********************************
import numpy as np
import numpy.linalg as la
from core import Element
import templates as tmp

class Spring(Element):
    """
    Bar element for finite element analysis
    
    *nodes* : :class:`~core.base.Node`
        Conectivity for element
    
    *E* : float
        Young modulus
        
    *A* : float
        Area of element
        
    *L* : float
        Length of element
    
    
    Example::
    
        n1 = Node((0,0),1)
        n2 = Node((1,0),1)
        e1 = Spring((n1,n2),200e9,0.02,1)
    
    """
    def __init__(self,nodes,ke,label=""):
        Element.__init__(self,etype="spring",label=label)
        self.nodes = nodes
        self.ke = ke
    
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
        Get stiffness matrix for this element
        
        The stiffness matrix for a spring element is defined by:
        
        .. math::
        
            [K]_e = \left[\begin{matrix}
            k  & -k \\
            -k &  k \\
            \end{matrix} \right]
        
        Return a numpy array.
        """
        self._KE = np.array([[self.ke,-self.ke],[-self.ke,self.ke]])
        return self._KE
    
    def getNodes(self):
        return self.nodes



class Bar(Element):
    """
    Bar element for finite element analysis
    
    *nodes* : :class:`~core.base.Node`
        Conectivity for element
    
    *E* : float
        Young modulus
        
    *A* : float
        Area of element
        
    *L* : float
        Length of element
    
    
    Example::
    
        n1 = Node((0,0))
        n2 = Node((1,0))
        e1 = Bar((n1,n2),200e9,0.02,1)
    
    """
    def __init__(self,nodes,E,A,L,label=""):
        Element.__init__(self,etype="bar",label=label)
        self.nodes = nodes
        self.E = E # Elastic modulus
        self.A = A # Cross-section
        self.L = L # Length
        
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
        """
        Get stiffness matrix for this element
        """
        self._KE = (self.A*self.E/self.L)*np.array([[1,-1],[-1,1]])
        return self._KE
        
    def getNodes(self):
        return self.nodes




class Beam(Element):
    """
    Beam element for finite element analysis
    
    *nodes* : :class:`~core.base.Node`
        Conectivity for element
    
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
    def __init__(self,nodes,E,I,L,label=""):
        Element.__init__(self,etype="beam",label=label)
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
        
    def getNodes(self):
        return self.nodes


class Truss(Element):
    """
    Truss element for finite element analysis
    
    *nodes* : :class:`~core.base.Node`
        Conectivity for element
    
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
