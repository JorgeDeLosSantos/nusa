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
from core import Element, Node, Model


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
    
    
    """
    def __init__(self,nodes,E,A,L,theta):
        Element.__init__(self,etype="truss")
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



#~ *********************************************************************
#~ ****************************  TrussModel ****************************
#~ *********************************************************************
class TrussModel(Model):
    """
    Model for finite element analysis
    """
    def __init__(self,name="Truss Model 01"):
        Model.__init__(self,name=name,mtype="truss")
        self.F = {} # Forces
        self.U = {} # Displacements
        self.dof = 2 # 2 DOF for truss element
        self.IS_KG_BUILDED = False
        
    def buildGlobalMatrix(self):
        msz = (self.dof)*self.getNumberOfNodes()
        self.KG = np.zeros((msz,msz))
        for element in self.elements.values():
            ku = element.getElementStiffness()
            n1,n2 = element.getNodes()
            self.KG[2*n1.label, 2*n1.label] += ku[0,0]
            self.KG[2*n1.label, 2*n1.label+1] += ku[0,1]
            self.KG[2*n1.label, 2*n2.label] += ku[0,2]
            self.KG[2*n1.label, 2*n2.label+1] += ku[0,3]
            
            self.KG[2*n1.label+1, 2*n1.label] += ku[1,0]
            self.KG[2*n1.label+1, 2*n1.label+1] += ku[1,1]
            self.KG[2*n1.label+1, 2*n2.label] += ku[1,2]
            self.KG[2*n1.label+1, 2*n2.label+1] += ku[1,3]
            
            self.KG[2*n2.label, 2*n1.label] += ku[2,0]
            self.KG[2*n2.label, 2*n1.label+1] += ku[2,1]
            self.KG[2*n2.label, 2*n2.label] += ku[2,2]
            self.KG[2*n2.label, 2*n2.label+1] += ku[2,3]
            
            self.KG[2*n2.label+1, 2*n1.label] += ku[3,0]
            self.KG[2*n2.label+1, 2*n1.label+1] += ku[3,1]
            self.KG[2*n2.label+1, 2*n2.label] += ku[3,2]
            self.KG[2*n2.label+1, 2*n2.label+1] += ku[3,3]
            
        self.buildForcesVector()
        self.buildDisplacementsVector()
        self.IS_KG_BUILDED = True
        
    def buildForcesVector(self):
        for node in self.nodes.values():
            self.F[node.label] = {"fx":0, "fy":0}
        
    def buildDisplacementsVector(self):
        for node in self.nodes.values():
            self.U[node.label] = {"ux":np.nan, "uy":np.nan}
    
    def addForce(self,node,force):
        node.setForces(fx=force[0], fy=force[1])
        
    def addConstraint(self,node,**constraint):
        cs = constraint
        if cs.has_key('ux') and cs.has_key("uy"):
            ux = cs.get('ux')
            uy = cs.get('uy')
            node.setDisplacements(ux=ux, uy=uy)
        elif cs.has_key('ux'):
            ux = cs.get('ux')
            node.setDisplacements(ux=ux)
        elif cs.has_key('uy'):
            uy = cs.get('uy')
            node.setDisplacements(uy=uy)
        
    def solve(self):
        self.buildDisplacementsVector()
        self.buildForcesVector()
        knw = [pos for pos,value in enumerate(self.U) if not np.isnan(value)]
        unknw = [pos for pos,value in enumerate(self.U) if np.isnan(value)]
        self.K2S = np.delete(np.delete(self.KG,knw,0),knw,1)
        self.F2S = np.delete(self.F,knw,0)
        # For displacements
        self.solved_u = la.solve(self.K2S,self.F2S)
        for k,idx in enumerate(unknw):
            nd = idx/2
            self.U[idx] = self.solved_u[k]
            if (-1)**idx > 0:
                self.nodes[nd].ux = self.solved_u[k]
            else:
                self.nodes[nd].uy = self.solved_u[k]
        #~ # For nodal forces/reactions
        self.NF = self.F.copy()
        nf_calc = np.dot(self.KG, self.U)
        for nd in self.nodes.keys():
            self.NF[2*nd] = nf_calc[2*nd]
            self.NF[2*nd+1] = nf_calc[2*nd+1]
        print self.NF

        

if __name__=='__main__':
    pass
