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
    def __init__(self,nodes,E,A):
        Element.__init__(self,etype="truss")
        self.nodes = nodes
        self.E = E
        self.A = A
        #~ self.theta = theta
        
    @property
    def L(self):
        ni,nj = self.getNodes()
        x0,x1,y0,y1 = ni.x, nj.x, ni.y, nj.y
        _l = np.sqrt( (x1-x0)**2 + (y1-y0)**2 )
        return _l
    
    @property
    def theta(self):
        ni,nj = self.getNodes()
        x0,x1,y0,y1 = ni.x, nj.x, ni.y, nj.y
        if x0==x1:
            theta = 90*(np.pi/180)
        else:
            theta = np.arctan((y1-y0)/(x1-x0))
        return theta
    
    @property
    def f(self):
        return self._compute_force()
        
    def _compute_force(self):
        theta = self.theta
        E, A, L = self.E, self.A, self.L
        C = np.cos(theta)
        S = np.sin(theta)
        ni, nj = self.getNodes()
        u = np.array([ni.ux, ni.uy, nj.ux, nj.uy]).T
        F = (E*A/L)*np.dot(np.array([-C, -S, C, S]), u)
        return F
        
        
        
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
        if not(self.IS_KG_BUILDED): self.buildGlobalMatrix()
        self.F[node.label]["fx"] = force[0]
        self.F[node.label]["fy"] = force[1]
        node.fx = force[0]
        node.fy = force[1]
        
    def addConstraint(self,node,**constraint):
        if not(self.IS_KG_BUILDED): self.buildGlobalMatrix()
        cs = constraint
        if cs.has_key('ux') and cs.has_key("uy"): # 
            ux = cs.get('ux')
            uy = cs.get('uy')
            node.setDisplacements(ux=ux, uy=uy)
            self.U[node.label]["ux"] = ux
            self.U[node.label]["uy"] = uy
        elif cs.has_key('uy'):
            uy = cs.get('uy')
            node.setDisplacements(uy=uy)
            self.U[node.label]["uy"] = uy
        
    def solve(self):
        # Solve LS
        self.VU = [node[key] for node in self.U.values() for key in ("ux","uy")]
        self.VF = [node[key] for node in self.F.values() for key in ("fx","fy")]
        knw = [pos for pos,value in enumerate(self.VU) if not value is np.nan]
        unknw = [pos for pos,value in enumerate(self.VU) if value is np.nan]
        self.K2S = np.delete(np.delete(self.KG,knw,0),knw,1)
        self.F2S = np.delete(self.VF,knw,0)
        
        # For displacements
        self.solved_u = la.solve(self.K2S,self.F2S)
        for k,ic in enumerate(unknw):
            nd, var = self.index2key(ic)
            self.U[nd][var] = self.solved_u[k]
            
        # Updating nodes displacements
        for nd in self.nodes.values():
            if np.isnan(nd.ux):
                nd.ux = self.U[nd.label]["ux"]
            if np.isnan(nd.uy):
                nd.uy = self.U[nd.label]["uy"]
                    
        # For nodal forces/reactions
        self.NF = self.F.copy()
        self.VU = [node[key] for node in self.U.values() for key in ("ux","uy")]
        nf_calc = np.dot(self.KG, self.VU)
        for k in range(2*self.getNumberOfNodes()):
            nd, var = self.index2key(k, ("fx","fy"))
            self.NF[nd][var] = nf_calc[k]
            cnlab = np.floor(k/float(self.dof))
            if var=="fx": 
                self.nodes[cnlab].fx = nf_calc[k]
            elif var=="fy":
                self.nodes[cnlab].fy = nf_calc[k]
                
    def index2key(self,idx,opts=("ux","uy")):
        """
        Index to key, where key can be ux or uy
        """
        node = idx//2
        var = opts[0] if ((-1)**idx)==1 else opts[1]
        return node,var
        
    def plot_model(self):
        """
        Plot the mesh model, including bcs
        """
        import matplotlib.pyplot as plt
        from matplotlib.patches import Polygon
        from matplotlib.collections import PatchCollection
        
        fig = plt.figure()
        ax = fig.add_subplot(111)
        
        for elm in self.getElements():
            ni, nj = elm.getNodes()
            ax.plot([ni.x,nj.x],[ni.y,nj.y],"k-")
            for nd in (ni,nj):
                if nd.fx != 0: self._draw_xforce(ax,nd.x,nd.y)
                if nd.fy != 0: self._draw_yforce(ax,nd.x,nd.y)
                if nd.ux == 0 and nd.uy == 0: self._draw_xyconstraint(ax,nd.x,nd.y)
        
        x0,x1,y0,y1 = self.rect_region()
        ax.set_xlim(x0,x1)
        ax.set_ylim(y0,y1)

    def _draw_xforce(self,axes,x,y):
        """
        Draw horizontal arrow -> Force in x-dir
        """
        dx, dy = self._calculate_arrow_size(), 0
        HW = dx/5.0
        HL = dx/3.0
        arrow_props = dict(head_width=HW, head_length=HL, fc='r', ec='r')
        axes.arrow(x, y, dx, dy, **arrow_props)
        
    def _draw_yforce(self,axes,x,y):
        """
        Draw vertical arrow -> Force in y-dir
        """
        dx,dy = 0, self._calculate_arrow_size()
        HW = dy/5.0
        HL = dy/3.0
        arrow_props = dict(head_width=HW, head_length=HL, fc='r', ec='r')
        axes.arrow(x, y, dx, dy, **arrow_props)
        
    def _draw_xyconstraint(self,axes,x,y):
        axes.plot(x, y, "gv", markersize=10, alpha=0.6)
        axes.plot(x, y, "g<", markersize=10, alpha=0.6)
        
    def _calculate_arrow_size(self):
        x0,x1,y0,y1 = self.rect_region(factor=10)
        sf = 8e-2
        kfx = sf*(x1-x0)
        kfy = sf*(y1-y0)
        return np.mean([kfx,kfy])
        

    def show(self):
        import matplotlib.pyplot as plt
        plt.show()
        
    def rect_region(self,factor=7.0):
        nx,ny = [],[]
        for n in self.getNodes():
            nx.append(n.x)
            ny.append(n.y)
        xmn,xmx,ymn,ymx = min(nx),max(nx),min(ny),max(ny)
        kx = (xmx-xmn)/factor
        ky = (ymx-ymn)/factor
        return xmn-kx, xmx+kx, ymn-ky, ymx+ky

        

if __name__=='__main__':
    pass
