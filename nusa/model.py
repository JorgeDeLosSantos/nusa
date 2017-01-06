# -*- coding: utf-8 -*-
# ***********************************
#  Author: Pedro Jorge De Los Santos     
#  E-mail: delossantosmfq@gmail.com 
#  Blog: numython.github.io
#  License: MIT License
# ***********************************
from __future__ import division
import re
import numpy as np
import numpy.linalg as la
import templates as tmp
from core import Model

#~ *********************************************************************
#~ ****************************  SpringModel ***************************
#~ *********************************************************************
class SpringModel(Model):
    """
    Spring Model for finite element analysis
    """
    def __init__(self,name="Spring Model 01"):
        Model.__init__(self,name=name,mtype="spring")
        self.F = {} # Forces
        self.U = {} # Displacements
        self.dof = 1 # 1 DOF per Node
        self.IS_KG_BUILDED = False

    def _buildGlobalMatrix(self):
        msz = (self.dof)*self.getNumberOfNodes() # Matrix size
        self.KG = np.zeros((msz,msz))
        for element in self.elements.values():
            ku = element.getElementStiffness()
            n1,n2 = element.getNodes()
            self.KG[n1.label, n1.label] += ku[0,0]
            self.KG[n1.label, n2.label] += ku[0,1]
            self.KG[n2.label, n1.label] += ku[1,0]
            self.KG[n2.label, n2.label] += ku[1,1]
        
        self.buildForcesVector()
        self.buildDisplacementsVector()
        self.IS_KG_BUILDED = True
        
    def buildGlobalMatrix(self):
        msz = (self.dof)*self.getNumberOfNodes() # Matrix size
        self.KG = np.zeros((msz,msz))
        for element in self.getElements():
            self.KG += element.get_global_stiffness(msz)
            
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
        
    def addConstraint(self,node,**constraint):
        """
        Only displacement in x-dir 
        """
        if not(self.IS_KG_BUILDED): self.buildGlobalMatrix()
        if constraint.has_key("ux"):
            ux = constraint.get("ux")
            node.setDisplacements(ux=ux)
            self.U[node.label]["ux"] = ux
        
    def solve(self):
        # known and unknown values
        self.VU = [node[key] for node in self.U.values() for key in ("ux",)]
        self.VF = [node[key] for node in self.F.values() for key in ("fx",)]
        knw = [pos for pos,value in enumerate(self.VU) if not value is np.nan]
        unknw = [pos for pos,value in enumerate(self.VU) if value is np.nan]
        # Matrices to solve
        self.K2S = np.delete(np.delete(self.KG,knw,0),knw,1)
        self.F2S = np.delete(self.VF,knw,0)
        # For displacements
        self.solved_u = la.solve(self.K2S,self.F2S)
        # Updating U (displacements vector)
        for k,ic in enumerate(unknw):
            nd, var = self.index2key(ic)
            self.U[nd][var] = self.solved_u[k]
            self.nodes[ic].ux = self.solved_u[k]
        # For nodal forces/reactions
        self.NF = self.F.copy()
        self.VU = [node[key] for node in self.U.values() for key in ("ux",)]
        nf_calc = np.dot(self.KG, self.VU)
        for k,ic in enumerate(range(self.getNumberOfNodes())):
            nd, var = self.index2key(ic, ("fx",))
            self.NF[nd][var] = nf_calc[k]
            self.nodes[ic].fx = nf_calc[k]
            
    def index2key(self,idx,opts=("ux",)):
        node = idx
        var = opts[0]
        return node,var



#~ *********************************************************************
#~ ****************************  BarModel ******************************
#~ *********************************************************************
class BarModel(Model):
    """
    Bar model for finite element analysis
    """
    def __init__(self,name="Bar Model 01"):
        Model.__init__(self,name=name,mtype="bar")
        self.F = {} # Forces
        self.U = {} # Displacements
        self.dof = 1 # 1 DOF for bar element (per node)
        self.IS_KG_BUILDED = False
        
    def buildForcesVector(self):
        for node in self.nodes.values():
            self.F[node.label] = {"fx":0, "fy":0}
        
    def buildGlobalMatrix(self):
        msz = (self.dof)*self.getNumberOfNodes()
        self.KG = np.zeros((msz,msz))
        for element in self.elements.values():
            ku = element.getElementStiffness()
            n1,n2 = element.getNodes()
            self.KG[n1.label, n1.label] += ku[0,0]
            self.KG[n1.label, n2.label] += ku[0,1]
            self.KG[n2.label, n1.label] += ku[1,0]
            self.KG[n2.label, n2.label] += ku[1,1]
        self.buildForcesVector()
        self.buildDisplacementsVector()
        self.IS_KG_BUILDED = True
        
    def buildDisplacementsVector(self):
        for node in self.nodes.values():
            self.U[node.label] = {"ux":np.nan, "uy":np.nan}
        
    def addForce(self,node,force):
        if not(self.IS_KG_BUILDED): self.buildGlobalMatrix()
        self.F[node.label]["fx"] = force[0]
        
    def addConstraint(self,node,**constraint):
        if not(self.IS_KG_BUILDED): self.buildGlobalMatrix()
        if constraint.has_key('ux'):
            ux = constraint.get('ux')
            node.setDisplacements(ux=ux)
            self.U[node.label]["ux"] = ux
        
    def solve(self):
        # known and unknown values
        self.VU = [node[key] for node in self.U.values() for key in ("ux",)]
        self.VF = [node[key] for node in self.F.values() for key in ("fx",)]
        knw = [pos for pos,value in enumerate(self.VU) if not value is np.nan]
        unknw = [pos for pos,value in enumerate(self.VU) if value is np.nan]
        
        if len(unknw)==1:
            _k = unknw[0]
            _rowtmp = self.KG[_k,:]
            _ftmp = self.VF[_k]
            _fk = _ftmp - np.dot(np.delete(_rowtmp,_k), np.delete(self.VU,_k))
            _uk = _fk / self.KG[_k, _k]
            # Then 
            self.solved_u = np.array([_uk])
        else: # "Normal" case
            self.K2S = np.delete(np.delete(self.KG,knw,0),knw,1)
            self.F2S = np.delete(self.VF,knw,0)
            self.solved_u = la.solve(self.K2S,self.F2S)
            
        # For displacements
        # Updating U (displacements vector)
        for k,ic in enumerate(unknw):
            nd, var = self.index2key(ic)
            self.U[nd][var] = self.solved_u[k]
            self.nodes[ic].ux = self.solved_u[k]
        # For nodal forces/reactions
        self.NF = self.F.copy()
        self.VU = [node[key] for node in self.U.values() for key in ("ux",)]
        nf_calc = np.dot(self.KG, self.VU)
        for k,ic in enumerate(range(self.getNumberOfNodes())):
            nd, var = self.index2key(ic, ("fx",))
            self.NF[nd][var] = nf_calc[k]
            self.nodes[ic].fx = nf_calc[k]

    def index2key(self,idx,opts=("ux",)):
        node = idx
        var = opts[0]
        return node,var



#~ *********************************************************************
#~ ****************************  BeamModel *****************************
#~ *********************************************************************    
class BeamModel(Model):
    """
    Model for finite element analysis
    """
    def __init__(self,name="Beam Model 01"):
        Model.__init__(self,name=name,mtype="beam")
        self.F = {} # Forces
        self.U = {} # Displacements
        self.dof = 2 # 2 DOF for beam element
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
            self.F[node.label] = {"fy":0.0, "m":0.0} # (fy, m)
            
    def buildDisplacementsVector(self):
        for node in self.nodes.values():
            self.U[node.label] = {"uy":np.nan, "ur":np.nan} # (uy, r)
    
    def addForce(self,node,force):
        if not(self.IS_KG_BUILDED): self.buildGlobalMatrix()
        self.F[node.label]["fy"] = force[0]
        node.fy = force[0]
        
    def addMoment(self,node,moment):
        if not(self.IS_KG_BUILDED): self.buildGlobalMatrix()
        self.F[node.label]["m"] = moment[0]
        node.m = moment[0]
        
    def addConstraint(self,node,**constraint):
        if not(self.IS_KG_BUILDED): self.buildGlobalMatrix()
        cs = constraint
        if cs.has_key('ux') and cs.has_key("uy") and cs.has_key('ur'): # 
            ux = cs.get('ux')
            uy = cs.get('uy')
            ur = cs.get('ur')
            node.setDisplacements(ux=ux, uy=uy, ur=ur)
            #~ print("Encastre")
            self.U[node.label]["uy"] = uy
            self.U[node.label]["ur"] = ur
        elif cs.has_key('ux') and cs.has_key("uy"): # 
            ux = cs.get('ux')
            uy = cs.get('uy')
            node.setDisplacements(ux=ux, uy=uy)
            #~ print("Fixed")
            self.U[node.label]["uy"] = uy
        elif cs.has_key('uy'):
            uy = cs.get('uy')
            node.setDisplacements(uy=uy)
            #~ print("Simple support")
            self.U[node.label]["uy"] = uy
        
    def solve(self):
        # Solve LS
        self.VU = [node[key] for node in self.U.values() for key in ("uy","ur")]
        self.VF = [node[key] for node in self.F.values() for key in ("fy","m")]
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
            if np.isnan(nd.uy):
                nd.uy = self.U[nd.label]["uy"]
            if np.isnan(nd.ur):
                nd.ur = self.U[nd.label]["ur"]
                    
        # For nodal forces/reactions
        self.NF = self.F.copy()
        self.VU = [node[key] for node in self.U.values() for key in ("uy","ur")]
        nf_calc = np.dot(self.KG, self.VU)
        for k in range(2*self.getNumberOfNodes()):
            nd, var = self.index2key(k, ("fy","m"))
            self.NF[nd][var] = nf_calc[k]
            cnlab = np.floor(k/float(self.dof))
            if var=="fy": 
                self.nodes[cnlab].fy = nf_calc[k]
            elif var=="m": 
                self.nodes[cnlab].m = nf_calc[k]
            
    def index2key(self,idx,opts=("uy","ur")):
        node = idx//2
        var = opts[0] if ((-1)**idx)==1 else opts[1]
        return node,var
        
    def plot_model(self):
        import matplotlib.pyplot as plt
        
        fig = plt.figure()
        ax = fig.add_subplot(111)
        
        for elm in self.getElements():
            ni,nj = elm.getNodes()
            xx = [ni.x, nj.x]
            yy = [ni.y, nj.y]
            ax.plot(xx, yy, "r.-")
            for nd in (ni,nj):
                if nd.fx > 0: self._draw_xforce(ax,nd.x,nd.y,1)
                if nd.fx < 0: self._draw_xforce(ax,nd.x,nd.y,-1)
                if nd.fy > 0: self._draw_yforce(ax,nd.x,nd.y,1)
                if nd.fy < 0: self._draw_yforce(ax,nd.x,nd.y,-1)
                if nd.ux == 0: self._draw_xconstraint(ax,nd.x,nd.y)
                if nd.uy == 0: self._draw_yconstraint(ax,nd.x,nd.y)
            
        ax.axis("equal")
        x0,x1,y0,y1 = self.rect_region()
        ax.set_xlim(x0,x1)
        ax.set_ylim(y0,y1)

    def _draw_xforce(self,axes,x,y,ddir=1):
        """
        Draw horizontal arrow -> Force in x-dir
        """
        dx, dy = self._calculate_arrow_size(), 0
        HW = dx/5.0
        HL = dx/3.0
        arrow_props = dict(head_width=HW, head_length=HL, fc='r', ec='r')
        axes.arrow(x, y, ddir*dx, dy, **arrow_props)
        
    def _draw_yforce(self,axes,x,y,ddir=1):
        """
        Draw vertical arrow -> Force in y-dir
        """
        dx,dy = 0, self._calculate_arrow_size()
        HW = dy/5.0
        HL = dy/3.0
        arrow_props = dict(head_width=HW, head_length=HL, fc='r', ec='r')
        axes.arrow(x, y, dx, ddir*dy, **arrow_props)
        
    def _draw_xconstraint(self,axes,x,y):
        axes.plot(x, y, "g<", markersize=10, alpha=0.6)
    
    def _draw_yconstraint(self,axes,x,y):
        axes.plot(x, y, "gv", markersize=10, alpha=0.6)
        
    def _calculate_arrow_size(self):
        x0,x1,y0,y1 = self.rect_region(factor=10)
        sf = 5e-2
        kfx = sf*(x1-x0)
        kfy = sf*(y1-y0)
        return np.mean([kfx,kfy])

    def rect_region(self,factor=7.0):
        nx,ny = [],[]
        for n in self.getNodes():
            nx.append(n.x)
            ny.append(n.y)
        xmn,xmx,ymn,ymx = min(nx),max(nx),min(ny),max(ny)
        kx = (xmx-xmn)/factor
        if ymx==0 and ymn==0:
            ky = 1.0/factor
        else:
            ky = (ymx-ymn)/factor
        return xmn-kx, xmx+kx, ymn-ky, ymx+ky
        
    def plot_disp(self):
        import matplotlib.pyplot as plt
        
        fig = plt.figure()
        ax = fig.add_subplot(111)
        
        df = 1000
        for elm in self.getElements():
            ni,nj = elm.getNodes()
            xx = [ni.x, nj.x]
            yy = [ni.y+ni.uy*df, nj.y+nj.uy*df]
            ax.plot(xx, yy, "ro--")
            
        ax.axis("equal")
        
    def plot_moment_diagram(self):
        import matplotlib.pyplot as plt
        
        fig = plt.figure()
        ax = fig.add_subplot(111)
        
        X,M = self._get_data_for_moment_diagram()
        ax.plot(X, M)
        
    def plot_shear_diagram(self):
        import matplotlib.pyplot as plt
        
        fig = plt.figure()
        ax = fig.add_subplot(111)
        
        X,S = self._get_data_for_shear_diagram()
        ax.plot(X, S)
        
    def _get_data_for_moment_diagram(self):
        cx = 0
        X, M = [], []
        for el in self.getElements():
            L = el.L
            X = np.concatenate((X, np.array([cx, cx+L])))
            mel = el.m.squeeze()
            mel[0] = - mel[0]
            M = np.concatenate((M, mel))
            cx = cx + L
        return X, M
        
    def _get_data_for_shear_diagram(self):
        cx = 0
        X, S = [], []
        for el in self.getElements():
            L = el.L # element length
            X = np.concatenate((X, np.array([cx, cx+L])))
            fel = el.fy.squeeze()
            fel[-1] = - fel[-1]
            S = np.concatenate((S, fel))
            cx = cx + L
        return X, S
    
    def show(self):
        import matplotlib.pyplot as plt
        plt.show()


#~ *********************************************************************
#~ ****************************  LinearTriangleModel *******************
#~ *********************************************************************    
class LinearTriangleModel(Model):
    """
    Model for finite element analysis
    """
    def __init__(self,name="LT Model 01"):
        Model.__init__(self,name=name,mtype="triangle")
        self.F = {} # Forces
        self.U = {} # Displacements
        self.dof = 2 # 2 DOF for triangle element (per node)
        self.IS_KG_BUILDED = False
        
    def buildGlobalMatrix(self):
        """
        Build global matrix -> KG
        """
        msz = (self.dof)*self.getNumberOfNodes()
        self.KG = np.zeros((msz,msz))
        for element in self.elements.values():
            ku = element.getElementStiffness()
            n1,n2,n3 = element.getNodes()
            i, j, m = n1.label, n2.label, n3.label
            self.KG[2*i,2*i] += ku[0,0]
            self.KG[2*i,2*i+1] += ku[0,1]
            self.KG[2*i,2*j] += ku[0,2]
            self.KG[2*i,2*j+1] += ku[0,3]
            self.KG[2*i,2*m] += ku[0,4]
            self.KG[2*i,2*m+1] += ku[0,5]
            self.KG[2*i+1,2*i] += ku[1,0]
            self.KG[2*i+1,2*i+1] += ku[1,1]
            self.KG[2*i+1,2*j] += ku[1,2]
            self.KG[2*i+1,2*j+1] += ku[1,3]
            self.KG[2*i+1,2*m] += ku[1,4]
            self.KG[2*i+1,2*m+1] += ku[1,5]
            self.KG[2*j,2*i] += ku[2,0]
            self.KG[2*j,2*i+1] += ku[2,1]
            self.KG[2*j,2*j] += ku[2,2]
            self.KG[2*j,2*j+1] += ku[2,3]
            self.KG[2*j,2*m] += ku[2,4]
            self.KG[2*j,2*m+1] += ku[2,5]
            self.KG[2*j+1,2*i] += ku[3,0]
            self.KG[2*j+1,2*i+1] += ku[3,1]
            self.KG[2*j+1,2*j] += ku[3,2]
            self.KG[2*j+1,2*j+1] += ku[3,3]
            self.KG[2*j+1,2*m] += ku[3,4]
            self.KG[2*j+1,2*m+1] += ku[3,5]
            self.KG[2*m,2*i] += ku[4,0]
            self.KG[2*m,2*i+1] += ku[4,1]
            self.KG[2*m,2*j] += ku[4,2]
            self.KG[2*m,2*j+1] += ku[4,3]
            self.KG[2*m,2*m] += ku[4,4]
            self.KG[2*m,2*m+1] += ku[4,5]
            self.KG[2*m+1,2*i] += ku[5,0]
            self.KG[2*m+1,2*i+1] += ku[5,1]
            self.KG[2*m+1,2*j] += ku[5,2]
            self.KG[2*m+1,2*j+1] += ku[5,3]
            self.KG[2*m+1,2*m] += ku[5,4]
            self.KG[2*m+1,2*m+1] += ku[5,5]
            
        self.buildForcesVector()
        self.buildDisplacementsVector()
        self.IS_KG_BUILDED = True
    
    def buildForcesVector(self):
        for node in self.nodes.values():
            self.F[node.label] = {"fx":0.0, "fy":0.0} # (fy, m)
            
    def buildDisplacementsVector(self):
        for node in self.nodes.values():
            self.U[node.label] = {"ux":np.nan, "uy":np.nan} # (uy, r)
    
    def addForce(self,node,force):
        if not(self.IS_KG_BUILDED): self.buildGlobalMatrix()
        self.F[node.label]["fx"] = force[0]
        self.F[node.label]["fy"] = force[1]
        node.fx = force[0]
        node.fy = force[1]
        
    def addMoment(self,node,moment):
        pass
        
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

        _x,_y = [],[]
        patches = []
        for k,elm in enumerate(self.getElements()):
            _x,_y,_ux,_uy = [],[],[],[]
            for nd in elm.nodes:
                if nd.fx != 0: self._draw_xforce(ax,nd.x,nd.y)
                if nd.fy != 0: self._draw_yforce(ax,nd.x,nd.y)
                if nd.ux == 0 and nd.uy == 0: self._draw_xyconstraint(ax,nd.x,nd.y)
                _x.append(nd.x)
                _y.append(nd.y)
            polygon = Polygon(zip(_x,_y), True)
            patches.append(polygon)

        pc = PatchCollection(patches, color="#7CE7FF", edgecolor="k", alpha=0.4)
        ax.add_collection(pc)
        x0,x1,y0,y1 = self.rect_region()
        ax.set_xlim(x0,x1)
        ax.set_ylim(y0,y1)
        ax.set_title("Model %s"%(self.name))
        ax.set_aspect("equal")

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
        
    def _get_tri(self):
        import matplotlib.tri as tri
        
        _x,_y = [],[]
        for n in self.getNodes():
            _x.append(n.x)
            _y.append(n.y)
            
        tg = []
        for e in self.getElements():
            ni,nj,nm = e.getNodes()
            tg.append([ni.label, nj.label, nm.label])
            
        tr = tri.Triangulation(_x,_y, triangles=tg)
        return tr
        
    def plot_nsol(self,var="ux"):
        import matplotlib.pyplot as plt
        import numpy as np
        
        fig = plt.figure()
        ax = fig.add_subplot(111)
        
        solutions = {
             "ux": [n.ux for n in self.getNodes()],
             "uy": [n.uy for n in self.getNodes()],
             "usum": [np.sqrt(n.ux**2 + n.uy**2) for n in self.getNodes()],
             "sxx": [n.sx for n in self.getNodes()],
             "syy": [n.sy for n in self.getNodes()],
             "sxy": [n.sxy for n in self.getNodes()],
             "seqv": [n.seqv for n in self.getNodes()],
             "exx": [n.ex for n in self.getNodes()],
             "eyy": [n.ey for n in self.getNodes()],
             "exy": [n.exy for n in self.getNodes()]
             }
        
        tr = self._get_tri()
        try:
            fsol = solutions.get(var)
        except:
            return None
        if isinstance(fsol,list): fsol = np.array(fsol)
        tp = ax.tricontourf(tr, fsol, cmap="jet")
        fig.colorbar(tp)
        x0,x1,y0,y1 = self.rect_region()
        ax.set_xlim(x0,x1)
        ax.set_ylim(y0,y1)
        ax.set_aspect("equal")
        ax_title = "{0} (Max:{1}, Min:{2})".format(var,fsol.max(),fsol.min())
        ax.set_title(ax_title)

    def plot_esol(self,var="ux"):
        import matplotlib.pyplot as plt
        import numpy as np
        from matplotlib.patches import Polygon
        from matplotlib.collections import PatchCollection
        
        fig = plt.figure()
        ax = fig.add_subplot(111)

        _x,_y = [],[]
        patches = []
        for k,elm in enumerate(self.getElements()):
            _x,_y,_ux,_uy = [],[],[],[]
            for nd in elm.nodes:
                _x.append(nd.x)
                _y.append(nd.y)
            polygon = Polygon(zip(_x,_y), True)
            patches.append(polygon)
            
        pc = PatchCollection(patches, cmap="jet", alpha=1)
        solutions = {
             "sxx": [e.sx for e in self.getElements()],
             "syy": [e.sy for e in self.getElements()],
             "sxy": [e.sxy for e in self.getElements()],
             "exx": [e.ex for e in self.getElements()],
             "eyy": [e.ey for e in self.getElements()],
             "exy": [e.exy for e in self.getElements()]
             }
        fsol = np.array(solutions.get(var.lower()))
        pc.set_array(fsol)
        ax.add_collection(pc)
        plt.colorbar(pc)
        x0,x1,y0,y1 = self.rect_region()
        ax.set_xlim(x0,x1)
        ax.set_ylim(y0,y1)
        ax.set_aspect("equal")
        ax_title = "{0} (Max:{1}, Min:{2})".format(var,fsol.max(),fsol.min())
        ax.set_title(ax_title)
        
    def show(self):
        """
        Show matplotlib plots
        """
        import matplotlib.pyplot as plt
        plt.show()
    
    def calculate_deformed_factor(self):
        x0,x1,y0,y1 = self.rect_region()
        ux = np.array([n.ux for n in self.getNodes()])
        uy = np.array([n.uy for n in self.getNodes()])
        sf = 1.5e-2
        kfx = sf*(x1-x0)/ux.max()
        kfy = sf*(y1-y0)/uy.max()
        return np.mean([kfx,kfy])
                
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
