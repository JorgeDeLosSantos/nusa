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

class LinearTriangle(Element):
    """
    Linear triangle element for finite element analysis
    
    *nodes* : :class:`~nusa.core.Node`
        Connectivity for element
    
    *E* : float
        Young's modulus
        
    *nu* : float
        Poisson ratio
        
    *t* : float
        Thickness
    
    Example::
        n1 = Node((0,0))
        n2 = Node((0.5,0))
        n3 = Node((0.5,0.25))
        e1 = LinearTriangle((n1,n2,n3),210e9, 0.3, 0.025)
    """
    def __init__(self,nodes,E,nu,t):
        Element.__init__(self,etype="triangle")
        self.nodes = nodes
        self.E = E
        self.nu = nu
        self.t = t
        self._sx = 0
        self._sy = 0
        self._sxy = 0
        
    @property
    def sx(self):
        _sx,_sy,_sxy = self.getElementStresses()
        self._sx = _sx
        return self._sx
    
    @sx.setter
    def sx(self,val):
        self._sx = val
        
    @property
    def sy(self):
        _sx,_sy,_sxy = self.getElementStresses()
        self._sy = _sy
        return self._sy
    
    @sy.setter
    def sy(self,val):
        self._sy = val
    
    @property
    def sxy(self):
        _sx,_sy,_sxy = self.getElementStresses()
        self._sxy = _sxy
        return self._sxy
    
    @sxy.setter
    def sxy(self,val):
        self._sxy = val
    
    @property
    def ex(self):
        ex,ey,exy = self.getElementStrains()
        self._ex = ex
        return self._ex
    
    @ex.setter
    def ex(self,val):
        self._ex = val
    
    @property
    def D(self):
        """
        Constitutive matrix 
        
        Currently only plane stress supported
        """
        nu, E = self.nu, self.E
        D = (E/(1-nu**2))*np.array([[1, nu, 0],
                            [nu, 1, 0],
                            [0, 0, (1-nu)/2]
                            ])
        return D
    
    @property
    def B(self):
        ni, nj, nm = self.nodes
        A = self.A
        betai = nj.y - nm.y
        betaj = nm.y - ni.y
        betam = ni.y - nj.y
        gammai = nm.x - nj.x
        gammaj = ni.x - nm.x
        gammam = nj.x - ni.x
        B = (1/(2*A))*np.array([[betai, 0, betaj, 0, betam, 0],
                                 [0, gammai, 0, gammaj, 0, gammam],
                                 [gammai, betai, gammaj, betaj, gammam, betam]
                                 ])
        return B
    
    @property
    def A(self):
        n1, n2, n3 = self.nodes
        xi, yi = n1.x, n1.y
        xj, yj = n2.x, n2.y
        xm, ym = n3.x, n3.y
        A = (xi*(yj-ym) + xj*(ym-yi) + xm*(yi-yj))/2
        return A
    
    def getElementStiffness(self):
        """
        Get stiffness matrix for this element
        """
        ni, nj, nm = self.nodes
        A, nu, t, E = self.A, self.nu, self.t, self.E
        B, D = self.B, self.D
        return t*A*np.dot(np.dot(B.T,D),B)
        
    def getElementStresses(self):
        ni, nj, nm = self.nodes
        A, nu, t, E = self.A, self.nu, self.t, self.E
        u = np.array([ni.ux,ni.uy,nj.ux,nj.uy,nm.ux,nm.uy])
        B, D = self.B, self.D
        return np.dot(np.dot(D,B),u)
        
    def getElementStrains(self):
        ni, nj, nm = self.nodes
        A, nu, t, E = self.A, self.nu, self.t, self.E 
        u = np.array([ni.ux,ni.uy,nj.ux,nj.uy,nm.ux,nm.uy])
        B = self.B
        return np.dot(B,u)
        
    def getNodes(self):
        return self.nodes


#~ *********************************************************************
#~ ****************************  BeamModel *****************************
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
       
    def draw_bcs(self,ax):
		pass
        
    def _draw_xforce(self,axes,x,y):
		dx, dy = self._calculate_arrow_size(), 0
		HW = dx/5.0
		HL = dx/3.0
		arrow_props = dict(head_width=HW, head_length=HL, fc='r', ec='r')
		axes.arrow(x, y, dx, dy, **arrow_props)
        
    def _draw_yforce(self,axes,x,y):
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
        
    def plot_ux(self):
        import matplotlib.pyplot as plt
        import numpy as np
        
        fig = plt.figure()
        ax = fig.add_subplot(111)
        
        tr = self._get_tri()
        ux = [n.ux for n in self.getNodes()]
        tp = ax.tricontourf(tr, ux, cmap="jet")
        fig.colorbar(tp)
        x0,x1,y0,y1 = self.rect_region()
        ax.set_xlim(x0,x1)
        ax.set_ylim(y0,y1)
        ax.set_aspect("equal")
        ax.set_title("UX")
        
    def plot_uy(self):
        import matplotlib.pyplot as plt
        import numpy as np
        
        fig = plt.figure()
        ax = fig.add_subplot(111)
        
        tr = self._get_tri()
        ux = [n.uy for n in self.getNodes()]
        tp = ax.tricontourf(tr, ux, cmap="jet")
        fig.colorbar(tp)
        x0,x1,y0,y1 = self.rect_region()
        ax.set_xlim(x0,x1)
        ax.set_ylim(y0,y1)
        ax.set_aspect("equal")
        ax.set_title("UY")
        
    def plot_usum(self):
        import matplotlib.pyplot as plt
        import numpy as np
        
        fig = plt.figure()
        ax = fig.add_subplot(111)
        
        tr = self._get_tri()
        ux = [np.sqrt(n.ux**2 + n.uy**2) for n in self.getNodes()]
        tp = ax.tricontourf(tr, ux, cmap="jet")
        fig.colorbar(tp)
        x0,x1,y0,y1 = self.rect_region()
        ax.set_xlim(x0,x1)
        ax.set_ylim(y0,y1)
        ax.set_aspect("equal")
        ax.set_title("USUM")
        
    def plot_sxx(self):
        import matplotlib.pyplot as plt
        import numpy as np
        
        fig = plt.figure()
        ax = fig.add_subplot(111)
        
        tr = self._get_tri()
        sxx = [n.sx for n in self.getNodes()]
        tp = ax.tricontourf(tr, sxx, cmap="jet")
        fig.colorbar(tp)
        x0,x1,y0,y1 = self.rect_region()
        ax.set_xlim(x0,x1)
        ax.set_ylim(y0,y1)
        ax.set_aspect("equal")
        ax.set_title("SXX")
        
    def plot_syy(self):
        import matplotlib.pyplot as plt
        import numpy as np
        
        fig = plt.figure()
        ax = fig.add_subplot(111)
        
        tr = self._get_tri()
        syy = [n.sy for n in self.getNodes()]
        tp = ax.tricontourf(tr, syy, cmap="jet")
        fig.colorbar(tp)
        x0,x1,y0,y1 = self.rect_region()
        ax.set_xlim(x0,x1)
        ax.set_ylim(y0,y1)
        ax.set_aspect("equal")
        ax.set_title("SYY")
        
    def plot_sxy(self):
        import matplotlib.pyplot as plt
        import numpy as np
        
        fig = plt.figure()
        ax = fig.add_subplot(111)
        
        tr = self._get_tri()
        sxy = [n.sxy for n in self.getNodes()]
        tp = ax.tricontourf(tr, sxy, cmap="jet")
        fig.colorbar(tp)
        x0,x1,y0,y1 = self.rect_region()
        ax.set_xlim(x0,x1)
        ax.set_ylim(y0,y1)
        ax.set_aspect("equal")
        ax.set_title("SXY")
        
    def plot_seqv(self):
        import matplotlib.pyplot as plt
        import numpy as np
        
        fig = plt.figure()
        ax = fig.add_subplot(111)
        
        tr = self._get_tri()
        sxx = np.array([n.sx for n in self.getNodes()])
        syy = np.array([n.sy for n in self.getNodes()])
        sxy = np.array([n.sx for n in self.getNodes()])
        seqv = np.sqrt(sxx**2 - sxx*syy + syy**2 + 3*sxy**2)
        tp = ax.tricontourf(tr, seqv, cmap="jet")
        fig.colorbar(tp)
        x0,x1,y0,y1 = self.rect_region()
        ax.set_xlim(x0,x1)
        ax.set_ylim(y0,y1)
        ax.set_aspect("equal")
        ax.set_title("SEQV")

    def plot_esxx(self):
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
        strss = [e.sx for e in self.getElements()]
        pc.set_array(np.array(strss))
        ax.add_collection(pc)
        plt.colorbar(pc)
        x0,x1,y0,y1 = self.rect_region()
        ax.set_xlim(x0,x1)
        ax.set_ylim(y0,y1)
        ax.set_aspect("equal")
        ax.set_title("SXX")
        
    def plot_esyy(self):
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
        strss = [e.sy for e in self.getElements()]
        pc.set_array(np.array(strss))
        ax.add_collection(pc)
        plt.colorbar(pc)
        x0,x1,y0,y1 = self.rect_region()
        ax.set_xlim(x0,x1)
        ax.set_ylim(y0,y1)
        ax.set_aspect("equal")
        ax.set_title("SYY")
        
    def plot_eseqv(self):
        import matplotlib.pyplot as plt
        import numpy as np
        from matplotlib.patches import Polygon
        from matplotlib.collections import PatchCollection
        
        fig = plt.figure()
        ax = fig.add_subplot(111)

        _x,_y = [],[]
        patches = []
        kf = self.calculate_deformed_factor()
        for k,elm in enumerate(self.getElements()):
            _x,_y,_ux,_uy = [],[],[],[]
            for nd in elm.nodes:
                _x.append(nd.x + nd.ux*kf)
                _y.append(nd.y + nd.uy*kf)
            polygon = Polygon(zip(_x,_y), True)
            patches.append(polygon)

        pc = PatchCollection(patches, cmap="jet", alpha=1)
        sxx = np.array([e.sx for e in self.getElements()])
        syy = np.array([e.sy for e in self.getElements()])
        sxy = np.array([e.sxy for e in self.getElements()])
        seqv = np.sqrt(sxx**2 - sxx*syy + syy**2 + 3*sxy**2)
        pc.set_array(seqv)
        ax.add_collection(pc)
        plt.colorbar(pc)
        x0,x1,y0,y1 = self.rect_region()
        ax.set_xlim(x0,x1)
        ax.set_ylim(y0,y1)
        ax.set_aspect("equal")
        ax.set_title("SEQV")
        
    def plot_eex(self):
        import matplotlib.pyplot as plt
        import numpy as np
        from matplotlib.patches import Polygon
        from matplotlib.collections import PatchCollection
        
        fig = plt.figure()
        ax = fig.add_subplot(111)

        _x,_y = [],[]
        patches = []
        kf = self.calculate_deformed_factor()
        for k,elm in enumerate(self.getElements()):
            _x,_y,_ux,_uy = [],[],[],[]
            for nd in elm.nodes:
                _x.append(nd.x + nd.ux*kf)
                _y.append(nd.y + nd.uy*kf)
            polygon = Polygon(zip(_x,_y), True)
            patches.append(polygon)

        pc = PatchCollection(patches, cmap="jet", alpha=1)
        ex = np.array([e.ex for e in self.getElements()])
        pc.set_array(ex)
        ax.add_collection(pc)
        plt.colorbar(pc)
        x0,x1,y0,y1 = self.rect_region()
        ax.set_xlim(x0,x1)
        ax.set_ylim(y0,y1)
        ax.set_aspect("equal")
        ax.set_title("EX")
        
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
    n1 = Node((0,0))
    n2 = Node((0.5,0))
    n3 = Node((0.5,0.25))
    n4 = Node((0,0.25))
    e1 = LinearTriangle((n1,n3,n4),210e6,0.3,0.025)
    e2 = LinearTriangle((n1,n2,n3),210e6,0.3,0.025)
    m = LinearTriangleModel()
    for node in (n1,n2,n3,n4): m.addNode(node)
    for elm in (e1,e2): m.addElement(elm)
    m.addConstraint(n1, ux=0, uy=0)
    m.addConstraint(n4, ux=0, uy=0)
    m.addForce(n2, (9.375,0))
    m.addForce(n3, (9.375,0))
    m.solve()
