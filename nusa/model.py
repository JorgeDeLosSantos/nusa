# ***********************************
#  Author: Pedro Jorge De Los Santos     
#  E-mail: delossantosmfq@gmail.com 
#  Blog: numython.github.io
#  License: MIT License
# ***********************************
import re
import numpy as np
import numpy.linalg as la
import nusa.templates as tmp
import matplotlib.pyplot as plt
from .core import Model


def _partition_system(K, F, U):
    """Build the reduced linear system for prescribed displacements."""
    U = np.asarray(U, dtype=float)
    F = np.asarray(F, dtype=float)

    known = np.flatnonzero(~np.isnan(U))
    unknown = np.flatnonzero(np.isnan(U))

    Kuu = K[np.ix_(unknown, unknown)]
    Fu = F[unknown].copy()
    if known.size:
        Kuk = K[np.ix_(unknown, known)]
        Fu -= np.dot(Kuk, U[known])

    return known.tolist(), unknown.tolist(), Kuu, Fu


def _element_dof_indices(model, element):
    """Return global DOF indices using model-owned contiguous node indices."""
    indices = []
    for node in element.get_nodes():
        node_index = model._get_node_index(node)
        base = model.dof * node_index
        indices.extend(base + component for component in range(model.dof))
    return indices


def _assemble_global_stiffness(model):
    """Assemble the dense global stiffness matrix from element matrices."""
    matrix_size = model.dof * model.n_nodes
    model.KG = np.zeros((matrix_size, matrix_size))

    for element in model.elements:
        element_stiffness = element.get_element_stiffness()
        global_dofs = _element_dof_indices(model, element)
        model.KG[np.ix_(global_dofs, global_dofs)] += element_stiffness

    model.build_forces_vector()
    model.build_displacements_vector()
    model._restore_input_state()
    model.IS_KG_BUILDED = True


def _solve_model_system(
    model,
    displacement_keys,
    force_keys,
    *,
    allow_lstsq=False,
):
    """Solve a model using its assembled stiffness matrix and DOF dictionaries."""
    if not model.IS_KG_BUILDED:
        model.build_global_matrix()

    model.VU = [
        node[key]
        for node in model.U.values()
        for key in displacement_keys
    ]
    model.VF = [
        node[key]
        for node in model.F.values()
        for key in force_keys
    ]

    _, unknown, model.K2S, model.F2S = _partition_system(
        model.KG, model.VF, model.VU
    )

    if allow_lstsq:
        try:
            model.solved_u = la.solve(model.K2S, model.F2S)
        except:
            print("Solved using LSTSQ")
            model.solved_u = la.lstsq(model.K2S, model.F2S)[0]
    else:
        model.solved_u = la.solve(model.K2S, model.F2S)

    for value, dof_index in zip(model.solved_u, unknown):
        node_index, variable = model.index2key(dof_index, displacement_keys)
        model.U[node_index][variable] = value
        setattr(model.nodes[node_index], variable, value)

    model.NF = model.F.copy()
    model.VU = [
        node[key]
        for node in model.U.values()
        for key in displacement_keys
    ]
    nodal_forces = np.dot(model.KG, model.VU)

    for dof_index, value in enumerate(nodal_forces):
        node_index, variable = model.index2key(dof_index, force_keys)
        model.NF[node_index][variable] = value
        setattr(model.nodes[node_index], variable, value)

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

    def build_global_matrix(self):
        _assemble_global_stiffness(self)
        
    def build_forces_vector(self):
        for node in self.nodes:
            self.F[self._get_node_index(node)] = {"fx":0, "fy":0}
        
    def build_displacements_vector(self):
        for node in self.nodes:
            self.U[self._get_node_index(node)] = {"ux":np.nan, "uy":np.nan}
        
    def add_force(self,node,force):
        if not(self.IS_KG_BUILDED): self.build_global_matrix()
        self.F[self._get_node_index(node)]["fx"] = force[0]
        self._record_applied_forces(node, fx=force[0])
        
    def add_constraint(self,node,**constraint):
        """
        Only displacement in x-dir 
        """
        if not(self.IS_KG_BUILDED): self.build_global_matrix()
        if "ux" in constraint:
            ux = constraint.get("ux")
            node.set_displacements(ux=ux)
            self.U[self._get_node_index(node)]["ux"] = ux
            self._record_prescribed_displacements(node, ux=ux)
        
    def solve(self):
        _solve_model_system(self, ("ux",), ("fx",))
            
    def index2key(self,idx,opts=("ux",)):
        node = idx
        var = opts[0]
        return node,var

    def simple_report(self,report_type="print",fname="nusa_rpt.txt"):
        from .templates import SPRING_SIMPLE_REPORT
        options = {"headers":"firstrow",
                   "tablefmt":"rst",
                   "numalign":"right"}
        _str = SPRING_SIMPLE_REPORT.format(
                model_name=self.name,
                nodes=self.n_nodes,
                elements=self.n_elements,
                nodal_displacements=self._get_ndisplacements(options),
                nodal_forces=self._get_nforces(options),
                element_forces=self._get_eforces(options),
                nodes_info=self._get_nodes_info(options),
                elements_info=self._get_elements_info(options))
        if report_type=="print": print(_str)
        elif report_type=="write": self._write_report(_str, fname)
        elif report_type=="string": return _str
        else: return _str

    def _get_eforces(self,options):
        from tabulate import tabulate
        F = [["Element","F"]]
        for elm in self.elements:
            F.append([elm.label+1, elm.fx])
        return tabulate(F, **options)
        


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
        
    def build_forces_vector(self):
        """
        Build forces vector, where each node has a dict with "fx" and "fy" keys, but only "fx" is used for bar model
        """
        for node in self.nodes:
            self.F[self._get_node_index(node)] = {"fx":0, "fy":0}
        
    def build_global_matrix(self):
        _assemble_global_stiffness(self)
        
    def build_displacements_vector(self):
        for node in self.nodes:
            self.U[self._get_node_index(node)] = {"ux":np.nan, "uy":np.nan}
        
    def add_force(self,node,force):
        if not(self.IS_KG_BUILDED): self.build_global_matrix()
        self.F[self._get_node_index(node)]["fx"] = force[0]
        self._record_applied_forces(node, fx=force[0])
        
    def add_constraint(self,node,**constraint):
        if not(self.IS_KG_BUILDED): self.build_global_matrix()
        if "ux" in constraint:
            ux = constraint.get('ux')
            node.set_displacements(ux=ux)
            self.U[self._get_node_index(node)]["ux"] = ux
            self._record_prescribed_displacements(node, ux=ux)
        
    def solve(self):
        _solve_model_system(self, ("ux",), ("fx",))

    def index2key(self,idx,opts=("ux",)):
        node = idx
        var = opts[0]
        return node,var



#~ *********************************************************************
#~ ****************************  TrussModel ****************************
#~ *********************************************************************
class TrussModel(Model):
    """
    Truss model for finite element analysis
    """
    def __init__(self,name="Truss Model 01"):
        Model.__init__(self,name=name,mtype="truss")
        self.F = {} # Forces
        self.U = {} # Displacements
        self.dof = 2 # 2 DOF for truss element
        self.IS_KG_BUILDED = False
        
    def build_global_matrix(self):
        _assemble_global_stiffness(self)
        
    def build_forces_vector(self):
        for node in self.nodes:
            self.F[self._get_node_index(node)] = {"fx":0, "fy":0}
        
    def build_displacements_vector(self):
        for node in self.nodes:
            self.U[self._get_node_index(node)] = {"ux":np.nan, "uy":np.nan}
    
    def add_force(self,node,force):
        if not(self.IS_KG_BUILDED): self.build_global_matrix()
        self.F[self._get_node_index(node)]["fx"] = force[0]
        self.F[self._get_node_index(node)]["fy"] = force[1]
        self._record_applied_forces(node, fx=force[0], fy=force[1])
        node.fx = force[0]
        node.fy = force[1]
        
    def add_constraint(self,node,**constraint):
        if not(self.IS_KG_BUILDED): self.build_global_matrix()
        cs = constraint
        if "ux" in cs and "uy" in cs: #
            ux = cs.get('ux')
            uy = cs.get('uy')
            node.set_displacements(ux=ux, uy=uy) # eqv to node.ux = ux, node.uy = uy
            self.U[self._get_node_index(node)]["ux"] = ux
            self.U[self._get_node_index(node)]["uy"] = uy
            self._record_prescribed_displacements(node, ux=ux, uy=uy)
        elif "ux" in cs:
            ux = cs.get('ux')
            node.set_displacements(ux=ux)
            self.U[self._get_node_index(node)]["ux"] = ux
            self._record_prescribed_displacements(node, ux=ux)
        elif "uy" in cs:
            uy = cs.get('uy')
            node.set_displacements(uy=uy)
            self.U[self._get_node_index(node)]["uy"] = uy
            self._record_prescribed_displacements(node, uy=uy)
        else: pass # todo
        
    def solve(self):
        _solve_model_system(self, ("ux", "uy"), ("fx", "fy"))
                
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
        
        fig = plt.figure()
        ax = fig.add_subplot(111)
        
        for elm in self.elements:
            ni, nj = elm.get_nodes()
            ax.plot([ni.x,nj.x],[ni.y,nj.y],"b-")
            for nd in (ni,nj):
                if nd.fx > 0: self._draw_xforce(ax,nd.x,nd.y,1)
                if nd.fx < 0: self._draw_xforce(ax,nd.x,nd.y,-1)
                if nd.fy > 0: self._draw_yforce(ax,nd.x,nd.y,1)
                if nd.fy < 0: self._draw_yforce(ax,nd.x,nd.y,-1)
                if nd.ux == 0: self._draw_xconstraint(ax,nd.x,nd.y)
                if nd.uy == 0: self._draw_yconstraint(ax,nd.x,nd.y)
        
        x0,x1,y0,y1 = self.rect_region()
        plt.axis('equal')
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
        x0,x1,y0,y1 = self.rect_region(factor=50)
        sf = 5e-2
        kfx = sf*(x1-x0)
        kfy = sf*(y1-y0)
        return np.mean([kfx,kfy])
        
    def plot_deformed_shape(self,dfactor=1.0):
        import matplotlib.pyplot as plt
        fig = plt.figure()
        ax = fig.add_subplot(111)
        
        df = dfactor*self._calculate_deformed_factor()
        
        for elm in self.elements:
            ni,nj = elm.get_nodes()
            x, y = [ni.x,nj.x], [ni.y,nj.y]
            xx = [ni.x+ni.ux*df, nj.x+nj.ux*df]
            yy = [ni.y+ni.uy*df, nj.y+nj.uy*df]
            ax.plot(x,y,'bo-')
            ax.plot(xx,yy,'ro--')

        x0,x1,y0,y1 = self.rect_region()
        plt.axis('equal')
        ax.set_xlim(x0,x1)
        ax.set_ylim(y0,y1)
        
    def _calculate_deformed_factor(self):
        x0,x1,y0,y1 = self.rect_region()
        ux = np.abs(np.array([n.ux for n in self.nodes]))
        uy = np.abs(np.array([n.uy for n in self.nodes]))
        sf = 1.5e-2
        if ux.max()==0 and uy.max()!=0:
            kfx = sf*(y1-y0)/uy.max()
            kfy = sf*(y1-y0)/uy.max()
        if uy.max()==0 and ux.max()!=0:
            kfx = sf*(x1-x0)/ux.max()
            kfy = sf*(x1-x0)/ux.max()
        if ux.max()!=0 and uy.max()!=0:
            kfx = sf*(x1-x0)/ux.max()
            kfy = sf*(y1-y0)/uy.max()
        return np.mean([kfx,kfy])

    def show(self):
        import matplotlib.pyplot as plt
        plt.show()
        
    def rect_region(self,factor=7.0):
        nx,ny = [],[]
        for n in self.nodes:
            nx.append(n.x)
            ny.append(n.y)
        xmn,xmx,ymn,ymx = min(nx),max(nx),min(ny),max(ny)
        kx = (xmx-xmn)/factor
        ky = (ymx-ymn)/factor
        return xmn-kx, xmx+kx, ymn-ky, ymx+ky
        
    def simple_report(self,report_type="print",fname="nusa_rpt.txt"):
        from .templates import TRUSS_SIMPLE_REPORT
        options = {"headers":"firstrow",
                   "tablefmt":"rst",
                   "numalign":"right"}
        _str = TRUSS_SIMPLE_REPORT.format(
                model_name=self.name,
                nodes=self.n_nodes,
                elements=self.n_elements,
                nodal_displacements=self._get_ndisplacements(options),
                nodal_forces=self._get_nforces(options),
                element_forces=self._get_eforces(options),
                element_stresses=self._get_estresses(options),
                nodes_info=self._get_nodes_info(options),
                elements_info=self._get_elements_info(options))
        if report_type=="print": print(_str)
        elif report_type=="write": self._write_report(_str, fname)
        elif report_type=="string": return _str
        else: return _str
        
    def _write_report(self,txt,fname):
        fobj = open(fname,"w")
        fobj.write(txt)
        fobj.close()
        
    def _get_ndisplacements(self,options):
        from tabulate import tabulate
        D = [["Node","UX","UY"]]
        for n in self.nodes:
            D.append([n.label,n.ux,n.uy])
        return tabulate(D, **options)
        
    def _get_nforces(self,options):
        from tabulate import tabulate
        F = [["Node","FX","FY"]]
        for n in self.nodes:
            F.append([n.label,n.fx,n.fy])
        return tabulate(F, **options)
        
    def _get_eforces(self,options):
        from tabulate import tabulate
        F = [["Element","F"]]
        for elm in self.elements:
            F.append([elm.label+1, elm.f])
        return tabulate(F, **options)
        
    def _get_estresses(self,options):
        from tabulate import tabulate
        S = [["Element","S"]]
        for elm in self.elements:
            S.append([elm.label+1, elm.s])
        return tabulate(S, **options)
    
    def _get_nodes_info(self,options):
        from tabulate import tabulate
        F = [["Node","X","Y"]]
        for n in self.nodes:
            F.append([n.label, n.x, n.y])
        return tabulate(F, **options)
    
    def _get_elements_info(self,options):
        from tabulate import tabulate
        S = [["Element","NI","NJ"]]
        for elm in self.elements:
            ni, nj = elm.get_nodes()
            S.append([elm.label+1, ni.label, nj.label])
        return tabulate(S, **options)



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
        
    def build_global_matrix(self):
        _assemble_global_stiffness(self)
    
    def build_forces_vector(self):
        for node in self.nodes:
            self.F[self._get_node_index(node)] = {"fy":0.0, "m":0.0} # (fy, m)
            
    def build_displacements_vector(self):
        for node in self.nodes:
            self.U[self._get_node_index(node)] = {"uy":np.nan, "ur":np.nan} # (uy, r)
    
    def add_force(self,node,force):
        if not(self.IS_KG_BUILDED): self.build_global_matrix()
        self.F[self._get_node_index(node)]["fy"] = force[0]
        self._record_applied_forces(node, fy=force[0])
        node.fy = force[0]
        
    def add_moment(self,node,moment):
        if not(self.IS_KG_BUILDED): self.build_global_matrix()
        self.F[self._get_node_index(node)]["m"] = moment[0]
        self._record_applied_forces(node, m=moment[0])
        node.m = moment[0]
        
    def add_constraint(self,node,**constraint):
        if not(self.IS_KG_BUILDED): self.build_global_matrix()
        cs = constraint
        if "ux" in cs and "uy" in cs and "ur" in cs: # 
            ux = cs.get('ux')
            uy = cs.get('uy')
            ur = cs.get('ur')
            node.set_displacements(ux=ux, uy=uy, ur=ur)
            #~ print("Encastre")
            self.U[self._get_node_index(node)]["uy"] = uy
            self.U[self._get_node_index(node)]["ur"] = ur
            self._record_prescribed_displacements(node, uy=uy, ur=ur)
        elif "ux" in cs and "uy" in cs: # 
            ux = cs.get('ux')
            uy = cs.get('uy')
            node.set_displacements(ux=ux, uy=uy)
            #~ print("Fixed")
            self.U[self._get_node_index(node)]["uy"] = uy
            self._record_prescribed_displacements(node, uy=uy)
        elif "uy" in cs:
            uy = cs.get('uy')
            node.set_displacements(uy=uy)
            #~ print("Simple support")
            self.U[self._get_node_index(node)]["uy"] = uy
            self._record_prescribed_displacements(node, uy=uy)
        
    def solve(self):
        _solve_model_system(self, ("uy", "ur"), ("fy", "m"))
            
    def index2key(self,idx,opts=("uy","ur")):
        node = idx//2
        var = opts[0] if ((-1)**idx)==1 else opts[1]
        return node,var
        
    def plot_model(self):
        import matplotlib.pyplot as plt
        
        fig = plt.figure()
        ax = fig.add_subplot(111)
        
        for elm in self.elements:
            ni,nj = elm.get_nodes()
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
        for n in self.nodes:
            nx.append(n.x)
            ny.append(n.y)
        xmn,xmx,ymn,ymx = min(nx),max(nx),min(ny),max(ny)
        kx = (xmx-xmn)/factor
        if ymx==0 and ymn==0:
            ky = 1.0/factor
        else:
            ky = (ymx-ymn)/factor
        return xmn-kx, xmx+kx, ymn-ky, ymx+ky
        
    def plot_disp(self, df = 1000, **kwargs):
        fig = plt.figure()
        ax = fig.add_subplot(111)
        
        xx = []
        yy = []
        for elm in self.elements:
            ni,nj = elm.get_nodes()
            xx.append( ni.x )
            xx.append( nj.x )
            yy.append( ni.y+ni.uy*df )
            yy.append( nj.y+nj.uy*df )
        
        ax.plot(xx, yy, "ro--", **kwargs)
            
        ax.axis("equal")
        
    def plot_moment_diagram(self):
        import matplotlib.pyplot as plt
        
        fig = plt.figure()
        ax = fig.add_subplot(111)
        
        X,M = self._get_data_for_moment_diagram()
        ax.plot(X, M, "r")
        ax.fill_between(X, M, facecolor="#EE5B5B")
        
    def plot_shear_diagram(self):
        import matplotlib.pyplot as plt
        
        fig = plt.figure()
        ax = fig.add_subplot(111)
        
        X,S = self._get_data_for_shear_diagram()
        ax.plot(X, S, "b")
        ax.fill_between(X, S, facecolor="#559EE5")
        
    def _get_data_for_moment_diagram(self):
        cx = 0
        X, M = [], []
        for el in self.elements:
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
        for el in self.elements:
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
        
    def build_global_matrix(self):
        """Build global stiffness matrix."""
        _assemble_global_stiffness(self)

    def build_forces_vector(self):
        for node in self.nodes:
            self.F[self._get_node_index(node)] = {"fx":0.0, "fy":0.0} # (fy, m)
            
    def build_displacements_vector(self):
        for node in self.nodes:
            self.U[self._get_node_index(node)] = {"ux":np.nan, "uy":np.nan} # (uy, r)
    
    def add_force(self,node,force):
        if not(self.IS_KG_BUILDED): self.build_global_matrix()
        self.F[self._get_node_index(node)]["fx"] = force[0]
        self.F[self._get_node_index(node)]["fy"] = force[1]
        self._record_applied_forces(node, fx=force[0], fy=force[1])
        node.fx = force[0]
        node.fy = force[1]
        
    def add_moment(self,node,moment):
        pass
        
    def add_constraint(self,node,**constraint):
        if not(self.IS_KG_BUILDED): self.build_global_matrix()
        cs = constraint
        if "ux" in cs and "uy" in cs: # 
            ux = cs.get('ux')
            uy = cs.get('uy')
            node.set_displacements(ux=ux, uy=uy)
            self.U[self._get_node_index(node)]["ux"] = ux
            self.U[self._get_node_index(node)]["uy"] = uy
            self._record_prescribed_displacements(node, ux=ux, uy=uy)
        elif "uy" in cs:
            uy = cs.get('uy')
            node.set_displacements(uy=uy)
            self.U[self._get_node_index(node)]["uy"] = uy
            self._record_prescribed_displacements(node, uy=uy)
        
    def _check_nodes(self):
        for node in self.nodes:
            if node._elements == []: self.add_constraint(node, ux=0, uy=0)
        
    def solve(self):
        self._check_nodes()
        _solve_model_system(
            self,
            ("ux", "uy"),
            ("fx", "fy"),
            allow_lstsq=True,
        )
                
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
        for k,elm in enumerate(self.elements):
            _x,_y,_ux,_uy = [],[],[],[]
            for nd in elm.nodes:
                if nd.fx != 0: self._draw_xforce(ax,nd.x,nd.y)
                if nd.fy != 0: self._draw_yforce(ax,nd.x,nd.y)
                if nd.ux == 0 and nd.uy == 0: self._draw_xyconstraint(ax,nd.x,nd.y)
                _x.append(nd.x)
                _y.append(nd.y)
            polygon = Polygon(list(zip(_x,_y)))
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
        # ~ df = 1
        for n in self.nodes:
            _x.append(n.x)
            # ~ _x.append(n.x + n.ux*df)
            _y.append(n.y)
            # ~ _y.append(n.y + n.uy*df)
            
        tg = []
        for e in self.elements:
            ni,nj,nm = e.get_nodes()
            tg.append([
                self._get_node_index(ni),
                self._get_node_index(nj),
                self._get_node_index(nm),
            ])
            
        tr = tri.Triangulation(_x,_y, triangles=tg)
        return tr


    def plot_nsol(self,var="ux"):
        import matplotlib.pyplot as plt
        import numpy as np
        
        fig = plt.figure()
        ax = fig.add_subplot(111)
        
        solutions = {
             "ux": (n.ux for n in self.nodes),
             "uy": (n.uy for n in self.nodes),
             "usum": (np.sqrt(n.ux**2 + n.uy**2) for n in self.nodes),
             "sxx": (n.sx for n in self.nodes),
             "syy": (n.sy for n in self.nodes),
             "sxy": (n.sxy for n in self.nodes),
             "seqv": (n.seqv for n in self.nodes),
             "exx": (n.ex for n in self.nodes),
             "eyy": (n.ey for n in self.nodes),
             "exy": (n.exy for n in self.nodes)
             }
        
        tr = self._get_tri()
        try:
            fsol = list(solutions.get(var))
        except:
            return None
        if isinstance(fsol,list): fsol = np.array(fsol)
        tp = ax.tricontourf(tr, fsol, cmap="jet")
        fig.colorbar(tp)
        x0,x1,y0,y1 = self.rect_region()
        ax.set_xlim(x0,x1)
        ax.set_ylim(y0,y1)
        ax.set_aspect("equal")
        ax_title = "{0} (Max:{1:0.3e}, Min:{2:0.3e})".format(var,fsol.max(),fsol.min())
        ax.set_title(ax_title, fontsize=8)


    def plot_esol(self,var="ux"):
        import matplotlib.pyplot as plt
        import numpy as np
        from matplotlib.patches import Polygon
        from matplotlib.collections import PatchCollection
        
        fig = plt.figure()
        ax = fig.add_subplot(111)

        _x,_y = [],[]
        patches = []
        for k,elm in enumerate(self.elements):
            _x,_y,_ux,_uy = [],[],[],[]
            for nd in elm.nodes:
                _x.append(nd.x)
                _y.append(nd.y)
            polygon = Polygon(list(zip(_x,_y)))
            patches.append(polygon)
            
        pc = PatchCollection(patches, cmap="jet", alpha=1)
        solutions = {
             "sxx": (e.sx for e in self.elements),
             "syy": (e.sy for e in self.elements),
             "sxy": (e.sxy for e in self.elements),
             "exx": (e.ex for e in self.elements),
             "eyy": (e.ey for e in self.elements),
             "exy": (e.exy for e in self.elements)
             }
        fsol = np.array(list(solutions.get(var.lower())))
        pc.set_array(fsol)
        ax.add_collection(pc)
        fig.colorbar(pc)
        x0,x1,y0,y1 = self.rect_region()
        ax.set_xlim(x0,x1)
        ax.set_ylim(y0,y1)
        ax.set_aspect("equal")
        ax_title = "{0} (Max:{1:0.3e}, Min:{2:0.3e})".format(var,fsol.max(),fsol.min())
        ax.set_title(ax_title, fontsize=8)
        
    def show(self):
        """
        Show matplotlib plots
        """
        import matplotlib.pyplot as plt
        plt.show()
    
    def calculate_deformed_factor(self):
        x0,x1,y0,y1 = self.rect_region()
        ux = np.array([n.ux for n in self.nodes])
        uy = np.array([n.uy for n in self.nodes])
        sf = 1.5e-2
        kfx = sf*(x1-x0)/ux.max()
        kfy = sf*(y1-y0)/uy.max()
        return np.mean([kfx,kfy])
                
    def rect_region(self,factor=7.0):
        nx,ny = [],[]
        for n in self.nodes:
            nx.append(n.x)
            ny.append(n.y)
        xmn,xmx,ymn,ymx = min(nx),max(nx),min(ny),max(ny)
        kx = (xmx-xmn)/factor
        ky = (ymx-ymn)/factor
        return xmn-kx, xmx+kx, ymn-ky, ymx+ky




if __name__=='__main__':
    pass
