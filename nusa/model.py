# ***********************************
#  Author: Pedro Jorge De Los Santos     
#  E-mail: delossantosmfq@gmail.com 
#  Blog: numython.github.io
#  License: MIT License
# ***********************************
import re
import numpy as np
import matplotlib.pyplot as plt
from .core import Model


#~ *********************************************************************
#~ ****************************  SpringModel ***************************
#~ *********************************************************************

class SpringModel(Model):
    """
    Spring Model for finite element analysis
    """
    displacement_dofs = ("ux",)
    force_dofs = ("fx",)

    def __init__(self,name="Spring Model 01"):
        Model.__init__(self,name=name,mtype="spring")
        self.dof = 1 # 1 DOF per Node

    def add_force(self,node,force):
        values = self._validated_component_vector(force, self.force_dofs, "force")
        self._record_applied_forces(node, **values)
        
    def add_constraint(self,node,**constraint):
        """Prescribe spring-model displacement components."""
        values = self._validated_named_components(
            constraint, self.displacement_dofs, "constraint"
        )
        for variable, value in values.items():
            setattr(node, variable, value)
        if values:
            self._record_prescribed_displacements(node, **values)
        

#~ *********************************************************************
#~ ****************************  BarModel ******************************
#~ *********************************************************************
class BarModel(Model):
    """
    Bar model for finite element analysis
    """
    displacement_dofs = ("ux",)
    force_dofs = ("fx",)

    def __init__(self,name="Bar Model 01"):
        Model.__init__(self,name=name,mtype="bar")
        self.dof = 1 # 1 DOF for bar element (per node)
        
    def add_force(self,node,force):
        values = self._validated_component_vector(force, self.force_dofs, "force")
        self._record_applied_forces(node, **values)
        
    def add_constraint(self,node,**constraint):
        values = self._validated_named_components(
            constraint, self.displacement_dofs, "constraint"
        )
        for variable, value in values.items():
            setattr(node, variable, value)
        if values:
            self._record_prescribed_displacements(node, **values)
        
#~ *********************************************************************
#~ ****************************  TrussModel ****************************
#~ *********************************************************************
class TrussModel(Model):
    """
    Truss model for finite element analysis
    """
    displacement_dofs = ("ux", "uy")
    force_dofs = ("fx", "fy")

    def __init__(self,name="Truss Model 01"):
        Model.__init__(self,name=name,mtype="truss")
        self.dof = 2 # 2 DOF for truss element
        
    def add_force(self,node,force):
        values = self._validated_component_vector(force, self.force_dofs, "force")
        self._record_applied_forces(node, **values)
        
    def add_constraint(self,node,**constraint):
        values = self._validated_named_components(
            constraint, self.displacement_dofs, "constraint"
        )
        for variable, value in values.items():
            setattr(node, variable, value)
        if values:
            self._record_prescribed_displacements(node, **values)
        
    def plot_model(self, show_reactions=False):
        """
        Plot model geometry, applied loads, constraints, and optional reactions.
        """
        import matplotlib.pyplot as plt
        
        fig = plt.figure()
        ax = fig.add_subplot(111)
        
        for elm in self.elements:
            ni, nj = elm.nodes
            ax.plot([ni.x,nj.x],[ni.y,nj.y],"b-")

        for nd in self.nodes:
            applied = self.applied_load(nd)
            if applied["fx"] > 0: self._draw_xforce(ax,nd.x,nd.y,1)
            if applied["fx"] < 0: self._draw_xforce(ax,nd.x,nd.y,-1)
            if applied["fy"] > 0: self._draw_yforce(ax,nd.x,nd.y,1)
            if applied["fy"] < 0: self._draw_yforce(ax,nd.x,nd.y,-1)

            if show_reactions:
                reaction = self.reaction(nd)
                if reaction["fx"] > 0: self._draw_xforce(ax,nd.x,nd.y,1,reaction=True)
                if reaction["fx"] < 0: self._draw_xforce(ax,nd.x,nd.y,-1,reaction=True)
                if reaction["fy"] > 0: self._draw_yforce(ax,nd.x,nd.y,1,reaction=True)
                if reaction["fy"] < 0: self._draw_yforce(ax,nd.x,nd.y,-1,reaction=True)

            if nd.ux == 0: self._draw_xconstraint(ax,nd.x,nd.y)
            if nd.uy == 0: self._draw_yconstraint(ax,nd.x,nd.y)
        
        x0,x1,y0,y1 = self._rect_region()
        plt.axis('equal')
        ax.set_xlim(x0,x1)
        ax.set_ylim(y0,y1)

    def _draw_xforce(self,axes,x,y,ddir=1,reaction=False):
        """
        Draw horizontal applied-load or reaction arrow.
        """
        dx, dy = self._calculate_arrow_size(), 0
        HW = dx/5.0
        HL = dx/3.0
        color = 'b' if reaction else 'r'
        arrow_props = dict(head_width=HW, head_length=HL, fc=color, ec=color)
        axes.arrow(x, y, ddir*dx, dy, **arrow_props)
        
    def _draw_yforce(self,axes,x,y,ddir=1,reaction=False):
        """
        Draw vertical applied-load or reaction arrow.
        """
        dx,dy = 0, self._calculate_arrow_size()
        HW = dy/5.0
        HL = dy/3.0
        color = 'b' if reaction else 'r'
        arrow_props = dict(head_width=HW, head_length=HL, fc=color, ec=color)
        axes.arrow(x, y, dx, ddir*dy, **arrow_props)
        
    def _draw_xconstraint(self,axes,x,y):
        axes.plot(x, y, "g<", markersize=10, alpha=0.6)
    
    def _draw_yconstraint(self,axes,x,y):
        axes.plot(x, y, "gv", markersize=10, alpha=0.6)
        
    def _calculate_arrow_size(self):
        x0,x1,y0,y1 = self._rect_region(factor=50)
        sf = 5e-2
        kfx = sf*(x1-x0)
        kfy = sf*(y1-y0)
        return np.mean([kfx,kfy])
        
    def plot_deformed_shape(self, scale=1.0, **kwargs):
        """Compatibility wrapper around StaticResult.plot_deformed_shape()."""
        if not hasattr(self, "_last_result"):
            raise RuntimeError(
                "plot_deformed_shape() is available only after solve()"
            )
        return self._last_result.plot_deformed_shape(scale=scale, **kwargs)

    def _calculate_deformed_factor(self):
        x0,x1,y0,y1 = self._rect_region()
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
        
    def _rect_region(self,factor=7.0):
        nx,ny = [],[]
        for n in self.nodes:
            nx.append(n.x)
            ny.append(n.y)
        xmn,xmx,ymn,ymx = min(nx),max(nx),min(ny),max(ny)
        kx = (xmx-xmn)/factor
        ky = (ymx-ymn)/factor
        if ky == 0:
            ky = 1.0/factor
        return xmn-kx, xmx+kx, ymn-ky, ymx+ky
        

#~ *********************************************************************
#~ ****************************  BeamModel *****************************
#~ *********************************************************************    
class BeamModel(Model):
    """
    Model for finite element analysis
    """
    displacement_dofs = ("uy", "ur")
    force_dofs = ("fy", "m")

    def __init__(self,name="Beam Model 01"):
        Model.__init__(self,name=name,mtype="beam")
        self.dof = 2 # 2 DOF for beam element
        
    def add_force(self,node,force):
        values = self._validated_component_vector(force, ("fy",), "force")
        self._record_applied_forces(node, **values)
        
    def add_moment(self,node,moment):
        values = self._validated_component_vector(moment, ("m",), "moment")
        self._record_applied_forces(node, **values)
        
    def add_constraint(self,node,**constraint):
        values = self._validated_named_components(
            constraint, self.displacement_dofs, "constraint"
        )
        for variable, value in values.items():
            setattr(node, variable, value)
        if values:
            self._record_prescribed_displacements(node, **values)
        
    def plot_model(self, show_reactions=False):
        """Plot beam geometry, applied transverse loads, and optional reactions."""
        import matplotlib.pyplot as plt
        
        fig = plt.figure()
        ax = fig.add_subplot(111)
        
        for elm in self.elements:
            ni,nj = elm.nodes
            xx = [ni.x, nj.x]
            yy = [ni.y, nj.y]
            ax.plot(xx, yy, "r.-")

        for nd in self.nodes:
            applied = self.applied_load(nd)
            if applied["fy"] > 0: self._draw_yforce(ax,nd.x,nd.y,1)
            if applied["fy"] < 0: self._draw_yforce(ax,nd.x,nd.y,-1)

            if show_reactions:
                reaction = self.reaction(nd)
                if reaction["fy"] > 0: self._draw_yforce(ax,nd.x,nd.y,1,reaction=True)
                if reaction["fy"] < 0: self._draw_yforce(ax,nd.x,nd.y,-1,reaction=True)

            if nd.ux == 0: self._draw_xconstraint(ax,nd.x,nd.y)
            if nd.uy == 0: self._draw_yconstraint(ax,nd.x,nd.y)
            
        ax.axis("equal")
        x0,x1,y0,y1 = self._rect_region()
        ax.set_xlim(x0,x1)
        ax.set_ylim(y0,y1)

    def _draw_xforce(self,axes,x,y,ddir=1,reaction=False):
        """
        Draw horizontal applied-load or reaction arrow.
        """
        dx, dy = self._calculate_arrow_size(), 0
        HW = dx/5.0
        HL = dx/3.0
        color = 'b' if reaction else 'r'
        arrow_props = dict(head_width=HW, head_length=HL, fc=color, ec=color)
        axes.arrow(x, y, ddir*dx, dy, **arrow_props)
        
    def _draw_yforce(self,axes,x,y,ddir=1,reaction=False):
        """
        Draw vertical applied-load or reaction arrow.
        """
        dx,dy = 0, self._calculate_arrow_size()
        HW = dy/5.0
        HL = dy/3.0
        color = 'b' if reaction else 'r'
        arrow_props = dict(head_width=HW, head_length=HL, fc=color, ec=color)
        axes.arrow(x, y, dx, ddir*dy, **arrow_props)
        
    def _draw_xconstraint(self,axes,x,y):
        axes.plot(x, y, "g<", markersize=10, alpha=0.6)
    
    def _draw_yconstraint(self,axes,x,y):
        axes.plot(x, y, "gv", markersize=10, alpha=0.6)
        
    def _calculate_arrow_size(self):
        x0,x1,y0,y1 = self._rect_region(factor=10)
        sf = 5e-2
        kfx = sf*(x1-x0)
        kfy = sf*(y1-y0)
        return np.mean([kfx,kfy])

    def _rect_region(self,factor=7.0):
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
        
    def plot_deformed_shape(self, scale=1000, **kwargs):
        """Compatibility wrapper around StaticResult.plot_deformed_shape()."""
        if not hasattr(self, "_last_result"):
            raise RuntimeError(
                "plot_deformed_shape() is available only after solve()"
            )
        return self._last_result.plot_deformed_shape(scale=scale, **kwargs)

    def plot_moment_diagram(self, **kwargs):
        """Compatibility wrapper around StaticResult.plot_moment_diagram()."""
        if not hasattr(self, "_last_result"):
            raise RuntimeError(
                "plot_moment_diagram() is available only after solve()"
            )
        return self._last_result.plot_moment_diagram(**kwargs)

    def plot_shear_diagram(self, **kwargs):
        """Compatibility wrapper around StaticResult.plot_shear_diagram()."""
        if not hasattr(self, "_last_result"):
            raise RuntimeError(
                "plot_shear_diagram() is available only after solve()"
            )
        return self._last_result.plot_shear_diagram(**kwargs)

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
    displacement_dofs = ("ux", "uy")
    force_dofs = ("fx", "fy")

    def __init__(self,name="LT Model 01"):
        Model.__init__(self,name=name,mtype="triangle")
        self.dof = 2 # 2 DOF for triangle element (per node)
        
    def add_force(self,node,force):
        values = self._validated_component_vector(force, self.force_dofs, "force")
        self._record_applied_forces(node, **values)
        
    def add_constraint(self,node,**constraint):
        values = self._validated_named_components(
            constraint, self.displacement_dofs, "constraint"
        )
        for variable, value in values.items():
            setattr(node, variable, value)
        if values:
            self._record_prescribed_displacements(node, **values)
        
    def plot_model(self, show_reactions=False):
        """
        Plot mesh geometry, applied loads, constraints, and optional reactions.
        """
        import matplotlib.pyplot as plt
        from matplotlib.patches import Polygon
        from matplotlib.collections import PatchCollection
        
        fig = plt.figure()
        ax = fig.add_subplot(111)

        patches = []
        for elm in self.elements:
            _x,_y = [],[]
            for nd in elm.nodes:
                _x.append(nd.x)
                _y.append(nd.y)
            polygon = Polygon(list(zip(_x,_y)))
            patches.append(polygon)

        for nd in self.nodes:
            applied = self.applied_load(nd)
            if applied["fx"] > 0: self._draw_xforce(ax,nd.x,nd.y,1)
            if applied["fx"] < 0: self._draw_xforce(ax,nd.x,nd.y,-1)
            if applied["fy"] > 0: self._draw_yforce(ax,nd.x,nd.y,1)
            if applied["fy"] < 0: self._draw_yforce(ax,nd.x,nd.y,-1)

            if show_reactions:
                reaction = self.reaction(nd)
                if reaction["fx"] > 0: self._draw_xforce(ax,nd.x,nd.y,1,reaction=True)
                if reaction["fx"] < 0: self._draw_xforce(ax,nd.x,nd.y,-1,reaction=True)
                if reaction["fy"] > 0: self._draw_yforce(ax,nd.x,nd.y,1,reaction=True)
                if reaction["fy"] < 0: self._draw_yforce(ax,nd.x,nd.y,-1,reaction=True)

            if nd.ux == 0 and nd.uy == 0:
                self._draw_xyconstraint(ax,nd.x,nd.y)

        pc = PatchCollection(patches, color="#7CE7FF", edgecolor="k", alpha=0.4)
        ax.add_collection(pc)
        x0,x1,y0,y1 = self._rect_region()
        ax.set_xlim(x0,x1)
        ax.set_ylim(y0,y1)
        ax.set_title("Model %s"%(self.name))
        ax.set_aspect("equal")

    def _draw_xforce(self,axes,x,y,ddir=1,reaction=False):
        """
        Draw horizontal applied-load or reaction arrow.
        """
        dx, dy = self._calculate_arrow_size(), 0
        HW = dx/5.0
        HL = dx/3.0
        color = 'b' if reaction else 'r'
        arrow_props = dict(head_width=HW, head_length=HL, fc=color, ec=color)
        axes.arrow(x, y, ddir*dx, dy, **arrow_props)
        
    def _draw_yforce(self,axes,x,y,ddir=1,reaction=False):
        """
        Draw vertical applied-load or reaction arrow.
        """
        dx,dy = 0, self._calculate_arrow_size()
        HW = dy/5.0
        HL = dy/3.0
        color = 'b' if reaction else 'r'
        arrow_props = dict(head_width=HW, head_length=HL, fc=color, ec=color)
        axes.arrow(x, y, dx, ddir*dy, **arrow_props)
        
    def _draw_xyconstraint(self,axes,x,y):
        axes.plot(x, y, "gv", markersize=10, alpha=0.6)
        axes.plot(x, y, "g<", markersize=10, alpha=0.6)
        
    def _calculate_arrow_size(self):
        x0,x1,y0,y1 = self._rect_region(factor=10)
        sf = 8e-2
        kfx = sf*(x1-x0)
        kfy = sf*(y1-y0)
        return np.mean([kfx,kfy])
        
    def plot_nodal_result(self, var="ux"):
        """Compatibility wrapper for result-based nodal-field plotting."""
        if not hasattr(self, "_last_result"):
            raise RuntimeError(
                "plot_nodal_result() is available only after solve()"
            )
        return self._last_result.plot_nodal_field(var)

    def plot_element_result(self, var="sxx"):
        """Compatibility wrapper for result-based element-field plotting."""
        if not hasattr(self, "_last_result"):
            raise RuntimeError(
                "plot_element_result() is available only after solve()"
            )
        return self._last_result.plot_element_field(var)

    def show(self):
        """
        Show matplotlib plots
        """
        import matplotlib.pyplot as plt
        plt.show()
    
    def _rect_region(self,factor=7.0):
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
