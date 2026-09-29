"""Regression tests for reporting and problem-plot force semantics."""

import matplotlib.pyplot as plt

from nusa import LinearTriangle, LinearTriangleModel, Node, Spring, SpringModel, Truss, TrussModel


def test_spring_result_report_separates_force_quantities():
    model=SpringModel("Report semantics")
    n1,n2=Node((0,0)),Node((0,0))
    model.add_nodes([n1,n2]); model.add_element(Spring((n1,n2),100)); model.add_constraint(n1,ux=0); model.add_force(n2,(50,))
    report=model.solve().simple_report(report_type="string")
    assert "APPLIED LOADS" in report
    assert "NODAL FORCES (K @ U)" in report
    assert "REACTIONS" in report
    assert "50" in report and "-50" in report


def test_truss_problem_plot_shows_applied_loads_only(monkeypatch):
    model=TrussModel("Plot semantics")
    n1,n2=Node((0,0)),Node((1,0))
    model.add_nodes([n1,n2]); model.add_element(Truss((n1,n2),100,1)); model.add_constraint(n1,ux=0,uy=0); model.add_constraint(n2,uy=0); model.add_force(n2,(10,0))

    arrows=[]
    monkeypatch.setattr(model,"_draw_xforce",lambda axes,x,y,ddir=1,reaction=False: arrows.append((x,y,ddir,reaction)))
    monkeypatch.setattr(model,"_draw_yforce",lambda *args,**kwargs:None)
    monkeypatch.setattr(model,"_draw_xconstraint",lambda *args,**kwargs:None)
    monkeypatch.setattr(model,"_draw_yconstraint",lambda *args,**kwargs:None)

    model.plot_model()
    assert arrows==[(1.0,0.0,1,False)]
    plt.close("all")


def test_triangle_problem_plot_preserves_negative_load_direction(monkeypatch):
    model=LinearTriangleModel("Negative load plot")
    n1,n2,n3=Node((0,0)),Node((1,0)),Node((0,1))
    model.add_nodes([n1,n2,n3]); model.add_element(LinearTriangle((n1,n2,n3),1000,.25,.5)); model.add_force(n2,(-10,-5))
    arrows=[]
    monkeypatch.setattr(model,"_draw_xforce",lambda axes,x,y,ddir=1,reaction=False: arrows.append(("x",x,y,ddir,reaction)))
    monkeypatch.setattr(model,"_draw_yforce",lambda axes,x,y,ddir=1,reaction=False: arrows.append(("y",x,y,ddir,reaction)))
    monkeypatch.setattr(model,"_draw_xyconstraint",lambda *args,**kwargs:None)
    model.plot_model()
    assert ("x",1.0,0.0,-1,False) in arrows
    assert ("y",1.0,0.0,-1,False) in arrows
    plt.close("all")
