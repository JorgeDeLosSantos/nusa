"""Regression tests for result-based reports and element result names."""

import pytest

from nusa import Bar, BarModel, Beam, BeamModel, LinearTriangle, LinearTriangleModel, Node, Spring, SpringModel, Truss, TrussModel

from nusa import Material, Section
from nusa import Material, Section

def _make_beam(nodes, E, I):
    return Beam(nodes, material=Material(E=E), section=Section(I=I))



def _make_bar(nodes, E, A):
    return Bar(nodes, material=Material(E=E), section=Section(A=A))

def _make_truss(nodes, E, A):
    return Truss(nodes, material=Material(E=E), section=Section(A=A))



def _result(kind):
    if kind=="spring":
        m=SpringModel("spring report"); a,b=Node((0,0)),Node((0,0)); e=Spring((a,b),100); m.add_nodes([a,b]); m.add_element(e); m.add_constraint(a,ux=0); m.add_force(b,(10,))
    elif kind=="bar":
        m=BarModel("bar report"); a,b=Node((0,0)),Node((2,0)); e=_make_bar((a,b),100,2); m.add_nodes([a,b]); m.add_element(e); m.add_constraint(a,ux=0); m.add_force(b,(10,))
    elif kind=="truss":
        m=TrussModel("truss report"); a,b=Node((0,0)),Node((2,0)); e=_make_truss((a,b),100,2); m.add_nodes([a,b]); m.add_element(e); m.add_constraint(a,ux=0,uy=0); m.add_constraint(b,uy=0); m.add_force(b,(10,0))
    elif kind=="beam":
        m=BeamModel("beam report"); a,b=Node((0,0)),Node((2,0)); e=_make_beam((a,b),100,1); m.add_nodes([a,b]); m.add_element(e); m.add_constraint(a,uy=0,ur=0); m.add_force(b,(-10,))
    else:
        m=LinearTriangleModel("triangle report"); a,b,c=Node((0,0)),Node((1,.5)),Node((0,1)); e=LinearTriangle((a,b,c),200e9,.3,.1); m.add_nodes([a,b,c]); m.add_element(e); m.add_constraint(a,ux=0,uy=0); m.add_constraint(c,ux=0,uy=0); m.add_force(b,(1000,0))
    return m.solve()


@pytest.mark.parametrize(("kind","headers"),[
    ("spring",("FORCE I","FORCE J")),
    ("bar",("FORCE I","FORCE J","AXIAL FORCE","AXIAL STRESS")),
    ("truss",("AXIAL FORCE","AXIAL STRESS")),
    ("beam",("SHEAR FORCE I","SHEAR FORCE J","BENDING MOMENT I","BENDING MOMENT J")),
    ("triangle",("STRESS XX","STRESS YY","STRESS XY","STRAIN XX","STRAIN YY","STRAIN XY")),
])
def test_element_report_uses_canonical_names(kind,headers):
    report=_result(kind).simple_report(report_type="string")
    for header in headers:
        assert header in report


def test_bar_report_exposes_physical_axial_force_and_stress():
    report=_result("bar").simple_report(report_type="string")
    assert "AXIAL FORCE" in report and "AXIAL STRESS" in report
    assert "10" in report and "5" in report
