"""Physical and assembly invariant tests for the 0.4 architecture."""

import numpy as np

from nusa import Bar, BarModel, Beam, BeamModel, LinearTriangle, LinearTriangleModel, Node, Spring, SpringModel, Truss, TrussModel
from nusa.analysis import _assemble_stiffness

from nusa import Material, Section
from nusa import Material, Section
from nusa import Material

def _make_triangle(nodes, E, nu, t):
    return LinearTriangle(nodes, material=Material(E=E, nu=nu), thickness=t)



def _make_beam(nodes, E, I):
    return Beam(nodes, material=Material(E=E), section=Section(I=I))



def _make_bar(nodes, E, A):
    return Bar(nodes, material=Material(E=E), section=Section(A=A))

def _make_truss(nodes, E, A):
    return Truss(nodes, material=Material(E=E), section=Section(A=A))



def _assert_symmetric(model):
    K = _assemble_stiffness(model)
    np.testing.assert_allclose(K, K.T, rtol=0.0, atol=1e-12)


def test_global_stiffness_is_symmetric_for_all_public_models():
    spring=SpringModel(); a,b=Node((0,0)),Node((0,0)); spring.add_nodes([a,b]); spring.add_element(Spring((a,b),100)); _assert_symmetric(spring)
    bar=BarModel(); a,b=Node((0,0)),Node((2,0)); bar.add_nodes([a,b]); bar.add_element(_make_bar((a,b),200,3)); _assert_symmetric(bar)
    truss=TrussModel(); a,b=Node((0,0)),Node((3,4)); truss.add_nodes([a,b]); truss.add_element(_make_truss((a,b),200,2)); _assert_symmetric(truss)
    beam=BeamModel(); a,b=Node((0,0)),Node((2,0)); beam.add_nodes([a,b]); beam.add_element(_make_beam((a,b),200,4)); _assert_symmetric(beam)
    tri=LinearTriangleModel(); a,b,c=Node((0,0)),Node((1,0)),Node((0,1)); tri.add_nodes([a,b,c]); tri.add_element(_make_triangle((a,b,c),1000,.25,.5)); _assert_symmetric(tri)


def test_spring_and_bar_global_equilibrium():
    spring=SpringModel(); a,b=Node((0,0)),Node((0,0)); spring.add_nodes([a,b]); spring.add_element(Spring((a,b),100)); spring.add_constraint(a,ux=0); spring.add_force(b,(25,))
    sr=spring.solve()
    assert np.isclose(sr.applied_loads.sum()+sr.reactions.sum(),0)

    bar=BarModel(); a,b,c=Node((0,0)),Node((1,0)),Node((2,0)); bar.add_nodes([a,b,c]); bar.add_elements([_make_bar((a,b),100,1),_make_bar((b,c),100,1)]); bar.add_constraint(a,ux=0); bar.add_force(c,(40,))
    br=bar.solve()
    assert np.isclose(br.applied_loads.sum()+br.reactions.sum(),0)


def test_truss_triangle_and_beam_equilibrium():
    tr=TrussModel(); n1,n2,n3=Node((0,0)),Node((1,1)),Node((2,0)); tr.add_nodes([n1,n2,n3]); tr.add_elements([_make_truss((n1,n2),1000,1),_make_truss((n2,n3),1000,1),_make_truss((n1,n3),1000,1)]); tr.add_constraint(n1,ux=0,uy=0); tr.add_constraint(n3,ux=0,uy=0); tr.add_force(n2,(10,-30))
    r=tr.solve()
    np.testing.assert_allclose(r.applied_loads.reshape(-1,2).sum(0)+r.reactions.reshape(-1,2).sum(0),[0,0],atol=1e-10)

    tri=LinearTriangleModel(); a,b,c=Node((0,0)),Node((1,.5)),Node((0,1)); tri.add_nodes([a,b,c]); tri.add_element(_make_triangle((a,b,c),1000,.25,.5)); tri.add_constraint(a,ux=0,uy=0); tri.add_constraint(c,ux=0,uy=0); tri.add_force(b,(12,-7))
    rr=tri.solve()
    np.testing.assert_allclose(rr.applied_loads.reshape(-1,2).sum(0)+rr.reactions.reshape(-1,2).sum(0),[0,0],atol=1e-10)

    beam=BeamModel(); a,b,c=Node((0,0)),Node((3,0)),Node((5,0)); beam.add_nodes([a,b,c]); beam.add_elements([_make_beam((a,b),200,4),_make_beam((b,c),200,4)]); beam.add_constraint(a,uy=0,ur=0); beam.add_force(b,(-10,)); beam.add_moment(c,(6,))
    rb=beam.solve()
    applied_force=sum(rb.applied_load(n)["fy"] for n in beam.nodes)
    reaction_force=sum(rb.reaction(n)["fy"] for n in beam.nodes)
    assert np.isclose(applied_force+reaction_force,0,atol=1e-10)
