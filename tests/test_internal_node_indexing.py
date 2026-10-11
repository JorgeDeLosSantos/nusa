"""Regression tests for public labels and internal contiguous indexing."""

import numpy as np

from nusa import Bar, BarModel, Beam, BeamModel, LinearTriangle, LinearTriangleModel, Node, Spring, SpringModel, Truss, TrussModel

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



def test_spring_solver_accepts_string_labels_and_label_mutation():
    model = SpringModel("labels")
    n1, n2 = Node((0,0)), Node((0,0))
    n1.label, n2.label = "fixed", "tip"
    model.add_nodes([n1,n2])
    model.add_element(Spring((n1,n2),100))
    n1.label, n2.label = "support-A", "load-point"
    model.add_constraint(n1,ux=0)
    model.add_force(n2,(50,))

    result = model.solve()

    assert result.node_labels == ("support-A", "load-point")
    np.testing.assert_allclose(result.displacements, [0.0, 0.5])


def test_sparse_bar_labels_do_not_affect_solver_order():
    model = BarModel("sparse")
    n1,n2,n3 = Node((0,0)),Node((1,0)),Node((2,0))
    n1.label,n2.label,n3.label = 10,30,80
    model.add_nodes([n1,n2,n3])
    model.add_elements([_make_bar((n1,n2),1,1),_make_bar((n2,n3),1,1)])
    model.add_constraint(n1,ux=0)
    model.add_force(n3,(1,))

    result = model.solve()

    np.testing.assert_allclose(result.displacements,[0,1,2])
    assert result.node_labels == (10,30,80)


def test_truss_and_beam_accept_string_labels():
    truss = TrussModel("truss")
    a,b = Node((0,0)),Node((2,0)); a.label,b.label="A","B"
    e=_make_truss((a,b),100,2)
    truss.add_nodes([a,b]); truss.add_element(e); truss.add_constraint(a,ux=0,uy=0); truss.add_constraint(b,uy=0); truss.add_force(b,(10,0))
    tr = truss.solve()
    assert np.isclose(tr.displacement(b)["ux"],0.1)

    beam = BeamModel("beam")
    f,t = Node((0,0)),Node((1,0)); f.label,t.label="fixed","tip"
    beam.add_nodes([f,t]); beam.add_element(_make_beam((f,t),1,1)); beam.add_constraint(f,uy=0,ur=0); beam.add_force(t,(-1,))
    br=beam.solve()
    assert np.isclose(br.displacement(t)["uy"],-1/3)
    assert np.isclose(br.displacement(t)["ur"],-0.5)


def test_triangle_result_connectivity_uses_internal_indices_not_labels():
    model = LinearTriangleModel("CST")
    n1,n2,n3=Node((0,0)),Node((1,.5)),Node((0,1))
    n1.label,n2.label,n3.label="left-bottom","loaded","left-top"
    e=_make_triangle((n1,n2,n3),200e9,.3,.1)
    model.add_nodes([n1,n2,n3]); model.add_element(e)
    model.add_constraint(n1,ux=0,uy=0); model.add_constraint(n3,ux=0,uy=0); model.add_force(n2,(1000,0))

    result=model.solve()

    assert result.node_labels == ("left-bottom","loaded","left-top")
    assert result.connectivity == ((0,1,2),)
