"""Regression tests for public API integrity and model invariants."""

import pytest

from nusa import (
    Bar, BarModel, Beam, BeamModel, LinearTriangle, LinearTriangleModel,
    Model, Node, Spring, SpringModel, Truss, TrussModel,
)
from nusa import Element

from nusa import Material, Section
from nusa import Material, Section

def _make_beam(nodes, E, I):
    return Beam(nodes, material=Material(E=E), section=Section(I=I))



def _make_bar(nodes, E, A):
    return Bar(nodes, material=Material(E=E), section=Section(A=A))

def _make_truss(nodes, E, A):
    return Truss(nodes, material=Material(E=E), section=Section(A=A))



class MockBarElement(Element):
    def __init__(self, nodes):
        super().__init__("bar")
        self.nodes = tuple(nodes)


def test_add_element_rejects_nodes_outside_model():
    model = Model("membership", "bar")
    n1, n2 = Node((0.0, 0.0)), Node((1.0, 0.0))
    model.add_node(n1)
    element = MockBarElement((n1, n2))

    with pytest.raises(ValueError, match="do not belong"):
        model.add_element(element)

    assert model.n_elements == 0


def test_element_labels_and_duplicate_object_rules():
    model = Model("labels", "bar")
    nodes = [Node((float(i), 0.0)) for i in range(4)]
    model.add_nodes(nodes)

    e0 = MockBarElement((nodes[0], nodes[1]))
    e2 = MockBarElement((nodes[1], nodes[2]))
    auto = MockBarElement((nodes[2], nodes[3]))
    e0.label = 0
    e2.label = 2

    model.add_elements([e0, e2, auto])
    assert [element.label for element in model.elements] == [0, 2, 1]

    with pytest.raises(ValueError, match="already belongs"):
        model.add_element(e0)


def _result(builder):
    model = builder()
    return model.solve()


def _spring():
    model = SpringModel("spring report")
    n1, n2 = Node((0, 0)), Node((1, 0))
    model.add_nodes([n1, n2])
    model.add_element(Spring((n1, n2), 100.0))
    model.add_constraint(n1, ux=0.0)
    model.add_force(n2, (10.0,))
    return model


def _bar():
    model = BarModel("bar report")
    n1, n2 = Node((0, 0)), Node((1, 0))
    model.add_nodes([n1, n2])
    model.add_element(_make_bar((n1, n2), E=100.0, A=2.0))
    model.add_constraint(n1, ux=0.0)
    model.add_force(n2, (10.0,))
    return model


def _truss():
    model = TrussModel("truss report")
    n1, n2 = Node((0, 0)), Node((1, 0))
    model.add_nodes([n1, n2])
    model.add_element(_make_truss((n1, n2), E=100.0, A=1.0))
    model.add_constraint(n1, ux=0.0, uy=0.0)
    model.add_constraint(n2, uy=0.0)
    model.add_force(n2, (10.0, 0.0))
    return model


def _beam():
    model = BeamModel("beam report")
    n1, n2 = Node((0, 0)), Node((1, 0))
    model.add_nodes([n1, n2])
    model.add_element(_make_beam((n1, n2), E=1.0, I=1.0))
    model.add_constraint(n1, uy=0.0, ur=0.0)
    model.add_force(n2, (-1.0,))
    return model


def _triangle():
    model = LinearTriangleModel("triangle report")
    n1, n2, n3 = Node((0, 0)), Node((1, 0.5)), Node((0, 1))
    model.add_nodes([n1, n2, n3])
    model.add_element(LinearTriangle((n1, n2, n3), E=1000.0, nu=0.25, t=0.5))
    model.add_constraint(n1, ux=0.0, uy=0.0)
    model.add_constraint(n3, ux=0.0, uy=0.0)
    model.add_force(n2, (10.0, 0.0))
    return model


@pytest.mark.parametrize(
    ("builder", "columns"),
    [
        (_spring, ("FORCE I", "FORCE J")),
        (_bar, ("FORCE I", "FORCE J", "AXIAL FORCE", "AXIAL STRESS")),
        (_truss, ("AXIAL FORCE", "AXIAL STRESS")),
        (_beam, ("SHEAR FORCE I", "SHEAR FORCE J", "BENDING MOMENT I", "BENDING MOMENT J")),
        (_triangle, ("STRESS XX", "STRESS YY", "STRESS XY", "STRAIN XX", "STRAIN YY", "STRAIN XY")),
    ],
)
def test_simple_report_is_owned_by_static_result(builder, columns):
    result = _result(builder)
    report = result.simple_report(report_type="string")

    for section in (
        "NODAL DISPLACEMENTS", "APPLIED LOADS", "NODAL FORCES (K @ U)",
        "REACTIONS", "ELEMENT RESULTS", "FINITE ELEMENT MODEL INFO",
    ):
        assert section in report
    for column in columns:
        assert column in report


def test_spring_has_no_element_level_global_assembly_method():
    element = Spring((Node((0, 0)), Node((1, 0))), 100.0)
    assert not hasattr(element, "get_global_stiffness")
