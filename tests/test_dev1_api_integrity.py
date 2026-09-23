"""Regression tests for dev1 public API integrity and model invariants."""

import pytest

from nusa.core import Element, Model, Node
from nusa.element import Bar, Beam, LinearTriangle, Spring, Truss
from nusa.model import BarModel, BeamModel, LinearTriangleModel, SpringModel, TrussModel


class MockBarElement(Element):
    def __init__(self, nodes):
        super().__init__("bar")
        self.nodes = nodes


def test_add_element_rejects_nodes_outside_model():
    model = Model("membership", "bar")
    n1 = Node((0.0, 0.0))
    n2 = Node((1.0, 0.0))
    model.add_node(n1)

    element = MockBarElement((n1, n2))

    with pytest.raises(ValueError, match="do not belong"):
        model.add_element(element)

    assert model.n_elements == 0
    assert element not in n1._elements


def test_add_element_rejects_duplicate_explicit_label():
    model = Model("labels", "bar")
    n1 = Node((0.0, 0.0))
    n2 = Node((1.0, 0.0))
    n3 = Node((2.0, 0.0))
    model.add_nodes([n1, n2, n3])

    e1 = MockBarElement((n1, n2))
    e2 = MockBarElement((n2, n3))
    e1.label = "member"
    e2.label = "member"

    model.add_element(e1)

    with pytest.raises(ValueError, match="already exists"):
        model.add_element(e2)

    assert model.elements == [e1]


def test_add_element_auto_label_uses_first_available_integer():
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


def test_add_element_rejects_same_object_twice():
    model = Model("duplicate object", "bar")
    n1 = Node((0.0, 0.0))
    n2 = Node((1.0, 0.0))
    model.add_nodes([n1, n2])
    element = MockBarElement((n1, n2))

    model.add_element(element)

    with pytest.raises(ValueError, match="already belongs"):
        model.add_element(element)

    assert model.n_elements == 1
    assert n1._elements.count(element) == 1
    assert n2._elements.count(element) == 1


def _solved_spring():
    model = SpringModel("spring report")
    n1 = Node((0.0, 0.0))
    n2 = Node((1.0, 0.0))
    model.add_nodes([n1, n2])
    model.add_element(Spring((n1, n2), 100.0))
    model.add_constraint(n1, ux=0.0)
    model.add_force(n2, (10.0,))
    model.solve()
    return model


def _solved_bar():
    model = BarModel("bar report")
    n1 = Node((0.0, 0.0))
    n2 = Node((1.0, 0.0))
    model.add_nodes([n1, n2])
    model.add_element(Bar((n1, n2), E=100.0, A=2.0))
    model.add_constraint(n1, ux=0.0)
    model.add_force(n2, (10.0,))
    model.solve()
    return model


def _solved_truss():
    model = TrussModel("truss report")
    n1 = Node((0.0, 0.0))
    n2 = Node((1.0, 0.0))
    model.add_nodes([n1, n2])
    model.add_element(Truss((n1, n2), E=100.0, A=1.0))
    model.add_constraint(n1, ux=0.0, uy=0.0)
    model.add_constraint(n2, uy=0.0)
    model.add_force(n2, (10.0, 0.0))
    model.solve()
    return model


def _solved_beam():
    model = BeamModel("beam report")
    n1 = Node((0.0, 0.0))
    n2 = Node((1.0, 0.0))
    model.add_nodes([n1, n2])
    model.add_element(Beam((n1, n2), E=1.0, I=1.0))
    model.add_constraint(n1, uy=0.0, ur=0.0)
    model.add_force(n2, (-1.0,))
    model.solve()
    return model


def _solved_triangle():
    model = LinearTriangleModel("triangle report")
    n1 = Node((0.0, 0.0))
    n2 = Node((1.0, 0.5))
    n3 = Node((0.0, 1.0))
    model.add_nodes([n1, n2, n3])
    model.add_element(LinearTriangle((n1, n2, n3), E=1000.0, nu=0.25, t=0.5))
    model.add_constraint(n1, ux=0.0, uy=0.0)
    model.add_constraint(n3, ux=0.0, uy=0.0)
    model.add_force(n2, (10.0, 0.0))
    model.solve()
    return model


@pytest.mark.parametrize(
    ("builder", "expected_columns"),
    [
        (_solved_spring, ("FORCE I", "FORCE J")),
        (
            _solved_bar,
            ("FORCE I", "FORCE J", "AXIAL FORCE", "AXIAL STRESS"),
        ),
        (_solved_truss, ("AXIAL FORCE", "AXIAL STRESS")),
        (
            _solved_beam,
            (
                "SHEAR FORCE I",
                "SHEAR FORCE J",
                "BENDING MOMENT I",
                "BENDING MOMENT J",
            ),
        ),
        (
            _solved_triangle,
            (
                "STRESS XX",
                "STRESS YY",
                "STRESS XY",
                "STRAIN XX",
                "STRAIN YY",
                "STRAIN XY",
            ),
        ),
    ],
)
def test_simple_report_is_available_for_all_public_models(builder, expected_columns):
    model = builder()
    report = model.simple_report(report_type="string")

    assert "NODAL DISPLACEMENTS" in report
    assert "APPLIED LOADS" in report
    assert "NODAL FORCES (K @ U)" in report
    assert "REACTIONS" in report
    assert "ELEMENT RESULTS" in report
    assert "FINITE ELEMENT MODEL INFO" in report
    for column in expected_columns:
        assert column in report


def test_simple_report_rejects_unknown_report_type():
    model = _solved_spring()

    with pytest.raises(ValueError, match="report_type"):
        model.simple_report(report_type="unknown")


def test_simple_report_requires_solved_model():
    model = SpringModel("unsolved")

    with pytest.raises(RuntimeError, match="after solve"):
        model.simple_report(report_type="string")


def test_simple_report_write_mode(tmp_path):
    model = _solved_truss()
    report_path = tmp_path / "report.txt"

    result = model.simple_report(report_type="write", fname=report_path)

    assert result is None
    text = report_path.read_text(encoding="utf-8")
    assert "truss report" in text
    assert "ELEMENT RESULTS" in text


def test_spring_has_no_element_level_global_assembly_method():
    n1 = Node((0.0, 0.0))
    n2 = Node((1.0, 0.0))
    element = Spring((n1, n2), 100.0)

    assert not hasattr(element, "get_global_stiffness")
