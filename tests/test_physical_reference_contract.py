"""Physical property reference and sharing contracts for NuSA 0.5."""

import numpy as np
import pytest

from nusa import Bar, Beam, LinearTriangle, Material, Node, Section, Truss


@pytest.fixture
def line_nodes():
    return Node((0.0, 0.0)), Node((2.0, 0.0))


@pytest.mark.parametrize(
    ("element_type", "section"),
    [
        (Bar, Section(A=2.0)),
        (Truss, Section(A=2.0)),
        (Beam, Section(I=3.0)),
    ],
)
def test_structural_element_references_are_read_only(line_nodes, element_type, section):
    material = Material(E=100.0)
    element = element_type(line_nodes, material=material, section=section)
    assert element.material is material
    assert element.section is section
    for name, replacement in (("material", Material(E=200.0)), ("section", section)):
        with pytest.raises(AttributeError):
            setattr(element, name, replacement)
    assert element.material is material
    assert element.section is section


def test_same_material_and_section_can_be_shared_across_axial_elements(line_nodes):
    material = Material(E=120.0, nu=0.25)
    section = Section(A=2.0, I=3.0)
    bar = Bar(line_nodes, material=material, section=section)
    truss = Truss(line_nodes, material=material, section=section)
    beam = Beam(line_nodes, material=material, section=section)

    assert all(element.material is material for element in (bar, truss, beam))
    assert all(element.section is section for element in (bar, truss, beam))
    assert bar.E == truss.E == beam.E == 120.0
    assert bar.A == truss.A == 2.0
    assert beam.I == 3.0
    np.testing.assert_allclose(
        bar.get_element_stiffness(),
        (120.0 * 2.0 / 2.0) * np.array([[1.0, -1.0], [-1.0, 1.0]]),
    )
    assert truss.get_element_stiffness().shape == (4, 4)
    assert beam.get_element_stiffness().shape == (4, 4)


def test_material_can_be_reused_by_triangle_and_beam(line_nodes):
    material = Material(E=210e9, nu=0.3)
    beam = Beam(line_nodes, material=material, section=Section(I=8e-6))
    n1, n2 = line_nodes
    n3 = Node((0.0, 1.0))
    triangle = LinearTriangle((n1, n2, n3), material=material, thickness=0.02)
    assert beam.material is triangle.material is material
    assert triangle.nu == 0.3
    assert triangle.E == beam.E
