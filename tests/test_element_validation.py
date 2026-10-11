"""Validation tests for finite-element connectivity and physical properties."""

import numpy as np
import pytest

from nusa import (
    Bar,
    Beam,
    LinearTriangle,
    Material,
    Node,
    Section,
    Spring,
    Truss,
)


@pytest.fixture
def line_nodes():
    return Node((0.0, 0.0)), Node((1.0, 0.0))


@pytest.mark.parametrize(
    ("factory", "kwargs"),
    [
        (Spring, {"ke": 0.0}),
        (Spring, {"ke": -1.0}),
    ],
)
def test_direct_element_properties_reject_nonpositive_values(
    line_nodes, factory, kwargs
):
    with pytest.raises(ValueError, match="positive"):
        factory(line_nodes, **kwargs)


@pytest.mark.parametrize("invalid_value", [np.nan, np.inf, -np.inf])
def test_spring_rejects_nonfinite_stiffness(line_nodes, invalid_value):
    with pytest.raises(ValueError, match="finite positive"):
        Spring(line_nodes, invalid_value)


def test_beam_requires_material_instance(line_nodes):
    with pytest.raises(TypeError, match="Material instance"):
        Beam(line_nodes, material=object(), section=Section(I=1.0))


def test_beam_requires_section_instance(line_nodes):
    with pytest.raises(TypeError, match="Section instance"):
        Beam(line_nodes, material=Material(E=1.0), section=object())


def test_beam_requires_section_second_moment(line_nodes):
    with pytest.raises(ValueError, match="requires section property 'I'"):
        Beam(
            line_nodes,
            material=Material(E=1.0),
            section=Section(A=1.0),
        )


def test_beam_properties_delegate_to_domain_objects(line_nodes):
    material = Material(E=200.0)
    section = Section(I=4.0)

    element = Beam(line_nodes, material=material, section=section)

    assert element.material is material
    assert element.section is section
    assert element.E == 200.0
    assert element.I == 4.0


def test_beam_material_and_section_are_keyword_only(line_nodes):
    material = Material(E=1.0)
    section = Section(I=1.0)

    with pytest.raises(TypeError):
        Beam(line_nodes, material, section)

@pytest.mark.parametrize("element_type", [Bar, Truss])
def test_axial_elements_require_material_instance(line_nodes, element_type):
    section = Section(A=1.0)

    with pytest.raises(TypeError, match="Material instance"):
        element_type(line_nodes, material=object(), section=section)


@pytest.mark.parametrize("element_type", [Bar, Truss])
def test_axial_elements_require_section_instance(line_nodes, element_type):
    material = Material(E=1.0)

    with pytest.raises(TypeError, match="Section instance"):
        element_type(line_nodes, material=material, section=object())


@pytest.mark.parametrize("element_type", [Bar, Truss])
def test_axial_elements_require_section_area(line_nodes, element_type):
    material = Material(E=1.0)
    section = Section(I=1.0)

    with pytest.raises(ValueError, match="requires section property 'A'"):
        element_type(line_nodes, material=material, section=section)


@pytest.mark.parametrize("element_type", [Bar, Truss])
def test_axial_element_properties_delegate_to_domain_objects(line_nodes, element_type):
    material = Material(E=200.0)
    section = Section(A=3.0)

    element = element_type(line_nodes, material=material, section=section)

    assert element.material is material
    assert element.section is section
    assert element.E == 200.0
    assert element.A == 3.0


@pytest.mark.parametrize("element_type", [Bar, Truss])
def test_axial_element_material_and_section_are_keyword_only(line_nodes, element_type):
    material = Material(E=1.0)
    section = Section(A=1.0)

    with pytest.raises(TypeError):
        element_type(line_nodes, material, section)


@pytest.mark.parametrize("element_type", [Bar, Truss])
def test_axial_elements_reject_coincident_nodes(element_type):
    n1 = Node((0.0, 0.0))
    n2 = Node((0.0, 0.0))

    with pytest.raises(ValueError, match="distinct node coordinates"):
        element_type(
            (n1, n2),
            material=Material(E=1.0),
            section=Section(A=1.0),
        )


def test_beam_rejects_coincident_nodes():
    n1 = Node((0.0, 0.0))
    n2 = Node((0.0, 0.0))

    with pytest.raises(ValueError, match="distinct node coordinates"):
        Beam(
            (n1, n2),
            material=Material(E=1.0),
            section=Section(I=1.0),
        )


def test_spring_allows_coincident_nodes():
    n1 = Node((0.0, 0.0))
    n2 = Node((0.0, 0.0))

    spring = Spring((n1, n2), 1000.0)

    np.testing.assert_allclose(
        spring.get_element_stiffness(),
        [[1000.0, -1000.0], [-1000.0, 1000.0]],
    )


@pytest.mark.parametrize("element_type", [Bar, Truss])
def test_axial_elements_require_exact_connectivity_size(element_type):
    nodes = (
        Node((0.0, 0.0)),
        Node((1.0, 0.0)),
        Node((2.0, 0.0)),
    )

    with pytest.raises(ValueError, match="exactly 2 nodes"):
        element_type(
            nodes,
            material=Material(E=1.0),
            section=Section(A=1.0),
        )


def test_spring_requires_exact_connectivity_size():
    nodes = (
        Node((0.0, 0.0)),
        Node((1.0, 0.0)),
        Node((2.0, 0.0)),
    )

    with pytest.raises(ValueError, match="exactly 2 nodes"):
        Spring(nodes, 1000.0)


def test_beam_requires_exact_connectivity_size():
    nodes = (
        Node((0.0, 0.0)),
        Node((1.0, 0.0)),
        Node((2.0, 0.0)),
    )

    with pytest.raises(ValueError, match="exactly 2 nodes"):
        Beam(
            nodes,
            material=Material(E=1.0),
            section=Section(I=1.0),
        )


def test_element_connectivity_requires_node_objects():
    n1 = Node((0.0, 0.0))

    with pytest.raises(ValueError, match="Node objects"):
        Bar(
            (n1, (1.0, 0.0)),
            material=Material(E=1.0),
            section=Section(A=1.0),
        )


def test_element_connectivity_requires_finite_coordinates():
    n1 = Node((0.0, 0.0))
    n2 = Node((1.0, 1.0))
    n2.coordinates[0] = np.nan

    with pytest.raises(ValueError, match="coordinates must be finite"):
        Truss(
            (n1, n2),
            material=Material(E=1.0),
            section=Section(A=1.0),
        )


@pytest.mark.parametrize("nu", [-1.0, 0.5, -1.1, 0.75, np.nan, np.inf])
def test_linear_triangle_rejects_invalid_poisson_ratio(nu):
    nodes = (
        Node((0.0, 0.0)),
        Node((1.0, 0.0)),
        Node((0.0, 1.0)),
    )

    with pytest.raises(ValueError, match="range -1 < nu < 0.5"):
        LinearTriangle(nodes, E=1.0, nu=nu, t=1.0)


def test_linear_triangle_accepts_auxetic_poisson_ratio():
    nodes = (
        Node((0.0, 0.0)),
        Node((1.0, 0.0)),
        Node((0.0, 1.0)),
    )

    element = LinearTriangle(nodes, E=1.0, nu=-0.2, t=1.0)

    assert element.nu == -0.2


@pytest.mark.parametrize(
    ("E", "t"),
    [
        (0.0, 1.0),
        (-1.0, 1.0),
        (1.0, 0.0),
        (1.0, -1.0),
        (np.inf, 1.0),
        (1.0, np.nan),
    ],
)
def test_linear_triangle_rejects_invalid_material_or_thickness(E, t):
    nodes = (
        Node((0.0, 0.0)),
        Node((1.0, 0.0)),
        Node((0.0, 1.0)),
    )

    with pytest.raises(ValueError, match="positive"):
        LinearTriangle(nodes, E=E, nu=0.3, t=t)


def test_linear_triangle_requires_exactly_three_nodes():
    nodes = (Node((0.0, 0.0)), Node((1.0, 0.0)))

    with pytest.raises(ValueError, match="exactly 3 nodes"):
        LinearTriangle(nodes, E=1.0, nu=0.3, t=1.0)
