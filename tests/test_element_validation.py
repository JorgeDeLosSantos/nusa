"""Validation tests for finite-element connectivity and physical properties."""

import numpy as np
import pytest

from nusa import Bar, Beam, LinearTriangle, Node, Spring, Truss


@pytest.fixture
def line_nodes():
    return Node((0.0, 0.0)), Node((1.0, 0.0))


@pytest.mark.parametrize(
    ("factory", "kwargs"),
    [
        (Spring, {"ke": 0.0}),
        (Spring, {"ke": -1.0}),
        (Bar, {"E": 0.0, "A": 1.0}),
        (Bar, {"E": 1.0, "A": 0.0}),
        (Truss, {"E": -1.0, "A": 1.0}),
        (Truss, {"E": 1.0, "A": -1.0}),
        (Beam, {"E": 0.0, "I": 1.0}),
        (Beam, {"E": 1.0, "I": 0.0}),
    ],
)
def test_two_node_elements_reject_nonpositive_properties(line_nodes, factory, kwargs):
    with pytest.raises(ValueError, match="positive"):
        factory(line_nodes, **kwargs)


@pytest.mark.parametrize("invalid_value", [np.nan, np.inf, -np.inf])
def test_elements_reject_nonfinite_properties(line_nodes, invalid_value):
    n1, n2 = line_nodes

    with pytest.raises(ValueError, match="finite positive"):
        Spring((n1, n2), invalid_value)

    with pytest.raises(ValueError, match="finite positive"):
        Bar((n1, n2), invalid_value, 1.0)

    with pytest.raises(ValueError, match="finite positive"):
        Truss((n1, n2), 1.0, invalid_value)

    with pytest.raises(ValueError, match="finite positive"):
        Beam((n1, n2), 1.0, invalid_value)


@pytest.mark.parametrize("element_type", [Bar, Truss, Beam])
def test_length_dependent_elements_reject_coincident_nodes(element_type):
    n1 = Node((0.0, 0.0))
    n2 = Node((0.0, 0.0))

    if element_type is Beam:
        kwargs = {"E": 1.0, "I": 1.0}
    else:
        kwargs = {"E": 1.0, "A": 1.0}

    with pytest.raises(ValueError, match="distinct node coordinates"):
        element_type((n1, n2), **kwargs)


def test_spring_allows_coincident_nodes():
    n1 = Node((0.0, 0.0))
    n2 = Node((0.0, 0.0))

    spring = Spring((n1, n2), 1000.0)

    np.testing.assert_allclose(
        spring.get_element_stiffness(),
        [[1000.0, -1000.0], [-1000.0, 1000.0]],
    )


@pytest.mark.parametrize("element_type", [Spring, Bar, Truss, Beam])
def test_two_node_elements_require_exact_connectivity_size(element_type):
    n1 = Node((0.0, 0.0))
    n2 = Node((1.0, 0.0))
    n3 = Node((2.0, 0.0))

    if element_type is Spring:
        args = (1000.0,)
    elif element_type is Beam:
        args = (1.0, 1.0)
    else:
        args = (1.0, 1.0)

    with pytest.raises(ValueError, match="exactly 2 nodes"):
        element_type((n1, n2, n3), *args)


def test_element_connectivity_requires_node_objects():
    n1 = Node((0.0, 0.0))

    with pytest.raises(ValueError, match="Node objects"):
        Bar((n1, (1.0, 0.0)), E=1.0, A=1.0)


def test_element_connectivity_requires_finite_coordinates():
    n1 = Node((0.0, 0.0))
    n2 = Node((np.nan, 1.0))

    with pytest.raises(ValueError, match="coordinates must be finite"):
        Truss((n1, n2), E=1.0, A=1.0)


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
