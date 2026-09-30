"""Tests for geometry/topology core objects."""

import numpy as np
import pytest

from nusa import Element, Model, Node


class TestNode:
    def test_node_creation_and_label(self):
        node = Node((1.0, 2.0))
        assert node.x == 1.0
        assert node.y == 2.0
        np.testing.assert_allclose(node.coordinates, [1.0, 2.0])
        assert node.label is None

        node.label = "A"
        assert node.label == "A"

    def test_node_rejects_invalid_coordinates(self):
        with pytest.raises(ValueError, match="exactly two"):
            Node((1.0,))
        with pytest.raises(ValueError, match="finite"):
            Node((np.nan, 0.0))

    def test_node_has_no_solved_state(self):
        node = Node((0.0, 0.0))
        for name in (
            "ux", "uy", "ur", "fx", "fy", "m",
            "sx", "sy", "sxy", "seqv", "ex", "ey", "exy",
            "_elements",
        ):
            assert not hasattr(node, name)

    def test_node_repr_is_geometry_only(self):
        node = Node((1.0, 2.0))
        node.label = 1
        assert "Node 1" in str(node)
        assert "(1.0,2.0)" in repr(node)


class TestElement:
    def test_element_creation_and_label(self):
        element = Element("test_type")
        assert element.etype == "test_type"
        assert element.label is None

        element.label = 3
        assert element.label == 3

    def test_element_base_has_no_solved_force_state(self):
        element = Element("test")
        for name in ("fx", "fy", "sx", "sy", "sxy"):
            assert not hasattr(element, name)


class TestModel:
    def test_model_creation(self):
        model = Model("test_model", "bar")
        assert model.name == "test_model"
        assert model.mtype == "bar"
        assert model.n_nodes == 0
        assert model.n_elements == 0

    def test_add_node_and_labels(self):
        model = Model("test", "bar")
        n1 = Node((0, 0))
        n2 = Node((1, 0))
        n3 = Node((2, 0))
        n1.label = 0
        n2.label = 2

        model.add_nodes([n1, n2, n3])

        assert [node.label for node in model.nodes] == [0, 2, 1]
        assert [model._get_node_index(node) for node in model.nodes] == [0, 1, 2]

    def test_duplicate_node_label_is_rejected(self):
        model = Model("test", "bar")
        n1 = Node((0, 0))
        n2 = Node((1, 0))
        n1.label = n2.label = "A"
        model.add_node(n1)

        with pytest.raises(ValueError, match="already exists"):
            model.add_node(n2)

    def test_node_index_is_independent_of_label_mutation(self):
        model = Model("test", "bar")
        node = Node((0, 0))
        node.label = "support"
        model.add_node(node)
        node.label = "renamed"

        assert model._get_node_index(node) == 0

    def test_add_element_validates_type_membership_and_connectivity(self):
        class MockElement(Element):
            def __init__(self, nodes, etype="bar"):
                super().__init__(etype)
                self.nodes = tuple(nodes)

        model = Model("test", "bar")
        n1, n2 = Node((0, 0)), Node((1, 0))
        model.add_nodes([n1, n2])

        element = MockElement((n1, n2))
        model.add_element(element)
        assert model.elements == [element]
        assert element.label == 0

        with pytest.raises(ValueError, match="incompatible"):
            model.add_element(MockElement((n1, n2), etype="truss"))

        foreign = Node((2, 0))
        with pytest.raises(ValueError, match="do not belong"):
            model.add_element(MockElement((n1, foreign)))

    def test_model_repr(self):
        model = Model("test_model", "bar")
        model.add_nodes([Node((0, 0)), Node((1, 0))])
        text = str(model)
        assert "Model: test_model" in text
        assert "Nodes: 2" in text
        assert "Elements: 0" in text
        assert repr(model) == text
