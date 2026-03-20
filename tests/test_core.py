"""
Tests for core.py module
"""
import numpy as np
import pytest
from nusa.core import Model, Element, Node


class TestNode:
    """Test Node class functionality"""
    
    def test_node_creation(self):
        """Test basic node creation with coordinates"""
        node = Node((1.0, 2.0))
        assert node.x == 1.0
        assert node.y == 2.0
        assert np.array_equal(node.coordinates, np.array([1.0, 2.0]))
    
    def test_node_label(self):
        """Test node label setting and getting"""
        node = Node((0, 0))
        assert node.label is None
        
        node.label = 5
        assert node.label == 5
    
    def test_node_displacements(self):
        """Test node displacement properties"""
        node = Node((0, 0))
        
        # Test initial values
        assert np.isnan(node.ux)
        assert np.isnan(node.uy)
        assert np.isnan(node.ur)
        
        # Test setting displacements
        node.ux = 0.01
        node.uy = 0.02
        node.ur = 0.005
        
        assert node.ux == 0.01
        assert node.uy == 0.02
        assert node.ur == 0.005
    
    def test_node_forces(self):
        """Test node force properties"""
        node = Node((0, 0))
        
        # Test initial values
        assert node.fx == 0.0
        assert node.fy == 0.0
        assert node.m == 0.0
        
        # Test setting forces
        node.fx = 100.0
        node.fy = 200.0
        node.m = 50.0
        
        assert node.fx == 100.0
        assert node.fy == 200.0
        assert node.m == 50.0
    
    def test_node_stresses(self):
        """Test node stress properties"""
        node = Node((0, 0))
        
        # Test initial values
        assert node.sx == 0.0
        assert node.sy == 0.0
        assert node.sxy == 0.0
        assert node.seqv == 0.0
    
    def test_node_str_repr(self):
        """Test string representation of node"""
        node = Node((1.0, 2.0))
        node.label = 1
        
        repr_str = repr(node)
        assert "Node 1" in repr_str
        assert "(1.0,2.0)" in repr_str
    
    def test_node_elements(self):
        """Test node element association"""
        node = Node((0, 0))
        assert len(node._elements) == 0
        
        # Add a mock element
        class MockElement:
            def __init__(self):
                self.sx = 100.0
                self.sy = 50.0
                self.sxy = 25.0
        
        mock_element = MockElement()
        node.add_element(mock_element)
        
        assert len(node._elements) == 1
        assert node._elements[0] == mock_element


class TestElement:
    """Test Element class functionality"""
    
    def test_element_creation(self):
        """Test basic element creation"""
        element = Element("test_type")
        assert element.etype == "test_type"
        assert element.label is None
    
    def test_element_label(self):
        """Test element label setting"""
        element = Element("test")
        element.set_label(1)
        assert element.label == 1
    
    def test_element_forces(self):
        """Test element force properties"""
        element = Element("test")
        
        # Test initial values
        assert element.fx == 0.0
        assert element.fy == 0.0
        
        # Test setting forces
        element.fx = 150.0
        element.fy = 75.0
        
        assert element.fx == 150.0
        assert element.fy == 75.0
        
        # Test get_element_forces method
        forces = element.get_element_forces()
        assert forces == (150.0, 75.0)
    
    def test_element_str(self):
        """Test string representation of element"""
        element = Element("test_type")
        str_repr = str(element)
        assert "Element" in str_repr


class TestModel:
    """Test Model class functionality"""
    
    def test_model_creation(self):
        """Test basic model creation"""
        model = Model("test_model", "bar")
        assert model.name == "test_model"
        assert model.mtype == "bar"
        assert model.n_nodes == 0
        assert model.n_elements == 0
    
    def test_add_node(self):
        """Test adding nodes to model"""
        model = Model("test", "bar")
        node = Node((0, 0))
        
        model.add_node(node)
        assert model.n_nodes == 1
        assert node in model.nodes
        assert node.label == 0  # Should be auto-assigned
    
    def test_add_nodes(self):
        """Test adding multiple nodes to model"""
        model = Model("test", "bar")
        nodes = [Node((i, i)) for i in range(3)]
        
        model.add_nodes(nodes)
        assert model.n_nodes == 3
        assert all(node in model.nodes for node in nodes)
    
    def test_add_element(self):
        """Test adding elements to model"""
        model = Model("test", "bar")
        
        # Create nodes
        n1 = Node((0, 0))
        n2 = Node((1, 0))
        
        # Create element
        class MockElement(Element):
            def __init__(self, nodes):
                super().__init__("bar")
                self.nodes = nodes
        
        element = MockElement([n1, n2])
        
        model.add_element(element)
        assert model.n_elements == 1
        assert element in model.elements
        assert element.label == 0  # Should be auto-assigned
    
    def test_add_element_wrong_type(self):
        """Test adding element with wrong type raises error"""
        model = Model("test", "bar")
        
        class MockElement(Element):
            def __init__(self):
                super().__init__("truss")  # Wrong type
        
        element = MockElement()
        
        with pytest.raises(ValueError, match="Element type 'truss' incompatible"):
            model.add_element(element)
    
    def test_add_elements(self):
        """Test adding multiple elements to model"""
        model = Model("test", "bar")
        
        class MockElement(Element):
            def __init__(self, label):
                super().__init__("bar")
                self.label = label
                self.nodes = []
        
        elements = [MockElement(i) for i in range(3)]
        
        model.add_elements(elements)
        assert model.n_elements == 3
        assert all(element in model.elements for element in elements)
    
    def test_model_str_repr(self):
        """Test string representation of model"""
        model = Model("test_model", "bar")
        n1 = Node((0, 0))
        n2 = Node((1, 0))
        model.add_node(n1)
        model.add_node(n2)
        
        str_repr = str(model)
        repr_repr = repr(model)
        
        assert "Model: test_model" in str_repr
        assert "Nodes: 2" in str_repr
        assert "Elements: 0" in str_repr
        assert str_repr == repr_repr
    
    def test_model_properties(self):
        """Test model property access"""
        model = Model("test", "bar")
        
        # Test empty model
        assert model.nodes == []
        assert model.elements == []
        
        # Add nodes and elements
        n1 = Node((0, 0))
        n2 = Node((1, 0))
        model.add_node(n1)
        model.add_node(n2)
        
        class MockElement(Element):
            def __init__(self):
                super().__init__("bar")
                self.nodes = [n1, n2]
        
        element = MockElement()
        model.add_element(element)
        
        # Test properties
        assert len(model.nodes) == 2
        assert len(model.elements) == 1
        assert model.n_nodes == 2
        assert model.n_elements == 1


class TestModelIntegration:
    """Integration tests for Model with real nodes and elements"""
    
    def test_full_model_workflow(self):
        """Test complete workflow of creating and populating a model"""
        # Create model
        model = Model("integration_test", "bar")
        
        # Create and add nodes
        nodes = [
            Node((0.0, 0.0)),
            Node((1.0, 0.0)),
            Node((2.0, 0.0))
        ]
        model.add_nodes(nodes)
        
        # Verify nodes
        assert model.n_nodes == 3
        for i, node in enumerate(nodes):
            assert node.label == i
            assert node in model.nodes
        
        # Create and add elements
        class BarElement(Element):
            def __init__(self, nodes):
                super().__init__("bar")
                self.nodes = nodes
        
        elements = [
            BarElement([nodes[0], nodes[1]]),
            BarElement([nodes[1], nodes[2]])
        ]
        model.add_elements(elements)
        
        # Verify elements
        assert model.n_elements == 2
        for i, element in enumerate(elements):
            assert element.label == i
            assert element in model.elements
        
        # Test element-node associations
        for element in elements:
            for node in element.nodes:
                assert element in node._elements
    
    def test_node_element_association(self):
        """Test that nodes properly track their associated elements"""
        model = Model("test", "bar")
        
        # Create nodes
        n1 = Node((0, 0))
        n2 = Node((1, 0))
        n3 = Node((2, 0))
        
        model.add_nodes([n1, n2, n3])
        
        # Create elements
        class MockElement(Element):
            def __init__(self, nodes, label):
                super().__init__("bar")
                self.label = label
                self.nodes = nodes
        
        e1 = MockElement([n1, n2], 0)
        e2 = MockElement([n2, n3], 1)
        
        model.add_elements([e1, e2])
        
        # Check associations
        assert len(n1._elements) == 1
        assert e1 in n1._elements
        
        assert len(n2._elements) == 2
        assert e1 in n2._elements
        assert e2 in n2._elements
        
        assert len(n3._elements) == 1
        assert e2 in n3._elements