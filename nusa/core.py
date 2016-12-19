# -*- coding: utf-8 -*-
# ***********************************
#  Author: Pedro Jorge De Los Santos    
#  E-mail: delossantosmfq@gmail.com 
#  Blog: numython.github.io
#  License: MIT License
# ***********************************

import numpy as np

#~ ===========================  MODEL  ===========================
class Model(object):
    """
    Superclass for all FEA models
    """
    def __init__(self,name,mtype):
        self.mtype = mtype # Model type
        self.name = name # Name 
        self.nodes = {} # Dictionary for nodes {number: NodeObject}
        self.elements = {} # Dictionary for elements {number: ElementObject}
        
    def addNode(self,node):
        """
        Add node to current model
        
        *node* :   :class:`~nusa.core.Node`
            Node instance
        
        """
        current_label = self.getNumberOfNodes()
        if node.label is "":
            node.setLabel(current_label)
        self.nodes[node.label] = node
        
    def addElement(self,element):
        """
        Add element to current model
        
        *element* :  :class:`~nusa.core.Element`
            Element instance 
        
        ::
        
            m1 = Model()
            e1 = Bar(Node((0,0),0),Node((1,0),0))
            m1.addElement(e1)
        
        """
        if self.mtype != element.etype:
            raise ValueError("Element type must be "+self.mtype)
        current_label = self.getNumberOfElements()
        if element.label is "":
            element.setLabel(current_label)
        self.elements[element.label] = element

    def getNumberOfNodes(self):
        """
        Returns the number of nodes
        """
        return len(self.nodes)
        
    def getNumberOfElements(self):
        """
        Returns the number of nodes
        """
        return len(self.elements)
        
    def getNodes(self):
        """
        Returns a list of Node objects
        """
        return self.nodes.values()
        
    def getElements(self):
        """
        Returns a list of Element objects
        """
        return self.elements.values()
    
    def __str__(self):
        custom_str = ("Model: "+self.name+"\nNodes: "+str(self.getNumberOfNodes())+
        "\nElements: "+str(self.getNumberOfElements()))
        return custom_str
            

#~ =========================== ELEMENT ===========================

class Element(object):
    """
    Superclass for all Elements
    """
    def __init__(self,etype):
        self.etype = etype # element type
        self.label = "" # label (reassignment -> Model.addElement)
        self._fx = 0.0
        self._fy = 0.0
        self._sx = 0.0
        self._sy = 0.0
        
    @property
    def fx(self):
        return self._fx
        
    @fx.setter
    def fx(self,val):
        self._fx = val
        
    @property
    def fy(self):
        return self._fy
        
    @fy.setter
    def fy(self,val):
        self._fy = val
        
    def setLabel(self,label):
        """
        Set the label property
        
        *label* : int
            Label, must be an integer
        """
        self.label = label
        
    def setElementForces(self,fx=0.0,fy=0.0):
        """
        Set element forces
        
        *fx* : float
            Force in x-dir
        *fy* : float
            Force in y-dir
        
        Normally this method is used by the `solve` method to 
        update computed element-forces.
        """
        self._fx = fx
        self._fy = fy
        
    def getElementForces(self):
        """
        Returns a tuple with element forces:  (fx, fy)
        """
        return self._fx, self._fy
        
    def __str__(self):
        _str = str(self.__class__)
        return _str


#~ =========================== NODE ===========================

class Node(object):
    """
    Class for node object.
    
    *coordinates* : `tuple`, `list`
        Coordinates of node
    
    *label* : int
        Label of node
        
    ::
    
        n1 = Node((0,0),0)
        n2 = Node((0,0),1)
    
    """
    def __init__(self,coordinates):
        self.coordinates = coordinates
        self.x = coordinates[0] # usable prop
        self.y = coordinates[1] # usable prop
        self._label = ""
        self._ux = np.nan
        self._uy = np.nan
        self._ur = np.nan
        self._fx = 0.0
        self._fy = 0.0
        self._m = 0.0
        
    @property
    def label(self):
        return self._label
        
    @label.setter
    def label(self,val):
        """
        Experimental setter for adjust range of labels: TO DO
        """
        self._label = val
        
    @property
    def ux(self):
        return self._ux
    
    @ux.setter
    def ux(self,val):
        if True:#type(val) in [int,float]:
            self._ux = val
        else:
            raise ValueError("Value must be float or int")
    
    @property
    def uy(self):
        return self._uy
    
    @uy.setter
    def uy(self,val):
        if True:#type(val) in [int,float]:
            self._uy = val
        else:
            raise ValueError("Value must be float or int")
    
    @property
    def ur(self):
        return self._ur
    
    @ur.setter
    def ur(self,val):
        if True:#type(val) in [int,float]:
            self._ur = val
        else:
            raise ValueError("Value must be float or int")
        
    @property
    def fx(self):
        return self._fx
    
    @fx.setter
    def fx(self,val):
        self._fx = val
    
    @property
    def fy(self):
        return self._fy
    
    @fy.setter
    def fy(self,val):
        self._fy = val
        
    @property
    def m(self):
        return self._m
    
    @m.setter
    def m(self,val):
        self._m = val
        
    def getLabel(self):
        return self._label
    
    def setLabel(self,label):
        self._label = label
    
    def getDisplacements(self):
        return self._ux,self._uy,self._ur
        
    def setDisplacements(self,ux=np.nan, uy=np.nan, ur=np.nan):
        self._ux = ux
        self._uy = uy
        self._ur = ur
    
    def getForces(self):
        return (self._fx,self._fy)
    
    def setForces(self,fx=np.nan,fy=np.nan):
        self._fx = fx
        self._fy = fy
        
    def __str__(self):
        _str = self.__class__
        _str = "%s\nU:(%g,%g)\n"%(_str,self.ux, self.uy)
        _str = "%sF:(%g,%g)"%(_str,self.fx,self.fy)
        return _str
        

if __name__=='__main__':
    pass
