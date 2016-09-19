# -*- coding: utf-8 -*-
# ***********************************
#  Author: Pedro Jorge De Los Santos     
#  E-mail: delossantosmfq@gmail.com 
#  Blog: numython.blogspot.mx
#  License: MIT License
# ***********************************

import numpy as np
import numpy.linalg as la
import templates as tmp

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
        Add element to current model
        
        *node* :   :class:`~nusa.core.Node`
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
        return len(self.nodes)
        
    def getNumberOfElements(self):
        return len(self.elements)
        
    def getNodes(self):
        return self.nodes.values()
        
    def getElements(self):
        return self.elements.values()
    
    def __str__(self):
        custom_str = ("Model: "+self.name+"\nNodes: "+str(self.getNumberOfNodes())+
        "\nElements: "+str(self.getNumberOfElements()))
        return custom_str
            

class Element(object):
    """
    Superclass for Element objects 
    """
    def __init__(self,etype,label):
        self.etype = etype
        self.label = label
        self._fx = 0.0
        self._fy = 0.0
        
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
        self.label = label
        
    def setElementForces(self,fx=0.0,fy=0.0):
        self._fx = fx
        self._fy = fy
        
    def getElementForces(self):
        return self._fx, self._fy
        
    def __str__(self):
        _str = str(self.__class__)
        return _str


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
    def __init__(self,coordinates,label=""):
        self.coordinates = coordinates
        self.x = coordinates[0] # usable prop
        self.y = coordinates[1] # usable prop
        self._label = label
        self._ux = np.nan
        self._uy = np.nan
        self._fx = 0.0
        self._fy = 0.0
        
    @property
    def label(self):
        return self._label
        
    @label.setter
    def label(self,val):
        """
        Experimental setter for adjust range of labels
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
        
    def getLabel(self):
        return self._label
    
    def setLabel(self,label):
        self._label = label
    
    def getDisplacements(self):
        return self._ux,self._uy
        
    def setDisplacements(self,ux=np.nan,uy=np.nan):
        self._ux = ux
        self._uy = uy
    
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
