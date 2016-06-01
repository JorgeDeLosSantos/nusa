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
    Model for finite element analysis
    """
    def __init__(self,name,mtype):
        self.mtype = mtype
        self.name = name
        self.nodes = {}
        self.elements = {}
        
    def addNode(self,node):
        current_label = self.getNumberOfNodes()
        if node.label is "":
            node.setLabel(current_label)
        self.nodes[node.label] = node
        
    def addElement(self,element):
        """
        Add element to current model
        
        *element* :  :class:`~core.base.Element`
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
        self.__fx = 0.0
        self.__fy = 0.0
        
    @property
    def fx(self):
        return self.__fx
        
    @fx.setter
    def fx(self,val):
        self.__fx = val
        
    @property
    def fy(self):
        return self.__fy
        
    @fy.setter
    def fy(self,val):
        self.__fy = val
        
    def setLabel(self,label):
        self.label = label
        
    def setElementForces(self,fx=0.0, fy=0.0):
        self.__fx = fx
        self.__fy = fy
        
    def getElementForces(self):
        return (self.__fx, self.__fy)
        
    def __str__(self):
        return self.etype, self.nodes, self.elements


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
        self.__label = label
        self.__ux = np.nan
        self.__uy = np.nan
        self.__fx = 0.0
        self.__fy = 0.0
        
    @property
    def label(self):
        return self.__label
        
    @label.setter
    def label(self,val):
        """
        Experimental setter for adjust range of labels
        """
        self.__label = val
        
    @property
    def ux(self):
        return self.__ux
    
    @ux.setter
    def ux(self,val):
        if True:#type(val) in [int,float]:
            self.__ux = val
        else:
            raise ValueError("Value must be float or int")
    
    @property
    def uy(self):
        return self.__uy
    
    @uy.setter
    def uy(self,val):
        if True:#type(val) in [int,float]:
            self.__uy = val
        else:
            raise ValueError("Value must be float or int")
        
    @property
    def fx(self):
        return self.__fx
    
    @fx.setter
    def fx(self,val):
        self.__fx = val
    
    @property
    def fy(self):
        return self.__fy
    
    @fy.setter
    def fy(self,val):
        self.__fy = val
        
    def getLabel(self):
        return self.__label
    
    def setLabel(self,label):
        self.__label = label
    
    def getDisplacements(self):
        return (self.__ux,self.__uy)
        
    def setDisplacements(self,ux=np.nan,uy=np.nan):
        self.__ux = ux
        self.__uy = uy
    
    def getForces(self):
        return (self.__fx,self.__fy)
    
    def setForces(self,fx=np.nan,fy=np.nan):
        self.__fx = fx
        self.__fy = fy
        
    def __str__(self):
        _str = "Node\nCoords:(%g,%g)\n"%(self.ux, self.uy)
        _str = "%sForces:(%g,%g)"%(_str,self.fx,self.fy)
        return _str
        

if __name__=='__main__':
    pass
