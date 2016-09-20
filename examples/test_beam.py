# -*- coding: utf-8 -*-
# ***********************************
#  Author: Pedro Jorge De Los Santos     
#  E-mail: delossantosmfq@gmail.com 
#  Web: labdls.blogspot.mx
#  License: MIT License
# ***********************************
import sys
sys.path.insert(0,'../') # Insert parent folder

import numpy as np
from nusa.core import *
from nusa.model import *
from nusa.element import *

def test1():
    """
    Logan, D. (2007). A first course in the finite element analysis.
    Example 4.2 , pp. 166.
    """
    # Input data 
    E = 30e6
    I = 500.0
    P = 10e3
    L = 10*(12.0)  # ft -> in
    # Model
    m1 = BeamModel("Beam Model")
    # Nodes
    n1 = Node((0,0))
    n2 = Node((0,10*12))
    n3 = Node((0,20*12))
    n4 = Node((0,30*12))
    n5 = Node((0,40*12))
    # Elements
    e1 = Beam((n1,n2),E,I,L)
    e2 = Beam((n2,n3),E,I,L)
    e3 = Beam((n3,n4),E,I,L)
    e4 = Beam((n4,n5),E,I,L)

    # Add elements 
    for nd in (n1,n2,n3,n4,n5): m1.addNode(nd)
    for el in (e1,e2,e3,e4): m1.addElement(el)

    m1.addForce(n2,(-P,))
    m1.addForce(n4,(-P,))
    m1.addConstraint(n1, ux=0,uy=0,ur=0) # fixed 
    m1.addConstraint(n5, ux=0,uy=0,ur=0) # fixed
    m1.addConstraint(n3, uy=0, ur=0) # fixed
    m1.addConstraint(n2, ur=0)
    m1.addConstraint(n4, ur=0)
    m1.solve() # Solve model
    print m1.U
    print m1.NF

def test2():
    """
    Logan, D. (2007). A first course in the finite element analysis.
    Example 4.4, pp. 171.
    """
    # Input data 
    E = 210e9
    I = 4e-4
    P = 10e3
    M = 20e3
    L = 3
    # Model
    m1 = BeamModel("Beam Model")
    # Nodes
    n1 = Node((0,0))
    n2 = Node((3,0))
    n3 = Node((6,0))
    # Elements
    e1 = Beam((n1,n2),E,I,L)
    e2 = Beam((n2,n3),E,I,L)

    # Add elements 
    for nd in (n1,n2,n3):
        m1.addNode(nd)
    for el in (e1,e2):
        m1.addElement(el)
        
    m1.addForce(n2, (-P,))
    m1.addMoment(n2, (M,))
    m1.addConstraint(n1, ux=0, uy=0, ur=0) # fixed 
    m1.addConstraint(n3, ux=0, uy=0, ur=0) # fixed
    m1.solve() # Solve model
    print m1.KG
    print m1.U
    print m1.NF


def test3():
    """
    Kattan, P. (XXXX).
    Example 7.1, pp. 109.
    """
    # Input data 
    E = 210e9
    I = 60e-6
    P = 20e3
    L = 2.0
    # Model
    m1 = BeamModel("Beam Model")
    # Nodes
    n1 = Node((0,0))
    n2 = Node((2,0))
    n3 = Node((4,0))
    # Elements
    e1 = Beam((n1,n2),E,I,L)
    e2 = Beam((n2,n3),E,I,L)

    # Add elements 
    for nd in (n1,n2,n3): m1.addNode(nd)
    for el in (e1,e2): m1.addElement(el)
        
    m1.addForce(n2, (-P,))
    m1.addConstraint(n1, ux=0, uy=0, ur=0) # fixed 
    m1.addConstraint(n3, ux=0, uy=0) # fixed
    m1.solve() # Solve model
    #~ print m1.KG
    print m1.U
    print m1.NF
    print n1
    #~ print m1.NF



if __name__ == '__main__':
    test3()
