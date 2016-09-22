# -*- coding: utf-8 -*-
# ***********************************
#  Author: Pedro Jorge De Los Santos     
#  E-mail: delossantosmfq@gmail.com 
#  Web: labdls.blogspot.mx
#  License: MIT License
# ***********************************
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
    #~ print m1.U



def test4():
    """
    Beer & Johnston. Mechanics of materials
    Problem 9.13 , pp. 568.
    """
    # Input data 
    E = 29e6
    I = 291 # W14x30 
    P = 35e3
    L1 = 5*12 # in
    L2 = 10*12 #in
    # Model
    m1 = BeamModel("Beam Model")
    # Nodes
    n1 = Node((0,0))
    n2 = Node((L1,0))
    n3 = Node((L1+L2,0))
    # Elements
    e1 = Beam((n1,n2),E,I,L1)
    e2 = Beam((n2,n3),E,I,L2)

    # Add elements 
    for nd in (n1,n2,n3): m1.addNode(nd)
    for el in (e1,e2): m1.addElement(el)
        
    m1.addForce(n2, (-P,))
    m1.addConstraint(n1, ux=0, uy=0) # fixed 
    m1.addConstraint(n3, uy=0) # fixed
    m1.solve() # Solve model
    #~ print m1.KG
    #~ print m1.U
    #~ print e1.m
    
def test5():
    """
    Beer & Johnston. Mechanics of materials
    Problem 9.75 , pp. 589.
    """
    # Input data 
    E = 29e6 # psi
    b, h = 2.0, 4.0  #in
    I = (1/12.0)*b*h**3
    w = 1e3/12.0 # lb/in
    L1 = 2*12 # in
    L2 = 3*12 #in
    P1 = -1e3 # lb
    P2 = -w*L2/2.0
    P3 = -w*L2/2.0
    M2 = -w*L2**2/12.0
    M3 = w*L2**2/12.0
    # Model
    m1 = BeamModel("Beam Model")
    
    # Nodes
    n1 = Node((0,0))
    n2 = Node((L1,0))
    n3 = Node((L1+L2,0))
    
    # Elements
    e1 = Beam((n1,n2),E,I,L1)
    e2 = Beam((n2,n3),E,I,L2)

    # Add elements 
    for nd in (n1,n2,n3): m1.addNode(nd)
    for el in (e1,e2): m1.addElement(el)
        
    m1.addForce(n1, (P1,))
    m1.addForce(n2, (P2,))
    m1.addForce(n3, (P3,))
    m1.addMoment(n2, (M2,))
    m1.addMoment(n3, (M3,))
    m1.addConstraint(n3, ux=0, uy=0, ur=0) # fixed
    m1.solve() # Solve model
    
    # Slope and deflection in n1
    print n1.uy, n1.ur


if __name__ == '__main__':
    test5()
