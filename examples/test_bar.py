# -*- coding: utf-8 -*-
# ***********************************
#  Author: Pedro Jorge De Los Santos     
#  E-mail: delossantosmfq@gmail.com 
#  Web: labdls.blogspot.mx
#  License: MIT License
# ***********************************
import sys
sys.path.insert(0,'../') # Insert parent folder

from nusa.core import *
from nusa.model import *
from nusa.element import *

def test1():
    """
    Logan, D. (2007). A first course in the finite element analysis.
    Example 3.1, pp. 70.
    """
    # Input data 
    E1 = 30e6
    A1 = 1.0
    E2 = 15e6
    A2 = 2.0
    L = 30.0
    P = 3000.0
    # Model
    m1 = BarModel("Bar Model")
    # Nodes
    n1 = Node((0,0))
    n2 = Node((30,0))
    n3 = Node((60,0))
    n4 = Node((90,0))
    # Elements
    e1 = Bar((n1,n2),E1,A1,L)
    e2 = Bar((n2,n3),E1,A1,L)
    e3 = Bar((n3,n4),E2,A2,L)

    # Add elements 
    for nd in (n1,n2,n3,n4):
        m1.addNode(nd)
    for el in (e1,e2,e3):
        m1.addElement(el)

    m1.addForce(n2,(P,))
    m1.addConstraint(n1,ux=0) # fixed 
    m1.addConstraint(n4,ux=0) # fixed
    m1.solve() # Solve model
    
    print("Node | Displacements | Forces")
    for node in m1.getNodes():
        print "{0}\t{1}\t{2}".format(node.label, node.ux, node.fx)
        

def test2():
    """
    Kattan, P.I. (2003). MATLAB guide to finite elements, an interactive approach.
    Example 3.1, pp. 29
    """
    E = 210e6
    A = 0.003
    P = -10
    L1, L2 = 1.5, 1
    UX3 = 0.002
    
    m1 = BarModel("Bar model 02")
    
    n1 = Node((0,0))
    n2 = Node((1.5,0))
    n3 = Node((2.5,0))
    
    e1 = Bar((n1,n2),E,A,L1)
    e2 = Bar((n2,n3),E,A,L2)
    
    for nd in (n1,n2,n3): m1.addNode(nd)
    for el in (e1,e2): m1.addElement(el)
    
    m1.addForce(n2, (P,))
    m1.addConstraint(n1, ux=0.0)
    m1.addConstraint(n3, ux=UX3)
    m1.solve()
    
    print m1.KG
    print n1.ux, n2.ux, n3.ux
    print n1.fx, n2.fx, n3.fx
    print e1.sx, e2.sx
    
    
if __name__ == '__main__':
    test2()
