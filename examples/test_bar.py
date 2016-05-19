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
    n2 = Node((0,30))
    n3 = Node((0,60))
    n4 = Node((0,90))
    # Elements
    e1 = Bar((n1,n2),E1,A1,L)
    e2 = Bar((n2,n3),E1,A1,L)
    e3 = Bar((n3,n4),E2,A2,L)

    # Add elements 
    for nd in (n1,n2,n3,n4):
        m1.addNode(nd)
    for el in (e1,e2,e3):
        m1.addElement(el)

    m1.buildGlobalMatrix()
    m1.addForce(n2,(P,)) 
    m1.addConstraint(n1,ux=0) # fixed 
    m1.addConstraint(n4,ux=0) # fixed
    m1.solve() # Solve model
    m1.report(out="cli") # Print report


if __name__ == '__main__':
    test1()
