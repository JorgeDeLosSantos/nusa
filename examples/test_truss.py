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
from nusa.graph import *

def test1():
    """
    Logan, D. (2007). A first course in the finite element analysis.
    Example 3.1, pp. 70.
    """
    # Input data 
    E = 30e6
    A = 2.0
    P = 10e3
    L = 10*(12.0)  # ft -> in
    L2 = np.sqrt(120**2 + 120**2)
    # Model
    m1 = TrussModel("Truss Model")
    # Nodes
    n1 = Node((0,0))
    n2 = Node((0,120))
    n3 = Node((120,120))
    n4 = Node((120,0))
    # Elements
    kdg = np.pi/180.0
    e1 = Truss((n1,n2),E,A,L,90*kdg)
    e2 = Truss((n1,n3),E,A,L2,45*kdg)
    e3 = Truss((n1,n4),E,A,L,0*kdg)

    # Add elements 
    for nd in (n1,n2,n3,n4):
        m1.addNode(nd)
    for el in (e1,e2,e3):
        m1.addElement(el)

    m1.buildGlobalMatrix()
    m1.addForce(n1,(0.0,-P))
    m1.addConstraint(n2,ux=0,uy=0) # fixed 
    m1.addConstraint(n3,ux=0,uy=0) # fixed
    m1.addConstraint(n4,ux=0,uy=0) # fixed
    m1.solve() # Solve model
    #~ plot_truss_model(m1)
    for node in m1.nodes.values():
        print node.uy
    #m1.report(out="cli")


if __name__ == '__main__':
    test1()
