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
    n2 = Node((10*12,0))
    n3 = Node((20*12,0))
    n4 = Node((30*12,0))
    n5 = Node((40*12,0))
    # Elements
    e1 = Beam((n1,n2),E,I,L)
    e2 = Beam((n2,n3),E,I,L)
    e3 = Beam((n3,n4),E,I,L)
    e4 = Beam((n4,n5),E,I,L)

    # Add elements 
    for nd in (n1,n2,n3,n4,n5): m1.add_node(nd)
    for el in (e1,e2,e3,e4): m1.add_element(el)

    m1.add_force(n2,(-P,))
    m1.add_force(n4,(-P,))
    m1.add_constraint(n1, ux=0,uy=0,ur=0) # fixed 
    m1.add_constraint(n5, ux=0,uy=0,ur=0) # fixed
    m1.add_constraint(n3, uy=0, ur=0) # fixed
    m1.add_constraint(n2, ur=0)
    m1.add_constraint(n4, ur=0)
    m1.solve() # Solve model
    m1.plot_disp()
    m1.show()


if __name__ == '__main__':
    test1()
