# -*- coding: utf-8 -*-
# ***********************************
#  Author: Pedro Jorge De Los Santos     
#  E-mail: delossantosmfq@gmail.com 
#  Web: labdls.blogspot.mx
#  License: MIT License
# ***********************************
from nusa.core import *
from nusa.model import *
from nusa.element import *

def test2():
    """
    Logan, D. (2007). A first course in the finite element analysis.
    Example 2.2, pp. 45.
    """
    P = 4e3
    k = 200e3
    # Model
    m2 = SpringModel("Spring Model 02")
    # Nodes
    n1 = Node((0,0),0)
    n2 = Node((0,0),1)
    n3 = Node((0,0),2)
    n4 = Node((0,0),3)
    n5 = Node((0,0),4)
    # Elements
    e1 = Spring((n1,n2),k)
    e2 = Spring((n2,n3),k)
    e3 = Spring((n3,n4),k)
    e4 = Spring((n4,n5),k)
    
    for nd in (n1,n2,n3,n4,n5):
        m2.addNode(nd)
    for el in (e1,e2,e3,e4):
        m2.addElement(el)
    
    m2.addForce(n4,(P,))
    m2.addConstraint(n1,ux=0)
    m2.addConstraint(n5,ux=0.02)
    
    m2.solve()


if __name__ == '__main__':
    test_2()
