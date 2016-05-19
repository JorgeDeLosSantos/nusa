#~ # -*- coding: utf-8 -*-
#~ # ***********************************
#~ #  Author: Pedro Jorge De Los Santos     
#~ #  E-mail: delossantosmfq@gmail.com 
#~ #  Web: labdls.blogspot.mx
#~ #  License: MIT License
#~ # ***********************************
#~ import sys
#~ sys.path.insert(0, '../')
#~ 
#~ import matplotlib.pyplot as plt
#~ from nusa.core import *
#~ from nusa.lib import STEEL_1018,ALUMINIUM_6061
#~ from nusa.model import SpringModel
#~ from nusa.element import Spring
#~ from nusa.graph import plot_spring_model
#~ 
#~ P = 5000.0
#~ # Model
#~ m1 = SpringModel("2D Model")
#~ # Nodes
#~ n1 = Node((0,0))
#~ n2 = Node((0,0))
#~ n3 = Node((0,0))
#~ n4 = Node((0,0))
#~ # Elements
#~ e1 = Spring((n1,n3),1000.0)
#~ e2 = Spring((n3,n4),2000.0)
#~ e3 = Spring((n4,n2),3000.0)
#~ 
#~ # Add elements 
#~ for nd in (n1,n2,n3,n4):
    #~ m1.addNode(nd)
#~ for el in (e1,e2,e3):
    #~ m1.addElement(el)
#~ 
#~ m1.buildGlobalMatrix()
#~ m1.addForce(n4,(P,))
#~ m1.addConstraint(n1,ux=0)
#~ m1.solve()
#~ print m1.U
