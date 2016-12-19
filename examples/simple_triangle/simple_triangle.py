# -*- coding: utf-8 -*-
# ***********************************
#  Author: Pedro Jorge De Los Santos     
#  E-mail: delossantosmfq@gmail.com 
#  Web: labdls.blogspot.mx
#  License: MIT License
# ***********************************
import numpy as np
from nusa._experimental import LinearTriangle, LinearTriangleModel, Node

# Creando el modelo
m = LinearTriangleModel()
# Creando nodos
n1 = Node((0,0))
n2 = Node((1,0.5))
n3 = Node((0,1))
# Creando elemento
e1 = LinearTriangle((n1,n2,n3),200e9,0.3,0.1)
for n in (n1,n2,n3): m.addNode(n)
m.addElement(e1)
# Agregando condiciones de frontera
m.addConstraint(n1, ux=0, uy=0)
m.addConstraint(n3, ux=0, uy=0)
m.addForce(n2, (1000,0))
#~ m.plot_model()
m.solve()
m.plot_ux()
m.plot_uy()
m.plot_usum()
import matplotlib.pyplot as plt
plt.show()
#~ for n in m.getNodes():
    #~ print n.ux, n.uy
