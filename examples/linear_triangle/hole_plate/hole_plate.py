# -*- coding: utf-8 -*-
# ***********************************
#  Author: Pedro Jorge De Los Santos     
#  E-mail: delossantosmfq@gmail.com 
#  Blog: numython.github.io
#  License: MIT License
# ***********************************
import numpy as np
from nusa import *
from nusa.mesh import *

m = Modeler()
a = m.add_rectangle((0,0),(2,1), 0.2)
b = m.add_rectangle((0.3,0.3),(0.7,0.7), 0.02)
#~ b = m.add_circle((0.5,0.5), 0.1, 0.02)
m.substract_surfaces(a,b)
nc, ec = m.generate_mesh()
x,y = nc[:,0], nc[:,1]

nodos = []
elementos = []

for k,nd in enumerate(nc):
    cn = Node((x[k],y[k]))
    nodos.append(cn)
    
for k,elm in enumerate(ec):
    i,j,m = int(elm[0]),int(elm[1]),int(elm[2])
    ni,nj,nm = nodos[i],nodos[j],nodos[m]
    ce = LinearTriangle((ni,nj,nm),200e9,0.3,1)
    elementos.append(ce)

m = LinearTriangleModel()
for node in nodos: m.addNode(node)
for elm in elementos: m.addElement(elm)

minx = min(x)
maxx = max(x)

for node in nodos:
    if node.x == minx:
        m.addConstraint(node, ux=0, uy=0)
    if node.x == maxx:
        m.addForce(node, (1e4,0))

m.plot_model()
m.solve()
# Plotting
m.plot_nsol("sxy")
m.show()
