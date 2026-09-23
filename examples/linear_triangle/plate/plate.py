# -*- coding: utf-8 -*-
# ***********************************
#  Author: Pedro Jorge De Los Santos     
#  E-mail: delossantosmfq@gmail.com 
#  Blog: numython.github.io
#  License: MIT License
# ***********************************
import numpy as np

from nusa import LinearTriangle, LinearTriangleModel, Node
from nusa.mesh import Modeler

modeler = Modeler()
modeler.add_rectangle((0,0),(0.3,0.3), esize=0.05)
nc, ec = modeler.generate_mesh()
x,y = nc[:,0], nc[:,1]

nodos = []
elementos = []

for k,nd in enumerate(nc):
    cn = Node((x[k],y[k]))
    nodos.append(cn)
    
for elm in ec:
    i, j, k = int(elm[0]), int(elm[1]), int(elm[2])
    ni, nj, nk = nodos[i], nodos[j], nodos[k]
    ce = LinearTriangle((ni,nj,nk),200e9, 0.3, 0.01)
    elementos.append(ce)

model = LinearTriangleModel()
for node in nodos: model.add_node(node)
for elm in elementos: model.add_element(elm)

minx = min(x)
maxx = max(x)

nnf = len([node for node in nodos if np.isclose(node.x, maxx)])
F = (6000./nnf)

for node in nodos:
    if np.isclose(node.x, minx):
        model.add_constraint(node, ux=0, uy=0)
    if np.isclose(node.x, maxx):
        model.add_force(node, (F,0))

model.plot_model()
model.solve()
# Plotting
model.plot_nodal_result("seqv") # von Mises Stress
#~ model.plot_nodal_result("exx")
model.show()


