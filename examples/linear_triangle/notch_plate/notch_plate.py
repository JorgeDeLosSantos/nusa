# -*- coding: utf-8 -*-
# ***********************************
#  Author: Pedro Jorge De Los Santos     
#  E-mail: delossantosmfq@gmail.com 
#  Blog: numython.github.io
#  License: MIT License
# ***********************************
import matplotlib.pyplot as plt
import numpy as np

from nusa import LinearTriangle, LinearTriangleModel, Material, Node, plot_model
from nusa.mesh import Modeler

modeler = Modeler()
a = modeler.add_poly((0,0),(1,0),(1,1),(0.6,1),(0.5,0.9),(0.4,1),(0,1), esize=0.08)
nc, ec = modeler.generate_mesh()
x,y = nc[:,0], nc[:,1]

nodos = []
elementos = []

for k,nd in enumerate(nc):
    cn = Node((x[k],y[k]))
    nodos.append(cn)
    
material = Material(E=200e9, nu=0.3)

for elm in ec:
    i, j, k = int(elm[0]), int(elm[1]), int(elm[2])
    ni, nj, nk = nodos[i], nodos[j], nodos[k]
    ce = LinearTriangle((ni,nj,nk), material=material, thickness=0.1)
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

plot_model(model)
result = model.solve()
result.plot_nodal_field("sxx")
plt.show()

