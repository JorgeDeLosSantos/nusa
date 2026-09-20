# -*- coding: utf-8 -*-
# ***********************************
#  Author: Pedro Jorge De Los Santos     
#  E-mail: delossantosmfq@gmail.com 
#  Blog: numython.github.io
#  License: MIT License
# ***********************************
from nusa import LinearTriangle, LinearTriangleModel, Node
from nusa.mesh import Modeler

modeler = Modeler()
a = modeler.add_rectangle((0,0),(1,1), esize=0.1)
b = modeler.add_circle((0.5,0.5),0.15, esize=0.02)
modeler.subtract_surfaces(a,b)
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
    ce = LinearTriangle((ni,nj,nk),200e9, 0.3, 0.1)
    elementos.append(ce)

model = LinearTriangleModel()
for node in nodos: model.add_node(node)
for elm in elementos: model.add_element(elm)

minx = min(x)
maxx = max(x)

nnf = len([node for node in nodos if node.x==maxx])
F = (6000./nnf)

for node in nodos:
    if node.x == minx:
        model.add_constraint(node, ux=0, uy=0)
    if node.x == maxx:
        model.add_force(node, (F,0))

model.plot_model()
model.solve()
# Plotting
model.plot_element_result("sxx") # element solution
model.plot_nodal_result("sxx") # nodal solution
model.show()

w,d,t = 1, 0.3, 0.1
s0 = 6000./((w-d)*t)
kt = 2.35
s = kt*s0
print("Theory: {0}".format(s))

