# -*- coding: utf-8 -*-
# ***********************************
#  Author: Pedro Jorge De Los Santos     
#  E-mail: delossantosmfq@gmail.com 
#  Blog: numython.github.io
#  License: MIT License
# ***********************************
import numpy as np
from nusa import *

nc = np.loadtxt("nodos")
ec = np.loadtxt("elementos")
x,y = nc[:,0], nc[:,1]

nodos = []
elementos = []

for k,nd in enumerate(nc):
    cn = Node((x[k],y[k]))
    nodos.append(cn)
    
for k,elm in enumerate(ec):
    i,j,m = int(elm[0]-1),int(elm[1]-1),int(elm[2]-1)
    ni,nj,nm = nodos[i],nodos[j],nodos[m]
    ce = LinearTriangle((ni,nj,nm),200e9,0.3,1)
    elementos.append(ce)

m = LinearTriangleModel()
for node in nodos: m.addNode(node)
for elm in elementos: m.addElement(elm)

for n in (1,22,32,33,34,35,36,37,38,39,40):
    m.addConstraint(nodos[n-1], ux=0, uy=0)
for n in [2]+range(12,22):
    m.addForce(nodos[n-1], (1000/11.,0))

#~ m.plot_model()
m.solve()
# Plotting
m.plot_nsol("sxx")
#~ m.plot_esol("sxx")
m.show()
