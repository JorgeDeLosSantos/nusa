# -*- coding: utf-8 -*-
# ***********************************
#  Author: Pedro Jorge De Los Santos     
#  E-mail: delossantosmfq@gmail.com 
#  Blog: numython.github.io
#  License: MIT License
# ***********************************
import numpy as np
from nusa import *
from nusa.io import *

nc, ec = read_msh("geom2.msh")
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

m.plot_model()
m.show()
