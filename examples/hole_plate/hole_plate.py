# -*- coding: utf-8 -*-
# ***********************************
#  Author: Pedro Jorge De Los Santos     
#  E-mail: delossantosmfq@gmail.com 
#  Web: labdls.blogspot.mx
#  License: MIT License
# ***********************************
import numpy as np
from nusa._experimental import *

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

#~ sx = []
#~ for n in m.getNodes():
    #~ sx.append(n.sx)
#~ print max(sx)
#~ m.plot_model()
#~ m.plot_ux()
#~ m.plot_uy()
#~ m.plot_usum()
#~ m.plot_sxx()
#~ m.plot_syy()
#~ m.plot_sxy()
#~ m.plot_seqv()
#~ m.plot_syy() 
#~ m.plot_seqv()
print elementos[3].ex
m.plot_eex()
m.show()
