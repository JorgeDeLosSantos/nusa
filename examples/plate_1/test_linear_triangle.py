# -*- coding: utf-8 -*-
# ***********************************
#  Author: Pedro Jorge De Los Santos     
#  E-mail: delossantosmfq@gmail.com 
#  Web: labdls.blogspot.mx
#  License: MIT License
# ***********************************
import numpy as np
from nusa._experimental import *


nc = np.loadtxt("nlist")
ec = np.loadtxt("elist")[:,6:-1]
x,y = nc[:,1], nc[:,2]

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

for n in (1,12,14,15,16):
    m.addConstraint(nodos[n-1], ux=0, uy=0)
for n in (6,10):
    m.addForce(nodos[n-1], (5000,1000))

m.solve()
#~ for n in m.getNodes():
    #~ print n.label+1, n.ux

#~ for e in m.getElements():
    #~ print e.sx, e.sy
#~ m.plot_sxx()
#~ m.plot_syy()
m.plot_seqv()

#~ n1 = Node((0,0))
#~ n2 = Node((0.5,0))
#~ n3 = Node((0.5,0.25))
#~ n4 = Node((0,0.25))
#~ n5 = Node((0,0.5))
#~ e1 = LinearTriangle((n1,n3,n4),210e6,0.3,0.025)
#~ e2 = LinearTriangle((n1,n2,n3),210e6,0.3,0.025)
#~ e3 = LinearTriangle((n4,n3,n5),210e6,0.3,0.025)
#~ m = LinearTriangleModel()
#~ for node in (n1,n2,n3,n4,n5): m.addNode(node)
#~ for elm in (e1,e2,e3): m.addElement(elm)
#~ m.addConstraint(n1, ux=0, uy=0)
#~ m.addConstraint(n4, ux=0, uy=0)
#~ m.addForce(n2, (9375,0))
#~ m.addForce(n3, (9375,0))
#~ m.solve()
#~ print e1.sx, e1.sy
#~ print e2.sx, e2.sy
#~ print e3.sx, e3.sy
#~ m.plot_stresses()
