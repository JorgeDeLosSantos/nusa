# -*- coding: utf-8 -*-
# ***********************************
#  Author: Pedro Jorge De Los Santos     
#  E-mail: delossantosmfq@gmail.com 
#  Web: labdls.blogspot.mx
#  License: MIT License
# ***********************************
from __future__ import division
import numpy as np
from nusa import *
import itertools
import matplotlib.pyplot as plt

def pairwise(iterable):
    #~ "s -> (s0,s1), (s1,s2), (s2, s3), ..."
    a, b = itertools.tee(iterable)
    next(b, None)
    return itertools.izip(a, b)


# Input data 
E = 29e6 # psi
I = 10
L = 10
P = 10e3

nelm = 3
parts = np.linspace(0,L,nelm)

nodos = []
for xc in parts:
    cn = Node((xc,0))
    nodos.append(cn)

elementos = []
for x in pairwise(nodos):
    ni,nj = x[0], x[1]
    ce = Beam((ni,nj),E,I,L/(nelm-1))
    elementos.append(ce)

m = BeamModel()

for n in nodos: m.addNode(n)
for e in elementos: m.addElement(e)

m.addConstraint(nodos[0], ux=0, uy=0, ur=0)
m.addForce(nodos[-1], (-P,))
m.solve()

m.plot_disp()

xx = np.linspace(0,L)
d = ((-P*xx**2.0)/(6.0*E*I))*(3*L - xx)
plt.plot(xx,d)
plt.axis("auto")
plt.xlim(0,L+1)

m.show()



