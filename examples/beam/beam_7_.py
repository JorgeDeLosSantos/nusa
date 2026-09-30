# -*- coding: utf-8 -*-
# ***********************************
#  Author: Pedro Jorge De Los Santos
#  E-mail: delossantosmfq@gmail.com
#  License: MIT License
# ***********************************

import itertools

import matplotlib.pyplot as plt
import numpy as np

from nusa import Beam, BeamModel, Node


def pairwise(iterable):
    a, b = itertools.tee(iterable)
    next(b, None)
    return zip(a, b)


E = 29e6
I = 10.0
L = 10.0
P = 10e3

nn = 20
parts = np.linspace(0.0, L, nn)
nodes = [Node((x, 0.0)) for x in parts]
elements = [Beam((ni, nj), E, I) for ni, nj in pairwise(nodes)]

model = BeamModel()
model.add_nodes(nodes)
model.add_elements(elements)
model.add_constraint(nodes[0], uy=0.0)
model.add_constraint(nodes[-1], uy=0.0)
model.add_force(nodes[5], (-P,))

result = model.solve()
result.plot_deformed_shape(scale=1.0)

xa = np.linspace(0.0, nodes[5].x)
xb = np.linspace(nodes[5].x, nodes[-1].x)
a = nodes[5].x
b = nodes[-1].x - nodes[5].x
da = ((-P * b * xa) / (6 * L * E * I)) * (L**2 - b**2 - xa**2)
db = ((-P * a * (xb - L)) / (6 * L * E * I)) * (
    a**2 - 2 * L * xb + xb**2
)
plt.plot(xa, da)
plt.plot(xb, db)
plt.axis("auto")
plt.xlim(-1, L + 1)
plt.show()

print(result.displacement(nodes[4])["uy"])
