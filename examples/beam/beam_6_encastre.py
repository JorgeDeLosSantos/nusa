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

nelm = 10
parts = np.linspace(0.0, L, nelm + 1)
nodes = [Node((x, 0.0)) for x in parts]
elements = [Beam((ni, nj), E, I) for ni, nj in pairwise(nodes)]

model = BeamModel()
model.add_nodes(nodes)
model.add_elements(elements)
model.add_constraint(nodes[0], uy=0.0, ur=0.0)
model.add_force(nodes[-1], (-P,))

result = model.solve()
result.plot_deformed_shape(scale=1.0, label="Approx.")

xx = np.linspace(0.0, L)
d = ((-P * xx**2.0) / (6.0 * E * I)) * (3.0 * L - xx)
plt.plot(xx, d, label="Classic")
plt.legend()
plt.axis("auto")
plt.xlim(0.0, L + 1.0)
plt.show()
