# -*- coding: utf-8 -*-
# ***********************************
#  Author: Pedro Jorge De Los Santos
#  E-mail: delossantosmfq@gmail.com
#  License: MIT License
# ***********************************

import matplotlib.pyplot as plt

from nusa import Beam, BeamModel, Node


def test1():
    """Logan (2007), Example 4.2."""
    E = 30e6
    I = 500.0
    P = 10e3
    L = 10 * 12.0

    model = BeamModel("Beam Model")
    nodes = [Node((k * L, 0.0)) for k in range(5)]
    elements = [
        Beam((nodes[k], nodes[k + 1]), E, I)
        for k in range(4)
    ]

    model.add_nodes(nodes)
    model.add_elements(elements)
    model.add_force(nodes[1], (-P,))
    model.add_force(nodes[3], (-P,))
    model.add_constraint(nodes[0], uy=0.0, ur=0.0)
    model.add_constraint(nodes[4], uy=0.0, ur=0.0)
    model.add_constraint(nodes[2], uy=0.0, ur=0.0)

    result = model.solve()

    print(result.displacement(nodes[1])["uy"])
    result.plot_deformed_shape(scale=100.0)
    plt.show()
    return result


if __name__ == "__main__":
    test1()
