# -*- coding: utf-8 -*-
# ***********************************
#  Author: Pedro Jorge De Los Santos
#  E-mail: delossantosmfq@gmail.com
#  License: MIT License
# ***********************************

import matplotlib.pyplot as plt

from nusa import Node, Truss, TrussModel, plot_model


def build_model():
    """Build Logan's three-member truss example."""
    E = 30e6
    A = 2.0
    P = 10e3

    model = TrussModel("Truss Model")
    n1 = Node((0.0, 0.0))
    n2 = Node((0.0, 120.0))
    n3 = Node((120.0, 120.0))
    n4 = Node((120.0, 0.0))

    model.add_nodes([n1, n2, n3, n4])
    model.add_elements([
        Truss((n1, n2), E, A),
        Truss((n1, n3), E, A),
        Truss((n1, n4), E, A),
    ])
    model.add_force(n1, (0.0, -P))
    model.add_constraint(n2, ux=0.0, uy=0.0)
    model.add_constraint(n3, ux=0.0, uy=0.0)
    model.add_constraint(n4, ux=0.0, uy=0.0)
    return model


def main():
    model = build_model()
    plot_model(model)

    result = model.solve()
    result.plot_deformed_shape()
    plt.show()


if __name__ == "__main__":
    main()
