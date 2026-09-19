# -*- coding: utf-8 -*-
# ***********************************
#  Author: Pedro Jorge De Los Santos
#  E-mail: delossantosmfq@gmail.com
#  Blog: numython.github.io
#  License: MIT License
# ***********************************
from nusa import Node, Truss, TrussModel


def build_model():
    """Kattan, Problem 5.1."""
    E = 210e9
    A = 0.005

    nodes = [
        Node((0, 0)),
        Node((5, 7)),
        Node((5, 0)),
        Node((10, 7)),
        Node((10, 0)),
        Node((15, 0)),
    ]
    n1, n2, n3, n4, n5, n6 = nodes

    elements = [
        Truss((n1, n2), E, A),
        Truss((n1, n3), E, A),
        Truss((n2, n3), E, A),
        Truss((n2, n4), E, A),
        Truss((n2, n5), E, A),
        Truss((n3, n5), E, A),
        Truss((n4, n5), E, A),
        Truss((n4, n6), E, A),
        Truss((n5, n6), E, A),
    ]

    model = TrussModel("Example 02")
    model.add_nodes(nodes)
    model.add_elements(elements)
    model.add_constraint(n1, ux=0, uy=0)
    model.add_constraint(n6, ux=0, uy=0)
    model.add_force(n2, (20e3, 0))
    model.solve()
    return model


def main():
    model = build_model()
    model.plot_model()
    model.plot_deformed_shape()
    model.show()


if __name__ == "__main__":
    main()
