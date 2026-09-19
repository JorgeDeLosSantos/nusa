# -*- coding: utf-8 -*-
# ***********************************
#  Author: Pedro Jorge De Los Santos
#  E-mail: delossantosmfq@gmail.com
#  Blog: numython.github.io
#  License: MIT License
# ***********************************
from nusa import Node, Truss, TrussModel


def build_model():
    """Logan, Example 3.1: three-member truss."""
    E = 30e6  # psi
    A = 2.0  # in^2
    P = 10e3  # lbf

    model = TrussModel("Truss Model")
    n1 = Node((0, 0))
    n2 = Node((0, 120))
    n3 = Node((120, 120))
    n4 = Node((120, 0))

    elements = [
        Truss((n1, n2), E, A),
        Truss((n1, n3), E, A),
        Truss((n1, n4), E, A),
    ]

    model.add_nodes([n1, n2, n3, n4])
    model.add_elements(elements)
    model.add_force(n1, (0, -P))
    model.add_constraint(n2, ux=0, uy=0)
    model.add_constraint(n3, ux=0, uy=0)
    model.add_constraint(n4, ux=0, uy=0)
    model.solve()
    return model


def main():
    model = build_model()
    model.plot_model()
    model.plot_deformed_shape()
    model.show()


if __name__ == "__main__":
    main()
