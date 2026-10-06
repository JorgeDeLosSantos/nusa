# -*- coding: utf-8 -*-
# ***********************************
#  Author: Pedro Jorge De Los Santos
#  E-mail: delossantosmfq@gmail.com
#  Blog: numython.github.io
#  License: MIT License
# ***********************************
import matplotlib.pyplot as plt

from nusa import LinearTriangle, LinearTriangleModel, Node, plot_model


def build_model():
    n1 = Node((0, 0))
    n2 = Node((0.5, 0))
    n3 = Node((0.5, 0.25))
    n4 = Node((0, 0.25))
    n5 = Node((0, 0.5))

    elements = [
        LinearTriangle((n1, n3, n4), 210e6, 0.3, 0.025),
        LinearTriangle((n1, n2, n3), 210e6, 0.3, 0.025),
        LinearTriangle((n4, n3, n5), 210e6, 0.3, 0.025),
    ]

    model = LinearTriangleModel("Simple plate")
    model.add_nodes([n1, n2, n3, n4, n5])
    model.add_elements(elements)
    model.add_constraint(n1, ux=0, uy=0)
    model.add_constraint(n4, ux=0, uy=0)
    model.add_constraint(n5, ux=0, uy=0)
    model.add_force(n2, (9375, 0))
    model.add_force(n3, (9375, 0))
    return model


def main():
    model = build_model()
    plot_model(model)
    result = model.solve()
    result.plot_nodal_field("seqv")
    plt.show()


if __name__ == "__main__":
    main()
