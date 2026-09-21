# -*- coding: utf-8 -*-
# ***********************************
#  Author: Pedro Jorge De Los Santos
#  E-mail: delossantosmfq@gmail.com
#  Blog: numython.github.io
#  License: MIT License
# ***********************************
from nusa import LinearTriangle, LinearTriangleModel, Node


def build_model():
    model = LinearTriangleModel("Single CST")

    n1 = Node((0, 0))
    n2 = Node((1, 0.5))
    n3 = Node((0, 1))
    element = LinearTriangle((n1, n2, n3), 200e9, 0.3, 0.1)

    model.add_nodes([n1, n2, n3])
    model.add_element(element)
    model.add_constraint(n1, ux=0, uy=0)
    model.add_constraint(n3, ux=0, uy=0)
    model.add_force(n2, (1000, 0))
    model.solve()
    return model


def main():
    model = build_model()
    model.plot_model()
    model.plot_nodal_result("ux")
    model.plot_nodal_result("sxx")
    model.show()


if __name__ == "__main__":
    main()
