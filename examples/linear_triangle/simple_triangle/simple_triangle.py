# -*- coding: utf-8 -*-
# ***********************************
#  Author: Pedro Jorge De Los Santos
#  E-mail: delossantosmfq@gmail.com
#  License: MIT License
# ***********************************

import matplotlib.pyplot as plt

from nusa import LinearTriangle, LinearTriangleModel, Material, Node, plot_model


def build_model():
    model = LinearTriangleModel("Single CST")

    n1 = Node((0.0, 0.0))
    n2 = Node((1.0, 0.5))
    n3 = Node((0.0, 1.0))

    material = Material(E=200e9, nu=0.3)
    model.add_nodes([n1, n2, n3])
    model.add_element(
        LinearTriangle((n1, n2, n3), material=material, thickness=0.1)
    )
    model.add_constraint(n1, ux=0.0, uy=0.0)
    model.add_constraint(n3, ux=0.0, uy=0.0)
    model.add_force(n2, (1000.0, 0.0))
    return model


def main():
    model = build_model()
    plot_model(model)

    result = model.solve()
    result.plot_nodal_field("ux")
    result.plot_nodal_field("stress_xx")
    plt.show()


if __name__ == "__main__":
    main()
