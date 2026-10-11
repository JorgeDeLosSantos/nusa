# -*- coding: utf-8 -*-
# ***********************************
#  Author: Pedro Jorge De Los Santos
#  E-mail: delossantosmfq@gmail.com
#  License: MIT License
# ***********************************

import matplotlib.pyplot as plt

from nusa import Beam, BeamModel, Material, Node, Section


def test4():
    """Beer & Johnston, Mechanics of Materials, Problem 9.13."""
    E = 29e6
    I = 291.0
    P = 35e3
    L1 = 5 * 12.0
    L2 = 10 * 12.0

    material = Material(E=E)
    section = Section(I=I)

    model = BeamModel("Beam Model")
    n1 = Node((0.0, 0.0))
    n2 = Node((L1, 0.0))
    n3 = Node((L1 + L2, 0.0))
    model.add_nodes([n1, n2, n3])
    model.add_elements([Beam((n1, n2), material=material, section=section), Beam((n2, n3), material=material, section=section)])

    model.add_force(n2, (-P,))
    model.add_constraint(n1, uy=0.0)
    model.add_constraint(n3, uy=0.0)

    result = model.solve()

    print(f"Displacement at point C: {result.displacement(n2)['uy']}")
    result.plot_moment_diagram()
    result.plot_shear_diagram()
    plt.show()
    return result


if __name__ == "__main__":
    test4()
