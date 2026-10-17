# -*- coding: utf-8 -*-
# ***********************************
#  Author: Pedro Jorge De Los Santos
#  E-mail: delossantosmfq@gmail.com
#  License: MIT License
# ***********************************

from nusa import Beam, BeamModel, Material, Node, Section


def test3():
    """Kattan, Example 7.1."""
    E = 210e9
    I = 60e-6
    P = 20e3

    material = Material(E=E)
    section = Section(I=I)

    model = BeamModel("Beam Model")
    n1 = Node((0.0, 0.0))
    n2 = Node((2.0, 0.0))
    n3 = Node((4.0, 0.0))
    elements = [Beam((n1, n2), material=material, section=section), Beam((n2, n3), material=material, section=section)]

    model.add_nodes([n1, n2, n3])
    model.add_elements(elements)
    model.add_force(n2, (-P,))
    model.add_constraint(n1, uy=0.0, ur=0.0)
    model.add_constraint(n3, uy=0.0)

    result = model.solve()

    print([result.displacement(node)["uy"] for node in model.nodes])
    return result


if __name__ == "__main__":
    test3()
