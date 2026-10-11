# -*- coding: utf-8 -*-
# ***********************************
#  Author: Pedro Jorge De Los Santos
#  E-mail: delossantosmfq@gmail.com
#  License: MIT License
# ***********************************

from nusa import Beam, BeamModel, Material, Node, Section


def test2():
    """Logan (2007), Example 4.4."""
    E = 210e9
    I = 4e-4
    P = 10e3
    M = 20e3
    L = 3.0

    material = Material(E=E)
    section = Section(I=I)

    model = BeamModel("Beam Model")
    n1 = Node((0.0, 0.0))
    n2 = Node((L, 0.0))
    n3 = Node((2.0 * L, 0.0))

    e1 = Beam((n1, n2), material=material, section=section)
    e2 = Beam((n2, n3), material=material, section=section)

    model.add_nodes([n1, n2, n3])
    model.add_elements([e1, e2])
    model.add_force(n2, (-P,))
    model.add_moment(n2, (M,))
    model.add_constraint(n1, uy=0.0, ur=0.0)
    model.add_constraint(n3, uy=0.0, ur=0.0)

    result = model.solve()

    print("Node 2 displacement:", result.displacement(n2))
    print("Nodal forces:", result.nodal_forces)
    print("Element 1 actions:", result.element_result(e1))
    print("Element 2 actions:", result.element_result(e2))

    return result


if __name__ == "__main__":
    test2()
