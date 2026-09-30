# -*- coding: utf-8 -*-
# ***********************************
#  Author: Pedro Jorge De Los Santos
#  E-mail: delossantosmfq@gmail.com
#  License: MIT License
# ***********************************

from nusa import Beam, BeamModel, Node


def test5():
    """Beer & Johnston, Mechanics of Materials, Problem 9.75."""
    E = 29e6
    b, h = 2.0, 4.0
    I = (1.0 / 12.0) * b * h**3
    w = 1e3 / 12.0
    L1 = 2 * 12.0
    L2 = 3 * 12.0
    P1 = -1e3
    P2 = -w * L2 / 2.0
    P3 = -w * L2 / 2.0
    M2 = -w * L2**2 / 12.0
    M3 = w * L2**2 / 12.0

    model = BeamModel("Beam Model")
    n1 = Node((0.0, 0.0))
    n2 = Node((L1, 0.0))
    n3 = Node((L1 + L2, 0.0))
    model.add_nodes([n1, n2, n3])
    model.add_elements([Beam((n1, n2), E, I), Beam((n2, n3), E, I)])

    model.add_force(n1, (P1,))
    model.add_force(n2, (P2,))
    model.add_force(n3, (P3,))
    model.add_moment(n2, (M2,))
    model.add_moment(n3, (M3,))
    model.add_constraint(n3, uy=0.0, ur=0.0)

    result = model.solve()
    displacement = result.displacement(n1)
    print(
        f"Displacement in node 1: {displacement['uy']}\n"
        f"Slope in node 1: {displacement['ur']}"
    )
    return result


if __name__ == "__main__":
    test5()
