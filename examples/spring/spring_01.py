# -*- coding: utf-8 -*-
# ***********************************
#  Author: Pedro Jorge De Los Santos
#  E-mail: delossantosmfq@gmail.com
#  License: MIT License
# ***********************************

from nusa import Node, Spring, SpringModel


def test1():
    """Logan (2007), Example 2.1."""
    P = 5000.0
    k1, k2, k3 = 1000.0, 2000.0, 3000.0

    model = SpringModel("2D Model")
    n1, n2, n3, n4 = [Node((0.0, 0.0)) for _ in range(4)]
    e1 = Spring((n1, n3), k1)
    e2 = Spring((n3, n4), k2)
    e3 = Spring((n4, n2), k3)

    model.add_nodes([n1, n2, n3, n4])
    model.add_elements([e1, e2, e3])
    model.add_force(n4, (P,))
    model.add_constraint(n1, ux=0.0)
    model.add_constraint(n2, ux=0.0)

    result = model.solve()

    print("Nodal displacements")
    print("UX3:", result.displacement(n3)["ux"])
    print("UX4:", result.displacement(n4)["ux"])

    print("\nSupport reactions")
    print("R1:", result.reaction(n1)["fx"])
    print("R2:", result.reaction(n2)["fx"])

    print("\nElement forces")
    for element in (e1, e2, e3):
        print(element.label, result.element_result(element))

    return result


if __name__ == "__main__":
    test1()
