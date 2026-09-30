# -*- coding: utf-8 -*-
# ***********************************
#  Author: Pedro Jorge De Los Santos
#  E-mail: delossantosmfq@gmail.com
#  License: MIT License
# ***********************************

from nusa import Node, Spring, SpringModel


def test3():
    """Logan (2007), Problem 2.8."""
    P = 500.0
    k = 500.0

    model = SpringModel("Spring Model 03")
    n1 = Node((0.0, 0.0))
    n2 = Node((0.0, 0.0))
    n3 = Node((0.0, 0.0))
    model.add_nodes([n1, n2, n3])
    model.add_elements([Spring((n1, n2), k), Spring((n2, n3), k)])
    model.add_force(n3, (P,))
    model.add_constraint(n1, ux=0.0)

    result = model.solve()

    for node in model.nodes:
        print(result.displacement(node)["ux"])
    return result


if __name__ == "__main__":
    test3()
