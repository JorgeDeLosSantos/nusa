# -*- coding: utf-8 -*-
# ***********************************
#  Author: Pedro Jorge De Los Santos
#  E-mail: delossantosmfq@gmail.com
#  License: MIT License
# ***********************************

from nusa import Node, Spring, SpringModel


def test2():
    """Logan (2007), Example 2.2."""
    P = 4e3
    k = 200e3

    model = SpringModel("Spring Model 02")
    nodes = [Node((0.0, 0.0)) for _ in range(5)]
    elements = [
        Spring((nodes[k], nodes[k + 1]), k)
        for k in range(4)
    ]

    model.add_nodes(nodes)
    model.add_elements(elements)
    model.add_force(nodes[3], (P,))
    model.add_constraint(nodes[0], ux=0.0)
    model.add_constraint(nodes[4], ux=0.02)

    result = model.solve()

    print(
        "Displacements of nodes 2, 3 and 4:",
        [result.displacement(node)["ux"] for node in nodes[1:4]],
    )
    print(
        "Nodal forces:",
        [result.nodal_force(node)["fx"] for node in nodes],
    )
    print(
        "Element forces:",
        [result.element_result(element) for element in elements],
    )
    return result


if __name__ == "__main__":
    test2()
