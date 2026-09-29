# -*- coding: utf-8 -*-
# ***********************************
#  Author: Pedro Jorge De Los Santos
#  E-mail: delossantosmfq@gmail.com
#  License: MIT License
# ***********************************

from nusa import Bar, BarModel, Node


def test1():
    """Logan (2007), Example 3.1."""
    model = BarModel("Bar Model")

    n1 = Node((0.0, 0.0))
    n2 = Node((30.0, 0.0))
    n3 = Node((60.0, 0.0))
    n4 = Node((90.0, 0.0))

    e1 = Bar((n1, n2), E=30e6, A=1.0)
    e2 = Bar((n2, n3), E=30e6, A=1.0)
    e3 = Bar((n3, n4), E=15e6, A=2.0)

    model.add_nodes([n1, n2, n3, n4])
    model.add_elements([e1, e2, e3])
    model.add_force(n2, (3000.0,))
    model.add_constraint(n1, ux=0.0)
    model.add_constraint(n4, ux=0.0)

    result = model.solve()

    print("Node | Displacement | Nodal force")
    for node in model.nodes:
        print(
            f"{node.label}\t"
            f"{result.displacement(node)['ux']:.6f}\t"
            f"{result.nodal_force(node)['fx']}"
        )

    return result


if __name__ == "__main__":
    test1()
