# -*- coding: utf-8 -*-
# ***********************************
#  Author: Pedro Jorge De Los Santos
#  E-mail: delossantosmfq@gmail.com
#  License: MIT License
# ***********************************

from nusa import Bar, BarModel, Material, Node, Section


def test1():
    """Logan (2007), Example 3.1."""
    model = BarModel("Bar Model")

    n1 = Node((0.0, 0.0))
    n2 = Node((30.0, 0.0))
    n3 = Node((60.0, 0.0))
    n4 = Node((90.0, 0.0))

    material_30 = Material(E=30e6)
    material_15 = Material(E=15e6)
    section_1 = Section(A=1.0)
    section_2 = Section(A=2.0)

    e1 = Bar((n1, n2), material=material_30, section=section_1)
    e2 = Bar((n2, n3), material=material_30, section=section_1)
    e3 = Bar((n3, n4), material=material_15, section=section_2)

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
