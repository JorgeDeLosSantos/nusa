# -*- coding: utf-8 -*-
# ***********************************
#  Author: Pedro Jorge De Los Santos
#  E-mail: delossantosmfq@gmail.com
#  License: MIT License
# ***********************************

from nusa import Bar, BarModel, Material, Node, Section


def test2():
    """Kattan (2003), Example 3.1."""
    E = 210e6
    A = 0.003
    P = -10.0
    UX3 = 0.002

    material = Material(E=E)
    section = Section(A=A)

    model = BarModel("Bar model 02")
    n1 = Node((0.0, 0.0))
    n2 = Node((1.5, 0.0))
    n3 = Node((2.5, 0.0))

    e1 = Bar((n1, n2), material=material, section=section)
    e2 = Bar((n2, n3), material=material, section=section)

    model.add_nodes([n1, n2, n3])
    model.add_elements([e1, e2])
    model.add_force(n2, (P,))
    model.add_constraint(n1, ux=0.0)
    model.add_constraint(n3, ux=UX3)

    result = model.solve()

    print(f"Displacement in node 2: {result.displacement(n2)['ux']:0.6f}")
    print(
        "Reactions at nodes 1 and 3: "
        f"R1={result.reaction(n1)['fx']:0.2f} "
        f"R3={result.reaction(n3)['fx']:0.2f}"
    )
    print(
        "Axial stress in each bar:\n"
        f"Element 1: {result.element_result(e1)['axial_stress']}\n"
        f"Element 2: {result.element_result(e2)['axial_stress']}"
    )
    return result


if __name__ == "__main__":
    test2()
