# -*- coding: utf-8 -*-
# ***********************************
#  Author: Pedro Jorge De Los Santos
#  E-mail: delossantosmfq@gmail.com
#  License: MIT License
# ***********************************

from nusa import Node, Spring, SpringModel


def simple_case():
    model = SpringModel("Simple")
    n1 = Node((0.0, 0.0))
    n2 = Node((0.0, 0.0))
    model.add_nodes([n1, n2])
    model.add_element(Spring((n1, n2), 300.0))
    model.add_force(n2, (750.0,))
    model.add_constraint(n1, ux=0.0)

    result = model.solve()
    result.simple_report()
    return result


if __name__ == "__main__":
    simple_case()
