# -*- coding: utf-8 -*-
# ***********************************
#  Author: Pedro Jorge De Los Santos
#  E-mail: delossantosmfq@gmail.com
#  License: MIT License
# ***********************************

import matplotlib.pyplot as plt
import numpy as np

from nusa import LinearTriangle, LinearTriangleModel, Material, Node, plot_model
from nusa.mesh import Modeler


E = 200e9
NU = 0.3
THICKNESS = 0.01
TOTAL_FORCE = 6000.0


def build_model(esize=0.05, gmsh_executable="gmsh"):
    """Build a meshed square plate using Gmsh-generated CST elements."""
    modeler = Modeler()
    modeler.add_rectangle((0.0, 0.0), (0.3, 0.3), esize=esize)
    coordinates, connectivity = modeler.generate_mesh(
        gmsh_executable=gmsh_executable,
    )

    material = Material(E=E, nu=NU)
    nodes = [Node(tuple(point[:2])) for point in coordinates]
    elements = [
        LinearTriangle(
            (nodes[int(i)], nodes[int(j)], nodes[int(k)]),
            material=material,
            thickness=THICKNESS,
        )
        for i, j, k in connectivity
    ]

    model = LinearTriangleModel("Meshed square plate")
    model.add_nodes(nodes)
    model.add_elements(elements)

    xmin = coordinates[:, 0].min()
    xmax = coordinates[:, 0].max()
    loaded_nodes = [node for node in nodes if np.isclose(node.x, xmax)]
    nodal_force = TOTAL_FORCE / len(loaded_nodes)

    for node in nodes:
        if np.isclose(node.x, xmin):
            model.add_constraint(node, ux=0.0, uy=0.0)
        if np.isclose(node.x, xmax):
            model.add_force(node, (nodal_force, 0.0))

    return model, modeler


def main():
    model, modeler = build_model()

    modeler.plot_mesh()
    plot_model(model)

    result = model.solve()
    result.plot_nodal_field("von_mises_stress")
    plt.show()
    return result


if __name__ == "__main__":
    main()
