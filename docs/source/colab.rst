Google Colab
============

NuSA can be used in Google Colab for small teaching and experimentation
workflows. Linear-static analyses use the same API as a local installation.
Mesh generation additionally requires the external Gmsh executable in the
current Colab runtime.

Install NuSA
------------

For a published stable release, install NuSA from PyPI:

.. code-block:: bash

   !pip install nusa

While testing the current development line, install ``develop`` directly from
GitHub instead:

.. code-block:: bash

   !pip install "nusa @ git+https://github.com/JorgeDeLosSantos/nusa.git@develop"

Verify the imported version:

.. code-block:: python

   import nusa

   print(nusa.__version__)

Basic solve smoke test
----------------------

The following case checks the standard ``Model -> solve -> StaticResult``
workflow without Gmsh:

.. code-block:: python

   from nusa import Bar, BarModel, Material, Node, Section

   n1 = Node((0.0, 0.0))
   n2 = Node((1.0, 0.0))
   material = Material(E=200e9)
   section = Section(A=1e-4)

   model = BarModel("Colab bar smoke test")
   model.add_nodes([n1, n2])
   model.add_element(Bar((n1, n2), material=material, section=section))
   model.add_constraint(n1, ux=0.0)
   model.add_force(n2, (1000.0,))

   result = model.solve()

   print(result.displacement(n2))
   print(result.reaction(n1))

The result should contain finite displacement and reaction values without
writing solved state back to the nodes.

Install Gmsh
------------

Colab runtimes are Ubuntu-based. Install Gmsh in the active notebook session:

.. code-block:: bash

   !apt-get update -qq
   !apt-get install -y gmsh
   !gmsh --version

The installation is temporary and must be repeated after creating a new Colab
runtime.

Mesh-to-result smoke test
-------------------------

This example exercises the complete preprocessing and solution path:

.. code-block:: text

   geometry -> Gmsh -> triangular mesh -> LinearTriangleModel
            -> solve -> recovered von Mises stress

.. code-block:: python

   import numpy as np

   from nusa import LinearTriangle, LinearTriangleModel, Material, Node
   from nusa.mesh import Modeler

   modeler = Modeler()
   modeler.add_rectangle((0.0, 0.0), (0.2, 0.1), esize=0.05)
   coordinates, connectivity = modeler.generate_mesh()

   nodes = [Node(tuple(point[:2])) for point in coordinates]
   material = Material(E=200e9, nu=0.3)
   elements = [
       LinearTriangle(
           (nodes[int(i)], nodes[int(j)], nodes[int(k)]),
           material=material,
           thickness=0.01,
       )
       for i, j, k in connectivity
   ]

   model = LinearTriangleModel("Colab meshed plate")
   model.add_nodes(nodes)
   model.add_elements(elements)

   xmin = coordinates[:, 0].min()
   xmax = coordinates[:, 0].max()
   loaded = [node for node in nodes if np.isclose(node.x, xmax)]
   force_per_node = 1000.0 / len(loaded)

   for node in nodes:
       if np.isclose(node.x, xmin):
           model.add_constraint(node, ux=0.0, uy=0.0)
       if np.isclose(node.x, xmax):
           model.add_force(node, (force_per_node, 0.0))

   result = model.solve()
   von_mises = result.nodal_field("von_mises_stress")

   assert np.all(np.isfinite(result.displacements))
   assert np.all(np.isfinite(von_mises))

   print(f"nodes: {len(nodes)}")
   print(f"elements: {len(elements)}")
   print(f"max von Mises: {von_mises.max():.6e}")

If this cell completes, the notebook has exercised NuSA, the external Gmsh
executable, ``meshio``, CST assembly and solution, and nodal stress recovery.

Plotting in Colab
-----------------

Matplotlib figures work normally in a notebook. For example:

.. code-block:: python

   import matplotlib.pyplot as plt

   modeler.plot_mesh()
   result.plot_nodal_field("von_mises_stress")
   plt.show()

Troubleshooting
---------------

If ``generate_mesh()`` reports that Gmsh cannot be found, run
``!gmsh --version`` in the same notebook runtime. If that command fails,
repeat the Gmsh installation cell.

If package code changes while a notebook remains open, reinstall the desired
NuSA revision and restart the runtime before interpreting unexpected behavior.
