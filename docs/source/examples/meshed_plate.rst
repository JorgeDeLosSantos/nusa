Meshed plate with Gmsh
======================

This example connects NuSA's preprocessing and analysis layers in one complete
workflow:

.. code-block:: text

   geometry -> Gmsh -> triangle mesh -> LinearTriangleModel
       -> solve -> StaticResult -> recovered stress field

It is the recommended reference when a problem contains more triangles than is
practical to create manually.

Problem
-------

The model is a square plate with side length ``0.3`` and thickness ``0.01``.
The left edge is fixed in both translational degrees of freedom and a total
horizontal load of ``6000`` is distributed equally across the nodes on the
right edge.

The material parameters are:

* Young's modulus ``E = 200e9``
* Poisson ratio ``nu = 0.3``
* thickness ``t = 0.01``

The plate is discretized with constant-strain triangular elements.

1. Create the geometry
----------------------

``Modeler`` stores a small 2D geometry description and delegates mesh
generation to the external Gmsh executable:

.. code-block:: python

   from nusa.mesh import Modeler

   modeler = Modeler()
   modeler.add_rectangle((0.0, 0.0), (0.3, 0.3), esize=0.05)
   coordinates, connectivity = modeler.generate_mesh()

``coordinates`` contains the compacted mesh points and ``connectivity`` contains
zero-based three-node triangle indices.

For example, one connectivity row

.. code-block:: text

   [4, 7, 3]

means that one finite element is formed by mesh points ``4``, ``7`` and ``3``.

2. Convert the mesh to NuSA objects
-----------------------------------

Mesh points become :class:`nusa.node.Node` objects:

.. code-block:: python

   nodes = [Node(tuple(point[:2])) for point in coordinates]

Each triangle is then converted to a :class:`nusa.element.LinearTriangle`:

.. code-block:: python

   elements = [
       LinearTriangle(
           (nodes[int(i)], nodes[int(j)], nodes[int(k)]),
           E=200e9,
           nu=0.3,
           t=0.01,
       )
       for i, j, k in connectivity
   ]

The mesher deliberately does not create a finite-element ``Model`` itself.
This keeps preprocessing independent from analysis and makes the conversion
step explicit.

3. Build the finite-element problem
-----------------------------------

.. code-block:: python

   model = LinearTriangleModel("Meshed square plate")
   model.add_nodes(nodes)
   model.add_elements(elements)

Boundary conditions and loads can be selected from geometry. Here, nodes on
the minimum x-coordinate are fixed and nodes on the maximum x-coordinate share
the total horizontal force:

.. code-block:: python

   xmin = coordinates[:, 0].min()
   xmax = coordinates[:, 0].max()

   loaded_nodes = [node for node in nodes if np.isclose(node.x, xmax)]
   force_per_node = 6000.0 / len(loaded_nodes)

   for node in nodes:
       if np.isclose(node.x, xmin):
           model.add_constraint(node, ux=0.0, uy=0.0)
       if np.isclose(node.x, xmax):
           model.add_force(node, (force_per_node, 0.0))

Using ``np.isclose`` is preferable to exact floating-point comparisons when
selecting mesh boundaries.

4. Solve
--------

The solve stage is unchanged from manually created models:

.. code-block:: python

   result = model.solve()

The mesh is now frozen into the returned :class:`nusa.result.StaticResult`
through its node coordinates and connectivity.

5. Post-process a stress field
------------------------------

A CST element produces element stress and strain components. NuSA can recover
those values to nodes and compute the plane-stress von Mises field:

.. code-block:: python

   von_mises = result.nodal_field("von_mises_stress")

or plot it directly:

.. code-block:: python

   result.plot_nodal_field("von_mises_stress")

The current nodal recovery policy is the arithmetic average of adjacent
element values. See :doc:`../postprocessing` for details.

Visualizing preprocessing and results
-------------------------------------

The three visualization stages remain separate:

.. code-block:: python

   modeler.plot_mesh()                  # generated mesh
   plot_model(model)                    # FE problem definition
   result.plot_nodal_field(
       "von_mises_stress"
   )                                    # solved field

This mirrors the ownership model used throughout NuSA:

.. code-block:: text

   Modeler -> mesh
   Model   -> finite-element problem
   Result  -> solved state

Complete executable example
---------------------------

The following script is the same meshed-plate example maintained under
``examples/``:

.. literalinclude:: ../../../examples/linear_triangle/plate/plate.py
   :language: python
   :linenos:

Gmsh requirement
----------------

Generating a new mesh requires the Gmsh command-line executable. Before running
this example, the following should succeed in the same terminal or notebook
environment:

.. code-block:: bash

   gmsh --version

Loading an existing triangular mesh with ``generate_mesh_from_file()`` only
requires ``meshio`` and does not invoke Gmsh.

See :doc:`../mesh` for installation notes, geometry helpers, custom executable
selection, and troubleshooting.
