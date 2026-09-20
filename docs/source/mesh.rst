Mesh utilities
==============

NuSA provides a small optional 2D geometry and triangular-mesh layer in
``nusa.mesh``. Install it with:

.. code-block:: bash

   pip install "nusa[mesh]"

``meshio`` is used to read mesh files. Generating a mesh from geometry also
requires the external Gmsh executable to be installed and available on
``PATH``.

Basic workflow
--------------

.. code-block:: python

   from nusa.mesh import Modeler

   modeler = Modeler()
   outer = modeler.add_rectangle((0.0, 0.0), (1.0, 1.0), esize=0.1)
   hole = modeler.add_circle((0.5, 0.5), 0.15, esize=0.05)
   modeler.subtract_surfaces(outer, hole)

   nodes, triangles = modeler.generate_mesh()

``nodes`` contains the point coordinates returned by ``meshio`` and
``triangles`` contains zero-based linear-triangle connectivity.

Geometry helpers
----------------

``Modeler`` currently supports:

* ``add_rectangle(p0, p1, esize=0.1)``
* ``add_poly(*points, esize=0.1)``
* ``add_circle(center, radius, esize=0.1)``
* ``subtract_surfaces(outer, inner)``

The circle helper is emitted as four quarter-circle Gmsh arcs so that the
generated curve loop is valid with Gmsh's built-in geometry kernel.

Loading existing meshes
-----------------------

Triangular meshes supported by ``meshio`` can be loaded without invoking
Gmsh:

.. code-block:: python

   modeler = Modeler()
   nodes, triangles = modeler.generate_mesh_from_file("mesh.msh")

If the file does not contain linear ``triangle`` cells, NuSA raises
``ValueError``.

Selecting a Gmsh executable
---------------------------

The default executable name is ``gmsh``. A custom executable can be supplied
when needed:

.. code-block:: python

   nodes, triangles = modeler.generate_mesh(
       gmsh_executable="/path/to/gmsh",
       verbose=True,
   )

NuSA reports a clear ``RuntimeError`` if the executable cannot be found or if
Gmsh fails while generating the mesh.
