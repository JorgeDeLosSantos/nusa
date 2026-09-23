Mesh utilities
==============

NuSA provides a small 2D geometry and triangular-mesh layer in ``nusa.mesh``.
``meshio`` is part of the standard NuSA installation, so loading existing
triangular meshes requires no extra installation:

.. code-block:: bash

   pip install nusa

Generating a new mesh from geometry additionally requires the external Gmsh
command-line application. See :doc:`installation` for platform-specific
installation notes.

Basic workflow
--------------

.. code-block:: python

   from nusa.mesh import Modeler

   modeler = Modeler()
   outer = modeler.add_rectangle((0.0, 0.0), (1.0, 1.0), esize=0.1)
   hole = modeler.add_circle((0.5, 0.5), 0.15, esize=0.05)
   modeler.subtract_surfaces(outer, hole)

   nodes, triangles = modeler.generate_mesh()

``nodes`` contains only points referenced by linear-triangle cells. Any
unused points present in the source mesh are removed, and ``triangles`` is
remapped to zero-based connectivity for the compacted point array.

Before using ``generate_mesh()``, the following command should succeed in the
same environment:

.. code-block:: bash

   gmsh --version

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

On Windows, conda installations can expose Gmsh through a ``.bat`` or ``.cmd``
launcher. NuSA resolves these launchers through ``cmd.exe`` automatically.

Troubleshooting
---------------

If NuSA reports that Gmsh was not found:

#. run ``gmsh --version`` in the same terminal or notebook environment;
#. if the command is missing, install Gmsh or add its directory to ``PATH``;
#. if Gmsh exists outside ``PATH``, pass ``gmsh_executable`` explicitly.

If Gmsh starts but rejects the generated geometry, call
``generate_mesh(verbose=True)`` to expose Gmsh's command-line diagnostics.

NuSA writes temporary ``.geo`` and ``.msh`` files for each generation call and
removes them automatically when the operation finishes.
