Mesh utilities
==============

NuSA provides a small 2D geometry and triangular-mesh layer in ``nusa.mesh``.
It is intended as a lightweight preprocessing bridge, not as a replacement for
a general-purpose meshing package.

The current workflow is deliberately explicit:

.. code-block:: text

   geometry -> Modeler -> Gmsh/meshio -> coordinates + connectivity
        -> Node + LinearTriangle -> LinearTriangleModel

The mesher does not create a finite-element ``Model`` automatically. This keeps
preprocessing independent from analysis and makes the conversion from mesh data
to finite-element objects visible to the user.

Installation requirements
-------------------------

``meshio`` is part of the standard NuSA installation, so loading existing
triangular meshes requires no extra Python dependency:

.. code-block:: bash

   pip install nusa

Generating a new mesh from geometry additionally requires the external Gmsh
command-line application. See :doc:`installation` for platform-specific
installation notes.

Before using ``generate_mesh()``, the following command should succeed in the
same terminal, notebook, or environment:

.. code-block:: bash

   gmsh --version

Basic geometry-to-mesh workflow
-------------------------------

.. code-block:: python

   from nusa.mesh import Modeler

   modeler = Modeler()
   outer = modeler.add_rectangle((0.0, 0.0), (1.0, 1.0), esize=0.1)
   hole = modeler.add_circle((0.5, 0.5), 0.15, esize=0.05)
   modeler.subtract_surfaces(outer, hole)

   coordinates, triangles = modeler.generate_mesh()

``coordinates`` is a NumPy array of mesh points. ``triangles`` is a NumPy array
with one zero-based three-node connectivity row per linear triangle.

Only points referenced by triangle cells are returned. Unused points present in
the source mesh are removed and connectivity is remapped to the compacted point
array.

Converting a mesh to a finite-element model
-------------------------------------------

A generated mesh becomes a NuSA model in two explicit steps.

First, convert mesh points to :class:`nusa.node.Node` objects:

.. code-block:: python

   from nusa import LinearTriangle, LinearTriangleModel, Node

   nodes = [Node(tuple(point[:2])) for point in coordinates]

Then convert each connectivity row to a
:class:`nusa.element.LinearTriangle`:

.. code-block:: python

   material = Material(E=200e9, nu=0.3)
   elements = [
       LinearTriangle(
           (nodes[int(i)], nodes[int(j)], nodes[int(k)]),
           material=material,
           thickness=0.01,
       )
       for i, j, k in triangles
   ]

   model = LinearTriangleModel("Meshed plate")
   model.add_nodes(nodes)
   model.add_elements(elements)

Loads and constraints can then be selected using the mesh geometry, after which
the model is solved exactly like a manually constructed model:

.. code-block:: python

   result = model.solve()

For a complete geometry-to-result walkthrough, including boundary selection and
von Mises post-processing, see :doc:`examples/meshed_plate`.

Geometry helpers
----------------

``Modeler`` currently supports:

* ``add_rectangle(p0, p1, esize=0.1)``
* ``add_poly(*points, esize=0.1)``
* ``add_circle(center, radius, esize=0.1)``
* ``subtract_surfaces(outer, inner)``

``esize`` controls the characteristic mesh size attached to generated geometry
points. Smaller values generally produce a finer mesh and therefore more finite
elements.

The circle helper is emitted as four quarter-circle Gmsh arcs so that the
generated curve loop is valid with Gmsh's built-in geometry kernel.

Visualizing the mesh
--------------------

After generating or loading a mesh, ``plot_mesh()`` displays the triangle
connectivity:

.. code-block:: python

   modeler.plot_mesh()

This visualization belongs to preprocessing. It is distinct from
``plot_model(model)``, which visualizes the finite-element problem, and from
``result.plot_nodal_field(...)``, which visualizes solved results.

Loading existing meshes
-----------------------

Triangular meshes supported by ``meshio`` can be loaded without invoking Gmsh:

.. code-block:: python

   modeler = Modeler()
   coordinates, triangles = modeler.generate_mesh_from_file("mesh.msh")

If the file does not contain linear ``triangle`` cells, NuSA raises
``ValueError``.

The same point-compaction behavior used for generated meshes is applied to
loaded meshes.

Selecting a Gmsh executable
---------------------------

The default executable name is ``gmsh``. A custom executable can be supplied
when needed:

.. code-block:: python

   coordinates, triangles = modeler.generate_mesh(
       gmsh_executable="/path/to/gmsh",
       verbose=True,
   )

On Windows, conda installations can expose Gmsh through a ``.bat`` or ``.cmd``
launcher. NuSA resolves these launchers through ``cmd.exe`` automatically.

Generated files
---------------

NuSA writes temporary ``.geo`` and ``.msh`` files for each mesh-generation
call. These files are removed automatically when the operation finishes.

The Gmsh output is requested in the MSH2 format and then read through
``meshio``. The public API exposes only the resulting coordinate and triangle
arrays.

Troubleshooting
---------------

If NuSA reports that Gmsh was not found:

#. run ``gmsh --version`` in the same terminal or notebook environment;
#. if the command is missing, install Gmsh or add its directory to ``PATH``;
#. if Gmsh exists outside ``PATH``, pass ``gmsh_executable`` explicitly.

If Gmsh starts but rejects the generated geometry, call
``generate_mesh(verbose=True)`` to expose Gmsh's command-line diagnostics.

If mesh generation succeeds but the finite-element model is singular, inspect
the applied constraints and verify that the selected boundary nodes actually
lie on the intended geometric edge. For generated coordinates, use
``numpy.isclose`` rather than exact floating-point equality when selecting
boundaries.

API reference
-------------

.. automodule:: nusa.mesh
   :members:
   :undoc-members:
   :show-inheritance:
