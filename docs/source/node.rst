Nodes
=====

Nodes are lightweight geometry and identity objects. Solved displacements,
forces, stresses, and strains are stored in analysis results rather than written
back to nodes.

Creating nodes
--------------

Create a two-dimensional node from ``(x, y)`` coordinates:

.. code-block:: python

   from nusa import Node

   node = Node((1.5, 2.0))

Coordinates are validated as finite two-component values and stored as a NumPy
array.

Geometry access
---------------

Use the convenience properties or the coordinate array:

.. code-block:: python

   node.x
   node.y
   node.coordinates

A node can be moved by editing its coordinate array before running a new
analysis. Existing ``StaticResult`` snapshots keep their own frozen coordinates.

Labels
------

``Node.label`` is a public identifier. A model assigns an unused label
automatically when a node is added without one, but user-supplied labels are
also supported.

For example:

.. code-block:: python

   node = Node((0.0, 0.0))
   node.label = "support-A"
   model.add_node(node)

The numerical solver does not use public labels as array indices. NuSA keeps an
internal contiguous node ordering, so sparse integers or string labels do not
change the assembly layout.

No solved state on Node
-----------------------

NuSA 0.4 deliberately does not create attributes such as:

.. code-block:: text

   node.ux
   node.uy
   node.fx
   node.fy

Read solved values from the analysis result instead:

.. code-block:: python

   result = model.solve()

   result.displacement(node)
   result.nodal_force(node)
   result.reaction(node)

This allows one node/model definition to participate in multiple analyses while
previous results remain independent snapshots.

API reference
-------------

.. automodule:: nusa.node
   :members:
   :undoc-members:
   :show-inheritance:
