Migrating from 0.3 to 0.4
=========================

NuSA 0.4 intentionally changes how solved state is represented. The central
change is simple:

.. code-block:: text

   NuSA 0.3: solve and read mutable state from Model / Node / Element
   NuSA 0.4: solve and read a stable StaticResult snapshot

The migration is therefore mostly a matter of keeping the object returned by
``solve()`` and reading solved quantities from that result.

Keep the solve result
---------------------

0.3-style code often ignored the return value:

.. code-block:: python

   model.solve()

In 0.4, keep it explicitly:

.. code-block:: python

   result = model.solve()

Node displacements
------------------

0.3:

.. code-block:: python

   model.solve()
   print(node.ux)
   print(node.uy)

0.4:

.. code-block:: python

   result = model.solve()
   displacement = result.displacement(node)
   print(displacement["ux"])
   print(displacement["uy"])

``Node`` no longer stores solved displacement attributes.

Nodal forces and reactions
--------------------------

0.3 code could read force-like values from mutable nodes or model state.

0.4 distinguishes applied loads, generalized nodal forces, and reactions:

.. code-block:: python

   result.applied_load(node)
   result.nodal_force(node)
   result.reaction(node)

Use ``reaction(node)`` for support reactions.

Element results
---------------

0.3:

.. code-block:: python

   print(element.sx)
   print(element.fx)

0.4:

.. code-block:: python

   values = result.element_result(element)
   print(values["axial_stress"])
   print(values["axial_force"])

Canonical result names depend on the element family. See :doc:`results` and
:doc:`elements` for the complete schemas.

Reports
-------

0.3:

.. code-block:: python

   model.simple_report()

0.4:

.. code-block:: python

   result = model.solve()
   result.simple_report()

Reports consume the solved snapshot rather than mutable model state.

Visualization
-------------

Problem visualization and solved visualization are now separated.

0.3-style solved plotting:

.. code-block:: python

   model.plot_deformed_shape()
   model.plot_nodal_result("sxx")

0.4:

.. code-block:: python

   from nusa import plot_model

   plot_model(model)                       # problem definition

   result = model.solve()
   result.plot_deformed_shape()            # solved geometry
   result.plot_nodal_field("stress_xx")   # solved field

Beam diagrams likewise move to the result:

.. code-block:: python

   result.plot_shear_diagram()
   result.plot_moment_diagram()

Global stiffness matrix and assembly state
------------------------------------------

The following 0.3-style public solved/assembly interfaces are no longer part of
the 0.4 public API:

.. code-block:: text

   model.assemble()
   model.stiffness_matrix
   model.displacements
   model.nodal_forces
   model.reactions

Assembly and reduced-system state are implementation details of the analysis.
The public solved interface is ``StaticResult``.

Module layout
-------------

The old catch-all ``nusa.core`` module has been removed.

Use the top-level public imports when possible:

.. code-block:: python

   from nusa import Node, Element, Model

The underlying modules are now:

.. code-block:: text

   nusa.node      -> Node
   nusa.element   -> Element and concrete elements
   nusa.model     -> Model and concrete model families

Typical before/after example
----------------------------

0.3:

.. code-block:: python

   model.solve()

   print(node.ux)
   print(model.reaction(node))
   model.simple_report()
   model.plot_deformed_shape()

0.4:

.. code-block:: python

   result = model.solve()

   print(result.displacement(node)["ux"])
   print(result.reaction(node))
   result.simple_report()
   result.plot_deformed_shape()

Summary table
-------------

.. list-table:: Common 0.3 to 0.4 replacements
   :header-rows: 1

   * - 0.3 pattern
     - 0.4 replacement
   * - ``node.ux``, ``node.uy``, ...
     - ``result.displacement(node)``
   * - node/model solved forces
     - ``result.nodal_force(node)`` or ``result.reaction(node)``
   * - ``element.sx``, ``element.fx``, ...
     - ``result.element_result(element)``
   * - ``model.simple_report()``
     - ``result.simple_report()``
   * - ``model.plot_deformed_shape()``
     - ``result.plot_deformed_shape()``
   * - ``model.plot_nodal_result(...)``
     - ``result.plot_nodal_field(...)``
   * - ``model.plot_element_result(...)``
     - ``result.plot_element_field(...)``
   * - ``model.assemble()``
     - internal analysis operation; no direct public replacement
   * - ``model.stiffness_matrix``
     - internal analysis data; no direct public replacement
   * - ``nusa.core``
     - ``nusa.node``, ``nusa.element``, ``nusa.model`` or top-level imports

For new code, start with :doc:`getting_started` and use
:doc:`how_nusa_works` as the conceptual reference for ownership of problem and
solution state.
