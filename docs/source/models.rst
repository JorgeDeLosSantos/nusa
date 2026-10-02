Models
======

Finite-element models define the problem: topology, active degrees of freedom,
applied loads, and prescribed displacements. Solved state is returned as a
separate :class:`nusa.result.StaticResult` snapshot.

Model responsibilities
----------------------

A model owns:

* the ordered node collection;
* the ordered element collection;
* element-family compatibility validation;
* active displacement and force DOFs;
* applied nodal loads;
* prescribed displacements.

A model does **not** own solved displacements, reactions, stresses, strains, or
assembled solver state.

Common workflow
---------------

.. code-block:: python

   model = TrussModel("Example")
   model.add_nodes([n1, n2])
   model.add_element(element)
   model.add_force(n2, (10.0, 0.0))
   model.add_constraint(n1, ux=0.0, uy=0.0)

   result = model.solve()

Nodes and elements
------------------

Use ``add_node`` / ``add_nodes`` and ``add_element`` / ``add_elements`` to build
the topology:

.. code-block:: python

   model.add_nodes([n1, n2, n3])
   model.add_elements([e1, e2])

NuSA assigns labels automatically when they are missing. Public labels may be
strings, sparse integers, or other user-facing identifiers; the solver uses an
internal contiguous ordering independent from those labels.

Loads and constraints
---------------------

The generic model API validates components against the DOFs declared by each
family.

.. list-table:: Model family DOFs
   :header-rows: 1

   * - Model
     - Displacements
     - Forces/actions
   * - ``SpringModel``
     - ``ux``
     - ``fx``
   * - ``BarModel``
     - ``ux``
     - ``fx``
   * - ``TrussModel``
     - ``ux``, ``uy``
     - ``fx``, ``fy``
   * - ``BeamModel``
     - ``uy``, ``ur``
     - ``fy``, ``m``
   * - ``LinearTriangleModel``
     - ``ux``, ``uy``
     - ``fx``, ``fy``

Examples:

.. code-block:: python

   truss.add_force(node, (1000.0, -500.0))
   truss.add_constraint(node, ux=0.0, uy=0.0)

   beam.add_force(node, (-1000.0,))
   beam.add_moment(node, (250.0,))
   beam.add_constraint(node, uy=0.0, ur=0.0)

Prescribed displacements are not limited to zero. For example:

.. code-block:: python

   model.add_constraint(node, ux=0.002)

Inspecting problem inputs
-------------------------

The model exposes the problem inputs before or after a solve:

.. code-block:: python

   model.applied_load(node)
   model.prescribed_displacement(node)
   model.applied_loads
   model.prescribed_displacements

Free displacement DOFs in ``prescribed_displacements`` are represented by
``nan``.

Solving
-------

``model.solve()`` is a convenience wrapper around
:class:`nusa.analysis.LinearStaticAnalysis`:

.. code-block:: python

   result = model.solve()

Calling it repeatedly creates independent results. The model itself remains a
problem definition.

Model families
--------------

The five public model families are intentionally thin, declarative subclasses
of :class:`nusa.model.Model`. They primarily declare the compatible element type
and active DOFs. ``BeamModel`` additionally provides ``add_moment`` as a
family-specific convenience.

API reference
-------------

.. automodule:: nusa.model
   :members:
   :undoc-members:
   :show-inheritance:
