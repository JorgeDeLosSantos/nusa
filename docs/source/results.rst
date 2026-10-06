Results
=======

``StaticResult`` is the numerical source of truth after a linear-static solve.
It stores a frozen snapshot of the problem and solution at the time the
analysis was performed.

Creating a result
-----------------

The most common form is:

.. code-block:: python

   result = model.solve()

The equivalent explicit analysis call is:

.. code-block:: python

   from nusa import LinearStaticAnalysis

   result = LinearStaticAnalysis().solve(model)

Both forms return :class:`nusa.result.StaticResult`.

Primary arrays
--------------

The result exposes defensive copies of the main solved vectors:

``applied_loads``
   External nodal loads specified in the model.

``prescribed_displacements``
   Prescribed displacement vector. Free DOFs are represented by ``nan``.

``displacements``
   Solved generalized displacement vector.

``nodal_forces``
   Generalized nodal force vector obtained from ``K @ u``.

``reactions``
   Constraint reactions, evaluated from the difference between the force
   required by the solved state and the applied load.

``element_results``
   Canonical result dictionaries in element order.

Example:

.. code-block:: python

   print(result.displacements)
   print(result.reactions)
   print(result.element_results)

Object-based queries
--------------------

For most user code, object-based accessors are easier to read than indexing the
global arrays manually:

.. code-block:: python

   result.applied_load(node)
   result.prescribed_displacement(node)
   result.displacement(node)
   result.nodal_force(node)
   result.reaction(node)
   result.element_result(element)

The returned dictionaries use the DOF names associated with the model family.
For example, a truss node returns ``ux``/``uy`` displacements and ``fx``/``fy``
forces, while a beam node returns ``uy``/``ur`` displacements and ``fy``/``m``
actions.

Element results
---------------

Each element family exposes a canonical result schema:

.. list-table:: Canonical element results
   :header-rows: 1

   * - Element
     - Result keys
   * - Spring
     - ``force_i``, ``force_j``
   * - Bar
     - ``force_i``, ``force_j``, ``axial_force``, ``axial_stress``
   * - Truss
     - ``axial_force``, ``axial_stress``
   * - Beam
     - ``shear_force_i``, ``shear_force_j``, ``bending_moment_i``, ``bending_moment_j``
   * - LinearTriangle
     - ``stress_xx``, ``stress_yy``, ``stress_xy``, ``strain_xx``, ``strain_yy``, ``strain_xy``

Example:

.. code-block:: python

   values = result.element_result(element)
   print(values["axial_force"])

Snapshot semantics
------------------

A result does not change if the model is later modified:

.. code-block:: python

   result1 = model.solve()

   model.add_force(node, (20.0,))
   result2 = model.solve()

``result1`` still represents the first analysis, while ``result2`` represents
the modified problem.

The result also stores frozen geometry, labels, connectivity, DOF metadata, and
model name. This allows reporting and visualization to remain stable even when
the original model objects are later edited.

Post-processing delegates
-------------------------

``StaticResult`` provides convenience methods that delegate to NuSA's
post-processing functions:

.. code-block:: python

   result.element_field("stress_xx")
   result.nodal_field("stress_xx")
   result.nodal_field("von_mises_stress")

See :doc:`postprocessing` for recovery policies and derived fields.

Reporting and visualization delegates
-------------------------------------

A result can also be used directly for solved reporting and visualization:

.. code-block:: python

   result.simple_report()
   result.plot_deformed_shape()
   result.plot_nodal_field("stress_xx")
   result.plot_element_field("stress_xx")

Beam results additionally expose:

.. code-block:: python

   result.plot_shear_diagram()
   result.plot_moment_diagram()

API reference
-------------

.. automodule:: nusa.result
   :members:
   :undoc-members:
   :show-inheritance:
