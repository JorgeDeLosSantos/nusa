Elements
========

NuSA provides five finite-element formulations in :mod:`nusa.element`.
Elements own connectivity and geometry-dependent formulation and reference
reusable physical-property objects when required. Solved response is evaluated explicitly from the element
displacement vector and stored in analysis results.

Common element contract
-----------------------

Every element has an element type, a label assigned by the model, and a tuple of
connected nodes.

The numerical analysis asks an element for its stiffness matrix and evaluates
its solved response from an explicit local displacement vector.

Conceptually:

.. code-block:: python

   Ke = element.get_element_stiffness()
   values = element.compute_results(u_e)

Elements do not read solved displacement attributes from nodes and do not cache
solved forces, stresses, or strains on themselves.

Spring
------

``Spring`` is a two-node, one-dimensional elastic element.

Required property:

``k``
   Positive spring stiffness.

Local displacement vector:

.. code-block:: text

   [ux_i, ux_j]

Canonical results:

.. code-block:: text

   force_i
   force_j

Bar
---

``Bar`` is a two-node axial element.

Requires ``material=Material(E=...)`` and ``section=Section(A=...)``.
Both references are read-only. The convenience properties ``E`` and ``A``
are delegated to these objects.

Local displacement vector:

.. code-block:: text

   [ux_i, ux_j]

Canonical results:

.. code-block:: text

   force_i
   force_j
   axial_force
   axial_stress

Positive axial force/stress denotes tension.

Truss
-----

``Truss`` is a two-dimensional, two-node axial element whose orientation maps
its axial response into the global coordinate system.

Requires ``material=Material(E=...)`` and ``section=Section(A=...)``.

Local/global element displacement vector used by the current formulation:

.. code-block:: text

   [ux_i, uy_i, ux_j, uy_j]

Canonical results:

.. code-block:: text

   axial_force
   axial_stress

Beam
----

``Beam`` is the current two-node Euler-Bernoulli beam formulation.

Requires ``material=Material(E=...)`` and ``section=Section(I=...)``.

Element displacement vector:

.. code-block:: text

   [uy_i, ur_i, uy_j, ur_j]

Canonical end actions:

.. code-block:: text

   shear_force_i
   shear_force_j
   bending_moment_i
   bending_moment_j

LinearTriangle
--------------

``LinearTriangle`` is a three-node constant-strain triangle (CST) for plane
stress.

Requires ``material=Material(E=..., nu=...)`` and a positive ``thickness``
parameter. ``thickness`` and its alias ``t`` are read-only.

Element displacement vector:

.. code-block:: text

   [ux_1, uy_1, ux_2, uy_2, ux_3, uy_3]

The element additionally exposes explicit strain/stress evaluation:

.. code-block:: python

   strain = element.compute_strain(u_e)
   stress = element.compute_stress(u_e)
   values = element.compute_results(u_e)

Canonical results:

.. code-block:: text

   stress_xx
   stress_yy
   stress_xy
   strain_xx
   strain_yy
   strain_xy

Validation
----------

NuSA validates material and section parameters when their domain objects are
created and checks required properties at element construction. ``thickness``
and spring stiffness must be finite, positive, and non-boolean.
Degenerate triangular geometry is also rejected.

For physical-property definitions see :doc:`materials_sections`.
For worked examples see :doc:`examples/index`.

API reference
-------------

.. automodule:: nusa.element
   :members:
   :undoc-members:
   :show-inheritance:
