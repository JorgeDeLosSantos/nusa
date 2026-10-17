Truss elements
==============

A two-dimensional truss element carries axial load while its orientation maps
that axial behavior into global ``x`` and ``y`` directions. Each node contributes
``ux`` and ``uy`` displacement DOFs, with corresponding ``fx`` and ``fy``
forces.

Use :class:`nusa.model.TrussModel` with :class:`nusa.element.Truss`.

Element properties
------------------

A truss receives reusable material and section objects. Its formulation requires
Young's modulus ``E`` from the material and cross-sectional area ``A`` from
the section:

.. code-block:: python

   material = Material(E=30e6)
   section = Section(A=2.0)

   element = Truss(
       (n1, n2),
       material=material,
       section=section,
   )

Its canonical result contains:

.. code-block:: text

   axial_force
   axial_stress

Example: three-member truss
---------------------------

The following example builds a three-member planar truss, applies a vertical
load, visualizes the problem definition, solves it, and plots the deformed
shape.

.. literalinclude:: ../../../examples/truss/truss_01.py
   :language: python
   :linenos:

Problem and solved visualization
--------------------------------

The undeformed problem definition is visualized from the model:

.. code-block:: python

   from nusa import plot_model

   plot_model(model)

The solved shape is visualized from the result:

.. code-block:: python

   result = model.solve()
   result.plot_deformed_shape()

Important result queries
------------------------

Nodal displacement:

.. code-block:: python

   result.displacement(n1)

Support reaction:

.. code-block:: python

   result.reaction(n2)

Axial member response:

.. code-block:: python

   result.element_result(model.elements[0])["axial_force"]

The model stores geometry and supports; the ``StaticResult`` stores the solved
state used for deformation plots and member forces.
