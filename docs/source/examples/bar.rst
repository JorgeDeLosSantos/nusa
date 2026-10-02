Bar elements
============

A bar element represents one-dimensional axial deformation. Each node contributes
one displacement degree of freedom, ``ux``, with corresponding nodal force
``fx``.

Use :class:`nusa.model.BarModel` with :class:`nusa.element.Bar`.

Element properties
------------------

A bar requires Young's modulus ``E`` and cross-sectional area ``A``:

.. code-block:: python

   e1 = Bar((n1, n2), E=30e6, A=1.0)

The element result contains end forces together with the physical axial force
and axial stress:

.. code-block:: python

   values = result.element_result(e1)
   values["axial_force"]
   values["axial_stress"]

NuSA uses positive axial force/stress for tension.

Example: three-bar assemblage
-----------------------------

The example below assembles three axial bars with different material/section
properties, applies a load at an internal node, and restrains the two ends.

.. literalinclude:: ../../../examples/bar/bar_1.py
   :language: python
   :linenos:

Important result queries
------------------------

Displacement at a node:

.. code-block:: python

   result.displacement(n2)["ux"]

Generalized nodal force:

.. code-block:: python

   result.nodal_force(n2)["fx"]

Element axial response:

.. code-block:: python

   result.element_result(e1)["axial_force"]
   result.element_result(e1)["axial_stress"]

For problems with prescribed nonzero axial displacement, use the same
``add_constraint(node, ux=value)`` interface used for zero supports.
