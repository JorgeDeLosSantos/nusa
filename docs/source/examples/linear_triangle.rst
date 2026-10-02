Linear triangle elements
========================

The linear triangle is a three-node constant-strain triangle (CST) for 2D plane
stress. Each node contributes ``ux`` and ``uy`` displacement DOFs, with
corresponding ``fx`` and ``fy`` nodal forces.

Use :class:`nusa.model.LinearTriangleModel` with
:class:`nusa.element.LinearTriangle`.

Element properties
------------------

A CST element requires Young's modulus ``E``, Poisson's ratio ``nu``, and
thickness ``t``:

.. code-block:: python

   element = LinearTriangle(
       (n1, n2, n3),
       E=200e9,
       nu=0.3,
       t=0.1,
   )

The canonical element result contains constant strain and stress components:

.. code-block:: text

   strain_xx
   strain_yy
   strain_xy
   stress_xx
   stress_yy
   stress_xy

Example: single CST element
---------------------------

The following example builds a single triangular element, constrains two nodes,
applies an in-plane load, solves the model, and plots nodal fields.

.. literalinclude:: ../../../examples/linear_triangle/simple_triangle/simple_triangle.py
   :language: python
   :linenos:

Element and nodal fields
------------------------

Element fields are read directly from the canonical element results:

.. code-block:: python

   result.element_field("stress_xx")

Nodal fields can be recovered from element values:

.. code-block:: python

   result.nodal_field("stress_xx")

For the current CST implementation, element scalar fields are recovered to
nodes using arithmetic averaging over adjacent elements.

Derived fields
--------------

NuSA also provides selected derived nodal fields:

.. code-block:: python

   result.nodal_field("displacement_magnitude")
   result.nodal_field("von_mises_stress")

Legacy short aliases such as ``sxx`` and ``seqv`` are still accepted by the
post-processing helpers, while the canonical names are preferred in new code.

Visualization
-------------

Problem geometry, loads, and constraints are plotted from the model:

.. code-block:: python

   from nusa import plot_model

   plot_model(model)

Solved scalar fields are plotted from the result:

.. code-block:: python

   result.plot_nodal_field("stress_xx")
   result.plot_element_field("stress_xx")

This distinction is especially useful for continuum problems because the mesh
belongs to the problem definition while stress/strain fields belong to a
particular analysis result.
