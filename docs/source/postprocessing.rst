Post-processing
===============

Post-processing consumes :class:`nusa.result.StaticResult` snapshots rather
than mutable model, node, or element solution state.

Element fields
--------------

Use ``element_field`` to extract one canonical scalar quantity from every
element result:

.. code-block:: python

   from nusa import element_field

   stress_xx = element_field(result, "stress_xx")

The equivalent result convenience method is:

.. code-block:: python

   stress_xx = result.element_field("stress_xx")

The requested name must be present in every element result included in the
field.

Nodal fields
------------

Use ``nodal_field`` to obtain a scalar field at nodes:

.. code-block:: python

   from nusa import nodal_field

   ux = nodal_field(result, "ux")
   stress_xx = nodal_field(result, "stress_xx")

or through the result:

.. code-block:: python

   ux = result.nodal_field("ux")
   stress_xx = result.nodal_field("stress_xx")

Displacement components are read directly from the solved displacement vector.
Element fields that require recovery are mapped to nodes using the selected
recovery policy.

Current recovery policy
-----------------------

NuSA 0.4 currently provides arithmetic averaging for scalar element-to-node
recovery:

.. code-block:: python

   result.nodal_field("stress_xx", recovery="average")

For each node, values from adjacent elements are averaged. This is intentionally
a simple first recovery policy; it should not be interpreted as a higher-order
stress-recovery technique.

Derived fields
--------------

Selected fields are computed from primary result quantities.

Displacement magnitude:

.. code-block:: python

   result.nodal_field("displacement_magnitude")

Plane-stress von Mises stress:

.. code-block:: python

   result.nodal_field("von_mises_stress")

For the current 2D stress convention:

.. math::

   \sigma_{vm} = \sqrt{\sigma_{xx}^2 - \sigma_{xx}\sigma_{yy}
   + \sigma_{yy}^2 + 3\tau_{xy}^2}

Field naming
------------

Canonical field names are preferred in new code:

.. code-block:: text

   stress_xx
   stress_yy
   stress_xy
   strain_xx
   strain_yy
   strain_xy
   displacement_magnitude
   von_mises_stress

Legacy short aliases used by older NuSA examples are still accepted by the
post-processing helpers, including ``sxx``, ``syy``, ``sxy``, ``exx``, ``eyy``,
``exy``, ``usum``, and ``seqv``.

Errors and applicability
------------------------

A field request raises a clear error when the requested quantity is unavailable
for the supplied result or when an unsupported recovery policy is requested.
For example, continuum stress fields are meaningful for
``LinearTriangle`` results but not for a spring-only model.

API reference
-------------

.. automodule:: nusa.post
   :members:
   :undoc-members:
   :show-inheritance:
