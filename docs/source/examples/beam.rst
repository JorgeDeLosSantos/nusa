Beam elements
=============

The current beam formulation is a two-node Euler-Bernoulli beam. Each node has
transverse displacement ``uy`` and rotation ``ur``. The corresponding generalized
nodal actions are transverse force ``fy`` and bending moment ``m``.

Use :class:`nusa.model.BeamModel` with :class:`nusa.element.Beam`.

Element properties
------------------

A beam element requires Young's modulus ``E`` and second moment of area ``I``:

.. code-block:: python

   e1 = Beam((n1, n2), E=210e9, I=4e-4)

Beam element results are reported as end actions:

.. code-block:: text

   shear_force_i
   shear_force_j
   bending_moment_i
   bending_moment_j

Loads and moments
-----------------

Transverse force and nodal moment are added separately:

.. code-block:: python

   model.add_force(n2, (-10e3,))
   model.add_moment(n2, (20e3,))

Example: two-element beam
-------------------------

The following example contains two beam elements, a nodal force, a nodal
moment, and fixed end conditions.

.. literalinclude:: ../../../examples/beam/beam_2.py
   :language: python
   :linenos:

Important result queries
------------------------

Displacement and rotation:

.. code-block:: python

   values = result.displacement(n2)
   values["uy"]
   values["ur"]

Support actions:

.. code-block:: python

   result.reaction(n1)

Element end actions:

.. code-block:: python

   result.element_result(e1)

Beam diagrams are result-owned visualizations:

.. code-block:: python

   result.plot_shear_diagram()
   result.plot_moment_diagram()

This keeps the beam model as the problem definition while the diagrams are
derived from the frozen element actions stored in ``StaticResult``.
