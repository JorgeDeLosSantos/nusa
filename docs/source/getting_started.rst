Getting started
===============

This guide walks through a complete NuSA analysis using a one-dimensional bar
model. The same workflow is used by the spring, truss, beam, and linear
triangle families.

The basic workflow is:

.. code-block:: text

   Node + Element + Model
            |
            v
          solve()
            |
            v
       StaticResult

A model describes the finite-element problem. Solving the model returns a
separate result object containing the numerical solution.

Installation
------------

Install the current stable release from PyPI:

.. code-block:: bash

   pip install nusa

For the current development version:

.. code-block:: bash

   pip install "nusa @ git+https://github.com/JorgeDeLosSantos/nusa.git@develop"

See :doc:`installation` for Gmsh setup, development installation, and platform
specific notes.

A first bar analysis
--------------------

Consider a bar assembled from three finite elements. The end nodes are fixed
and a horizontal force is applied at the second node.

The complete executable example used by NuSA's test suite is shown below:

.. literalinclude:: ../../examples/bar/bar_1.py
   :language: python
   :linenos:

The important steps are discussed below.

1. Create nodes
~~~~~~~~~~~~~~~

Nodes define geometry and identity:

.. code-block:: python

   from nusa import Node

   n1 = Node((0.0, 0.0))
   n2 = Node((30.0, 0.0))

A ``Node`` does not store solved displacements, forces, stresses, or strains.
Those quantities belong to the result of an analysis.

2. Create elements
~~~~~~~~~~~~~~~~~~

Elements connect nodes and define the finite-element formulation. Reusable
material and section objects provide their physical properties:

.. code-block:: python

   from nusa import Bar, Material, Section

   material = Material(E=30e6)
   section = Section(A=1.0)

   e1 = Bar(
       (n1, n2),
       material=material,
       section=section,
   )

For a bar element, the material provides Young's modulus ``E`` and the
section provides cross-sectional area ``A``. Both objects can be reused by
multiple elements.

3. Build the model
~~~~~~~~~~~~~~~~~~

A model contains the finite-element problem definition:

.. code-block:: python

   from nusa import BarModel

   model = BarModel("Bar Model")
   model.add_nodes([n1, n2])
   model.add_element(e1)

The model owns topology, loads, and prescribed displacements. It does not own
the solved state.

4. Apply loads and constraints
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Loads and constraints are expressed using the active degrees of freedom for the
model family:

.. code-block:: python

   model.add_force(n2, (3000.0,))
   model.add_constraint(n1, ux=0.0)

For a ``BarModel`` the only displacement degree of freedom is ``ux`` and the
corresponding force component is ``fx``.

5. Solve the problem
~~~~~~~~~~~~~~~~~~~~

The simplest public entry point is:

.. code-block:: python

   result = model.solve()

This is equivalent to:

.. code-block:: python

   from nusa import LinearStaticAnalysis

   result = LinearStaticAnalysis().solve(model)

Both forms return a :class:`nusa.result.StaticResult`.

6. Read results
~~~~~~~~~~~~~~~

Use the result object to query solved quantities:

.. code-block:: python

   result.displacement(n2)
   result.reaction(n1)
   result.element_result(e1)

The corresponding global arrays are also available:

.. code-block:: python

   result.displacements
   result.applied_loads
   result.nodal_forces
   result.reactions
   result.element_results

The distinction between these quantities is important:

``applied_loads``
   Loads explicitly specified in the model.

``nodal_forces``
   Generalized nodal forces computed as ``K @ u``.

``reactions``
   Forces associated with prescribed degrees of freedom.

7. Report and visualize
~~~~~~~~~~~~~~~~~~~~~~~

Solved reporting and visualization consume the result object:

.. code-block:: python

   result.simple_report()
   result.plot_deformed_shape()

Problem visualization consumes the model instead:

.. code-block:: python

   from nusa import plot_model

   plot_model(model)

This distinction reflects the central NuSA 0.4 design: models define problems;
results describe solved analyses.

Next steps
----------

Continue with :doc:`how_nusa_works` for the architecture and ownership model,
then explore the element families and result API in the user guide and API
reference.
