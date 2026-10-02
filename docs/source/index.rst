NuSA documentation
==================

NuSA is a small Python finite-element library for teaching, experimentation,
and structural-analysis workflows.

The public workflow in NuSA 0.4 is intentionally explicit:

.. code-block:: text

   define a Model -> solve the analysis -> work with a StaticResult

If you are new to NuSA, start with :doc:`getting_started`. If you already know
the basic API and want to understand the design and ownership of solved state,
read :doc:`how_nusa_works`.

Quick example
-------------

.. code-block:: python

   from nusa import Bar, BarModel, Node

   n1 = Node((0.0, 0.0))
   n2 = Node((1.0, 0.0))

   model = BarModel("Simple bar")
   model.add_nodes([n1, n2])
   model.add_element(Bar((n1, n2), E=200e9, A=1e-4))
   model.add_constraint(n1, ux=0.0)
   model.add_force(n2, (1000.0,))

   result = model.solve()

   print(result.displacement(n2))
   print(result.reaction(n1))

Getting started
---------------

.. toctree::
   :maxdepth: 2
   :caption: Getting started

   installation
   getting_started
   how_nusa_works

Examples
--------

.. toctree::
   :maxdepth: 2
   :caption: Examples

   examples/index

Migration
---------

.. toctree::
   :maxdepth: 1
   :caption: Migration

   migration_0_3_to_0_4

API reference
-------------

.. toctree::
   :maxdepth: 2
   :caption: API reference

   node
   elements
   models
   analysis
   results
   postprocessing
   reporting
   visualization
   mesh

Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
