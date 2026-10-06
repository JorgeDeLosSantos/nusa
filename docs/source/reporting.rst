Reporting
=========

Text reports are generated from completed :class:`nusa.result.StaticResult`
snapshots. Reporting therefore reflects one specific solved analysis rather
than mutable state on a model.

Quick usage
-----------

Print a report:

.. code-block:: python

   result = model.solve()
   result.simple_report()

Return the report as a string:

.. code-block:: python

   report = result.simple_report(report_type="string")
   print(report)

Write the report to a file:

.. code-block:: python

   result.simple_report(
       report_type="write",
       fname="analysis.txt",
   )

The top-level function is equivalent:

.. code-block:: python

   from nusa import simple_report

   report = simple_report(result, report_type="string")

Report contents
---------------

The simple report includes the major result categories available in
``StaticResult``:

.. code-block:: text

   NODAL DISPLACEMENTS
   APPLIED LOADS
   NODAL FORCES (K @ U)
   REACTIONS
   ELEMENT RESULTS
   FINITE ELEMENT MODEL INFO
   NODES
   ELEMENTS

Element-result columns are generated from the canonical element result names.
For example, a bar report includes axial force/stress while a beam report
includes shear forces and bending moments.

Why reporting is result-owned
-----------------------------

In NuSA 0.4, a model can be changed and solved again while old results remain
valid snapshots:

.. code-block:: python

   result1 = model.solve()
   report1 = result1.simple_report(report_type="string")

   model.add_force(node, (20.0,))
   result2 = model.solve()

``report1`` continues to describe ``result1``. It is not affected by the later
model mutation.

This behavior is one of the reasons reporting consumes ``StaticResult`` rather
than reading solved values from ``Node`` or ``Model``.

API reference
-------------

.. automodule:: nusa.reporting
   :members:
   :undoc-members:
   :show-inheritance:
