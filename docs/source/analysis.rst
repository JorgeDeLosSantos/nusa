Analysis
========

Linear-static analyses transform a finite-element problem definition into a
stable result snapshot.

The convenience function::

   result = nusa.solve(model)

is equivalent to::

   result = nusa.LinearStaticAnalysis().solve(model)

.. automodule:: nusa.analysis
   :members:
   :undoc-members:
   :show-inheritance:
