Models
======

Finite-element models define the problem: topology, active degrees of freedom,
applied loads, and prescribed displacements. Solved state is returned as a
separate :class:`nusa.result.StaticResult` snapshot.

The five current public model families are intentionally thin declarations over
the common :class:`nusa.model.Model` behavior.

.. automodule:: nusa.model
   :members:
   :undoc-members:
   :show-inheritance:
