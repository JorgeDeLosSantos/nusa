Analysis
========

Linear-static analyses transform a finite-element problem definition into a
stable result snapshot.

Public entry points
-------------------

The simplest form is the model convenience method:

.. code-block:: python

   result = model.solve()

The top-level function is equivalent:

.. code-block:: python

   from nusa import solve

   result = solve(model)

For explicit analysis objects:

.. code-block:: python

   from nusa import LinearStaticAnalysis

   analysis = LinearStaticAnalysis()
   result = analysis.solve(model)

All three forms produce :class:`nusa.result.StaticResult`.

Analysis stages
---------------

A linear-static solve performs the following high-level operations:

.. code-block:: text

   validate model topology
        |
        v
   assemble global stiffness matrix
        |
        v
   build applied-load and prescribed-displacement vectors
        |
        v
   partition free and prescribed DOFs
        |
        v
   solve the reduced linear system
        |
        v
   compute K @ u and reactions
        |
        v
   compute canonical element results
        |
        v
   freeze a StaticResult snapshot

The global stiffness matrix, reduced matrix, right-hand side, and solver
partition are internal analysis details. They are intentionally not stored as
public solved state on ``Model``.

Prescribed displacements
------------------------

The analysis supports both zero and nonzero prescribed displacement values.
For example:

.. code-block:: python

   model.add_constraint(n1, ux=0.0)
   model.add_constraint(n3, ux=0.002)

The contribution of prescribed nonzero DOFs is included in the reduced-system
right-hand side before solving the free DOFs.

Singular systems
----------------

A model that is not sufficiently constrained can produce a singular reduced
system. NuSA checks the reduced system and raises an error rather than silently
returning an invalid solution.

Typical causes include:

* rigid-body modes;
* missing supports;
* disconnected nodes or elements;
* invalid topology.

When this occurs, verify the physical constraints and element connectivity
before changing numerical tolerances.

Fresh solve semantics
---------------------

Each analysis assembles and solves the current problem definition from scratch.
There is no persistent solved-state cache on ``Model``.

This means the following is valid and produces two independent snapshots:

.. code-block:: python

   result1 = model.solve()

   model.add_force(node, (20.0,))
   result2 = model.solve()

``result1`` remains unchanged after the model is modified.

Supported analysis type
-----------------------

NuSA 0.4 currently exposes one public analysis type: linear static analysis.
The separation between ``Model``, ``Analysis``, and ``Result`` is designed so
future analysis types can be added without changing the meaning of nodes,
elements, or the existing static-result contract.

API reference
-------------

.. automodule:: nusa.analysis
   :members:
   :undoc-members:
   :show-inheritance:
