How NuSA works
==============

NuSA 0.4 separates the finite-element problem definition from the numerical
solution. The main public workflow is:

.. code-block:: text

   Node / Element / Model
            |
            v
   LinearStaticAnalysis
            |
            v
       StaticResult
            |
      +-----+------+-------------+
      |            |             |
      v            v             v
   reporting   post-processing  visualization

This separation is intentional. It keeps problem definition, numerical
analysis, and presentation independent from one another.

Problem definition
------------------

A finite-element problem is constructed from three kinds of objects.

Node
~~~~

A :class:`nusa.node.Node` stores geometry and identity. In the current 2D
models a node contains two coordinates and a label assigned by the model.

Nodes do not store solved values. In particular, attributes such as ``ux``,
``uy``, ``fx``, or stress components are not written back to node instances.

Element
~~~~~~~

An :class:`nusa.element.Element` connects nodes and implements a finite-element
formulation. Concrete elements store the physical and geometric parameters
needed to compute their stiffness matrix and response.

Element response is evaluated explicitly from the local displacement vector:

.. code-block:: python

   values = element.compute_results(u_e)

The solved response is then stored in the analysis result rather than on the
element itself.

Model
~~~~~

A :class:`nusa.model.Model` collects the complete finite-element problem:

* nodes and elements;
* topology and internal ordering;
* active displacement and force degrees of freedom;
* applied nodal loads;
* prescribed displacements.

The public model families are thin declarations of these DOF contracts:

.. code-block:: text

   SpringModel          ux       <-> fx
   BarModel             ux       <-> fx
   TrussModel           ux, uy   <-> fx, fy
   BeamModel            uy, ur   <-> fy, m
   LinearTriangleModel  ux, uy   <-> fx, fy

A model remains a problem definition before and after a solve.

Analysis
--------

A linear-static analysis transforms a ``Model`` into a ``StaticResult``.

The convenience form is:

.. code-block:: python

   result = model.solve()

The equivalent explicit form is:

.. code-block:: python

   from nusa import LinearStaticAnalysis

   analysis = LinearStaticAnalysis()
   result = analysis.solve(model)

Internally the analysis performs the usual finite-element steps:

.. code-block:: text

   validate topology
        |
        v
   assemble global stiffness matrix
        |
        v
   build load and prescribed-displacement vectors
        |
        v
   partition constrained/free DOFs
        |
        v
   solve reduced linear system
        |
        v
   compute nodal forces and reactions
        |
        v
   evaluate element results
        |
        v
   create StaticResult snapshot

The assembled stiffness matrix and reduced solver state are implementation
details and are not part of the public solved-state API.

StaticResult
------------

:class:`nusa.result.StaticResult` is the numerical source of truth after a
solve. It stores a frozen snapshot of the problem and solution at the time the
analysis was performed.

Primary solved quantities include:

.. code-block:: python

   result.applied_loads
   result.prescribed_displacements
   result.displacements
   result.nodal_forces
   result.reactions
   result.element_results

Convenience queries accept the original model objects:

.. code-block:: python

   result.applied_load(node)
   result.prescribed_displacement(node)
   result.displacement(node)
   result.nodal_force(node)
   result.reaction(node)
   result.element_result(element)

Snapshot semantics
~~~~~~~~~~~~~~~~~~

Results are intentionally independent from later model mutation:

.. code-block:: python

   result1 = model.solve()

   model.add_force(node, (20.0,))
   result2 = model.solve()

``result1`` still represents the first analysis, while ``result2`` represents
the modified problem.

This is different from designs where solved values are written into mutable
nodes or elements.

Applied loads, nodal forces, and reactions
------------------------------------------

These three quantities represent different concepts and should not be used
interchangeably.

Applied loads
~~~~~~~~~~~~~

``result.applied_loads`` contains the external nodal loads explicitly defined
in the model.

Nodal forces
~~~~~~~~~~~~

``result.nodal_forces`` contains the generalized nodal force vector obtained
from the solved displacement vector:

.. math::

   \mathbf{f}_{\mathrm{nodal}} = \mathbf{K}\mathbf{u}

Reactions
~~~~~~~~~

``result.reactions`` contains the constraint reactions. Conceptually, NuSA
computes them from the difference between the force required by the solved
state and the externally applied load:

.. math::

   \mathbf{r} = \mathbf{K}\mathbf{u} - \mathbf{f}_{\mathrm{applied}}

At unconstrained DOFs these values are normally zero up to numerical precision.

Post-processing
---------------

Post-processing consumes a completed ``StaticResult``. It does not inspect
mutable solved state on nodes or elements.

For example:

.. code-block:: python

   result.element_field("stress_xx")
   result.nodal_field("stress_xx")
   result.nodal_field("von_mises_stress")

For constant-strain triangles, element fields can be recovered to nodes using
the current arithmetic-average recovery policy.

Visualization and reporting
---------------------------

NuSA distinguishes visualization of the problem from visualization of the
solution.

Problem visualization:

.. code-block:: python

   from nusa import plot_model

   plot_model(model)

Solved visualization:

.. code-block:: python

   result.plot_deformed_shape()
   result.plot_nodal_field("stress_xx")
   result.plot_element_field("stress_xx")

Reporting follows the same rule:

.. code-block:: python

   result.simple_report()

Why this architecture?
----------------------

The separation provides several practical benefits:

* a model can be solved repeatedly after changing loads or constraints;
* previous results remain stable snapshots;
* reporting and plotting do not mutate or depend on model solution state;
* element formulations can be tested directly using explicit displacement
  vectors;
* future analysis types can return their own result objects without changing
  the meaning of ``Node`` or ``Model``.

The important mental model for users is therefore simple:

.. code-block:: text

   define the problem -> solve it -> work with the result

Continue with :doc:`getting_started` for a complete executable example or the
API reference for detailed class and function documentation.
