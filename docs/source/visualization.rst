Visualization
=============

NuSA 0.4 separates visualization of the problem definition from visualization
of solved results.

Problem visualization
---------------------

Use :func:`nusa.visualization.plot_model` to visualize geometry, applied
translational loads, and prescribed constraints:

.. code-block:: python

   from nusa import plot_model

   ax = plot_model(model)

This function consumes a ``Model`` because it represents the finite-element
problem before or after any analysis is performed.

Solved geometry
---------------

Use a completed ``StaticResult`` to plot the deformed shape:

.. code-block:: python

   result = model.solve()
   ax = result.plot_deformed_shape(scale=10.0)

The top-level function is equivalent:

.. code-block:: python

   from nusa import plot_deformed_shape

   ax = plot_deformed_shape(result, scale=10.0)

The plot uses the frozen geometry and displacement vector stored in the result,
so it remains stable if the original model is later modified.

Scalar fields
-------------

Triangular continuum results can be plotted as nodal or element scalar fields.

Nodal field:

.. code-block:: python

   result.plot_nodal_field("stress_xx")

Element field:

.. code-block:: python

   result.plot_element_field("stress_xx")

The corresponding top-level functions are:

.. code-block:: python

   from nusa import plot_element_field, plot_nodal_field

   plot_nodal_field(result, "stress_xx")
   plot_element_field(result, "stress_xx")

Nodal field plotting uses the same recovery behavior documented in
:doc:`postprocessing`.

Derived fields
--------------

Derived nodal fields can be plotted directly:

.. code-block:: python

   result.plot_nodal_field("displacement_magnitude")
   result.plot_nodal_field("von_mises_stress")

Beam diagrams
-------------

Beam results expose shear-force and bending-moment diagrams based on the frozen
end actions stored in the result:

.. code-block:: python

   result.plot_shear_diagram()
   result.plot_moment_diagram()

Equivalent top-level functions are also available:

.. code-block:: python

   from nusa import plot_moment_diagram, plot_shear_diagram

   plot_shear_diagram(result)
   plot_moment_diagram(result)

Axes and Matplotlib integration
-------------------------------

Visualization functions return Matplotlib axes, allowing callers to add titles,
labels, annotations, or combine NuSA output with other plotting code:

.. code-block:: python

   ax = result.plot_deformed_shape(scale=20.0)
   ax.set_title("Deformed shape")

Several functions also accept an existing ``ax`` argument when embedding NuSA
visualizations in a larger Matplotlib workflow.

Ownership rule
--------------

A useful rule is:

.. code-block:: text

   plot_model(model)      -> what was defined
   plot_*(result)         -> what was solved

This keeps presentation independent from mutable model solution state.

API reference
-------------

.. automodule:: nusa.visualization
   :members:
   :undoc-members:
   :show-inheritance:
