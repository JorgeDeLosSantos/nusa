Examples
========

The examples in this section use the same scripts that are exercised by NuSA's
test suite. The documentation therefore stays close to executable, regression-
tested code.

Each page focuses on one finite-element family and follows the same workflow:

.. code-block:: text

   define nodes -> create elements -> build model -> apply loads/constraints
        -> solve -> inspect StaticResult

.. toctree::
   :maxdepth: 1

   spring
   bar
   truss
   beam
   linear_triangle
