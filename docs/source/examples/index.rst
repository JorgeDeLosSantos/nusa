Examples
========

The examples in this section use the same scripts that are exercised by NuSA's
test suite. The documentation therefore stays close to executable, regression-
tested code.

Most family-specific examples follow the same workflow:

.. code-block:: text

   define nodes -> create elements -> build model -> apply loads/constraints
        -> solve -> inspect StaticResult

The meshed-plate example extends that sequence with preprocessing:

.. code-block:: text

   geometry -> Gmsh -> triangle mesh -> LinearTriangleModel
        -> solve -> post-process

.. toctree::
   :maxdepth: 1

   spring
   bar
   truss
   beam
   linear_triangle
   meshed_plate
