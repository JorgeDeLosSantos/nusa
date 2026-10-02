Spring elements
===============

A spring element represents a one-dimensional elastic connection with stiffness
``k``. Each spring node contributes one displacement degree of freedom,
``ux``, with corresponding nodal force ``fx``.

Use :class:`nusa.model.SpringModel` together with
:class:`nusa.element.Spring`.

Element definition
------------------

.. code-block:: python

   e1 = Spring((n1, n2), k=1000.0)

For a two-node spring, the local solved response contains the two end forces:

.. code-block:: python

   result.element_result(e1)

which provides ``force_i`` and ``force_j``.

Example: spring assemblage
--------------------------

This example contains three springs and four nodes. Two nodes are restrained and
a horizontal force is applied to the internal assemblage.

.. literalinclude:: ../../../examples/spring/spring_01.py
   :language: python
   :linenos:

Important result queries
------------------------

Nodal displacement:

.. code-block:: python

   result.displacement(n3)["ux"]

Support reaction:

.. code-block:: python

   result.reaction(n1)["fx"]

Element end forces:

.. code-block:: python

   result.element_result(e1)

A spring model is a useful first NuSA problem because it exposes the complete
``Model -> solve -> StaticResult`` workflow with only one DOF per node.
