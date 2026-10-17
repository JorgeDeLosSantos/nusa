Materials and sections
======================

NuSA 0.5 introduces reusable physical-property objects. A structural element
references a :class:`nusa.material.Material` and, where appropriate, a
:class:`nusa.section.Section`, rather than accepting individual physical
parameters in its constructor.

Material
--------

``Material(E, nu=None, name=None)`` defines an isotropic linear-elastic
material. Young's modulus ``E`` is required and must be finite and strictly
positive. Poisson's ratio ``nu`` is optional; if provided, it must satisfy
``-1 < nu < 0.5``. A triangle requires a material with ``nu`` defined.

.. code-block:: python

   from nusa import Material

   steel = Material(E=210e9, nu=0.3, name="Steel")

The material exposes read-only ``E``, ``nu``, and ``name`` properties.

Section
-------

``Section(A=None, I=None, name=None)`` stores geometric section properties:

* ``A``: cross-sectional area, required by ``Bar`` and ``Truss``.
* ``I``: second moment of area, required by ``Beam``.

At least one is required. Both must be finite and strictly positive when
provided. A section may contain both properties, but each element validates
only the property it needs.

.. code-block:: python

   from nusa import Section

   section = Section(A=2.0e-3, I=8.0e-6, name="Main member")

Properties ``A``, ``I``, and ``name`` are read-only.

Sharing physical properties
---------------------------

The same material and section can be reused across multiple elements:

.. code-block:: python

   from nusa import Bar, Material, Node, Section

   n1 = Node((0.0, 0.0))
   n2 = Node((1.0, 0.0))
   n3 = Node((2.0, 0.0))

   steel = Material(E=210e9)
   section = Section(A=1.0e-4)

   e1 = Bar((n1, n2), material=steel, section=section)
   e2 = Bar((n2, n3), material=steel, section=section)

The element's ``material`` and ``section`` references and convenience
properties such as ``E`` and ``A`` are read-only.

Element requirements
--------------------

.. list-table::
   :header-rows: 1

   * - Element
     - Construction
   * - ``Bar``
     - ``Bar(nodes, material=material, section=section)``; requires ``section.A``
   * - ``Truss``
     - ``Truss(nodes, material=material, section=section)``; requires ``section.A``
   * - ``Beam``
     - ``Beam(nodes, material=material, section=section)``; requires ``section.I``
   * - ``LinearTriangle``
     - ``LinearTriangle(nodes, material=material, thickness=t)``; requires ``material.nu``
   * - ``Spring``
     - ``Spring(nodes, ke)``; unchanged

All named physical-object arguments above are keyword-only, except the
``Spring`` stiffness argument. Boolean values (including NumPy booleans)
are not valid physical scalars. Units are chosen consistently by the user;
NuSA does not perform unit conversions.

API reference
-------------

.. autoclass:: nusa.material.Material
   :members:
   :undoc-members:

.. autoclass:: nusa.section.Section
   :members:
   :undoc-members:
