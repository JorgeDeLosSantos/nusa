Migrating from NuSA 0.4 to 0.5
================================

NuSA 0.5 introduces explicit, reusable ``Material`` and ``Section``
objects for physical properties. The change is intentionally incompatible
with the old element constructors.

Constructor changes
-------------------

.. code-block:: python

   from nusa import Bar, Truss, Beam, LinearTriangle, Material, Section

   # In 0.4 (no longer supported in 0.5):
   # Bar((n1, n2), 210e9, 1e-4)
   # Truss((n1, n2), 210e9, 1e-4)
   # Beam((n1, n2), 210e9, 8e-6)
   # LinearTriangle((n1, n2, n3), 210e9, 0.3, 0.01)

   # In 0.5:
   material = Material(E=210e9, nu=0.3)
   section = Section(A=1e-4, I=8e-6)

   bar = Bar((n1, n2), material=material, section=section)
   truss = Truss((n1, n2), material=material, section=section)
   beam = Beam((n1, n2), material=material, section=section)
   triangle = LinearTriangle(
       (n1, n2, n3),
       material=material,
       thickness=0.01,
   )

``Spring(nodes, ke)`` remains unchanged.

Validation and property access
------------------------------

``Material.E`` must be strictly positive and finite. ``Material.nu``,
when supplied, must satisfy ``-1 < nu < 0.5``. A ``LinearTriangle``
requires it. A ``Section`` requires ``A`` or ``I``, depending on the
element family, while ``thickness`` remains element-specific.

Read-only convenience properties retain their usual names:
``Bar.E``/``Bar.A``, ``Truss.E``/``Truss.A``,
``Beam.E``/``Beam.I`` and
``LinearTriangle.E``/``LinearTriangle.nu``.
The triangle accepts ``thickness`` in its constructor and exposes both
``thickness`` and the legacy read-only alias ``t``.

Materials and sections can be reused by multiple elements, and element
references to these objects cannot be reassigned.

The numerical formulations, model/analysis/results separation, and solution
result APIs introduced in 0.4 remain unchanged.

See :doc:`materials_sections` and :doc:`examples/index` for current examples.
