Installation
============

NuSA requires Python 3.10 or newer.

Core installation
-----------------

Install the current release from PyPI:

.. code-block:: bash

   pip install nusa

The core package includes the finite-element models, elements, plotting, and reporting
functionality.

Mesh utilities
--------------

Mesh helpers are provided as an optional dependency:

.. code-block:: bash

   pip install "nusa[mesh]"

The mesh extra installs ``meshio``. Generating meshes also requires the external
`Gmsh <https://gmsh.info/>`_ executable to be installed and available on ``PATH``.

Development version
-------------------

Install the current ``develop`` branch directly from GitHub:

.. code-block:: bash

   pip install "nusa[mesh] @ git+https://github.com/JorgeDeLosSantos/nusa.git@develop"

For contributor/development work:

.. code-block:: bash

   git clone https://github.com/JorgeDeLosSantos/nusa.git
   cd nusa
   git checkout develop
   python -m pip install -e ".[test]"
   python -m pytest
