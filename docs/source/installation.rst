Installation
============

NuSA requires Python 3.10 or newer.

Standard installation
---------------------

Install the current release from PyPI:

.. code-block:: bash

   pip install nusa

The default package includes the finite-element models, elements, plotting,
reporting, and ``meshio`` for reading triangular mesh files.

The historical ``nusa[mesh]`` extra remains accepted during the 0.3.0
transition, but ``meshio`` is now installed by default and the extra is no
longer necessary.

Gmsh for mesh generation
------------------------

Generating meshes from geometry with ``nusa.mesh.Modeler`` requires the
external `Gmsh <https://gmsh.info/>`_ executable. Gmsh is not a Python
dependency of NuSA and is not installed by ``pip install nusa``.

After installing Gmsh, verify that the executable is visible to NuSA:

.. code-block:: bash

   gmsh --version

If that command works, ``Modeler.generate_mesh()`` should be able to find the
same executable on ``PATH``.

Windows
~~~~~~~

Two practical options are:

* install a precompiled Gmsh application from the official Gmsh downloads and
  add its executable directory to ``PATH``;
* install the conda-forge package:

  .. code-block:: powershell

     conda install -c conda-forge gmsh
     gmsh --version

Some conda installations expose Gmsh through ``gmsh.bat`` or ``gmsh.cmd``.
NuSA detects these Windows launchers and runs them through ``cmd.exe``.

macOS
~~~~~

The official Gmsh binaries can be used directly. Homebrew users can install the
command-line application with:

.. code-block:: bash

   brew install gmsh
   gmsh --version

Linux
~~~~~

Use the package supplied by your distribution or an official Gmsh binary. On
Debian/Ubuntu systems, for example:

.. code-block:: bash

   sudo apt-get update
   sudo apt-get install gmsh
   gmsh --version

Google Colab
~~~~~~~~~~~~

Colab runs on an Ubuntu-based environment, so install Gmsh in the notebook
session before generating meshes:

.. code-block:: bash

   !apt-get update -qq
   !apt-get install -y gmsh
   !gmsh --version

Then install NuSA normally:

.. code-block:: bash

   !pip install git+https://github.com/JorgeDeLosSantos/nusa.git@develop

The Gmsh installation is specific to the current Colab runtime and must be
repeated when a new runtime is created.

Custom executable location
~~~~~~~~~~~~~~~~~~~~~~~~~~

If Gmsh is installed but is not on ``PATH``, pass its executable explicitly:

.. code-block:: python

   nodes, elements = modeler.generate_mesh(
       gmsh_executable="/path/to/gmsh"
   )

Development version
-------------------

Install the current ``develop`` branch directly from GitHub:

.. code-block:: bash

   pip install "nusa @ git+https://github.com/JorgeDeLosSantos/nusa.git@develop"

For contributor/development work:

.. code-block:: bash

   git clone https://github.com/JorgeDeLosSantos/nusa.git
   cd nusa
   git checkout develop
   python -m pip install -e ".[test]"
   python -m pytest
