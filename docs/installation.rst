Installation
============

Install ASPCAP and FERRE
------------------------

Use Python 3.9 or later. Install into the Python environment in which you
intend to run your analysis:

.. code-block:: sh

   python -m pip install --upgrade aspcap

Starting with ASPCAP 0.2.1, this also installs ``ferre>=1.0.4`` automatically,
alongside NumPy and SciPy. FERRE supplies the compiled fitting engine;
ASPCAP supplies the Python analysis interfaces and command-line workflow.
You do not need to install FERRE separately for a normal installation.

If you still have ASPCAP 0.2.0, upgrade with the command above. To upgrade
FERRE independently later, run:

.. code-block:: sh

   python -m pip install --upgrade ferre

An isolated environment is recommended. For example:

.. code-block:: sh

   python -m venv aspcap-env
   source aspcap-env/bin/activate
   python -m pip install --upgrade pip
   python -m pip install --upgrade aspcap

The activation command above is for macOS/Linux shells. On Windows, a virtual
environment normally uses ``aspcap-env\Scripts\activate``; see the platform
limitations below before installing FERRE there.

Precompiled FERRE and platform support
--------------------------------------

FERRE 1.0.4 provides precompiled wheels for:

* macOS 13 or later, Intel x86-64.
* macOS 13 or later, Apple Silicon arm64.
* Linux x86-64 with glibc 2.27 or later.

On these platforms, pip can install FERRE and its bundled runtime libraries
without a local Fortran compiler. Keep pip up to date so it recognizes the
available wheels. NumPy and SciPy also need compatible distributions for your
Python version and platform.

On platforms without a matching FERRE wheel, pip may attempt a source build.
That requires GNU Fortran and an appropriate build environment. Native Windows
and Linux ARM wheels are not supplied in FERRE 1.0.4; do not assume the ordinary
install command will work there without additional setup. For Windows, a
compatible Linux x86-64 environment under WSL is an alternative to investigate.
See the `FERRE repository <https://github.com/dnidever/ferre>`_ for build details.

Older macOS Python distributions can report macOS 11 or later as 10.16,
causing pip to miss a compatible wheel. If pip unexpectedly starts building
FERRE on a supported Mac, check ``python -c "import platform; print(platform.mac_ver())"``
and use a current Python distribution that reports the actual OS version.

Verify the installation
-----------------------

.. code-block:: sh

   python -m pip check
   python -c "from importlib.metadata import version; print('ASPCAP', version('aspcap')); print('FERRE', version('ferre'))"
   python -c "from ferre import get_executable; print(get_executable())"
   aspcap --help

The executable is named ``ferre``, without an ``.exe`` suffix in ASPCAP's
default configuration. Installing FERRE supplies that command in the active
environment's scripts directory. Activate the environment before running fits.

If ASPCAP cannot find it on PATH, obtain the full path using ``get_executable``
and pass it explicitly:

.. code-block:: python

   from aspcap import Grid
   from ferre import get_executable

   grid = Grid("/path/to/library.hdr", executable=get_executable())

For the command-line workflow, use ``--ferre /full/path/to/ferre``. Installing
the engine does not automatically load a production stellar library or correct
an observation's velocity, LSF, or normalization. See :doc:`grid`, :doc:`ferre`,
and :doc:`quickstart` for those next steps. The small bundled APOGEE grid is
included as package data for examples and tests.

Install the development version
-------------------------------

To get the current GitHub code:

.. code-block:: sh

   python -m pip install --upgrade 'git+https://github.com/dnidever/aspcap.git@main'

To modify the package locally:

.. code-block:: sh

   git clone https://github.com/dnidever/aspcap.git
   cd aspcap
   python -m pip install -e .

Both routes install the declared dependencies, including FERRE. FERRE's
Python distribution version and the compiled engine's version are different
version numbers; do not interpret ``ferre==1.0.4`` as an engine version of 1.0.4.

Build the documentation
-----------------------

From a repository checkout:

.. code-block:: sh

   python -m pip install -r docs/requirements.txt
   python -m sphinx -W -b html docs docs/_build/html
