Quick Start
===========

Install Python 3.9 or later and a working
`FERRE executable <https://github.com/callendeprieto/ferre>`_.
The worked example was verified with FERRE 5.1.3.

From the repository checkout, install ASPCAP and run the synthetic Mg example.

.. code-block:: sh

   python -m pip install -e .
   aspcap abundances examples/one-star.json --output abundance.csv --pixels pixels.csv

ASPCAP looks for ``ferre`` on your PATH. If it is elsewhere, append
``--ferre /path/to/ferre`` to the command.

Inspect the result
------------------

The example should recover a raw ``ALPHA`` coordinate of about 0.20 using four
pixels. The fixed global ``METALS`` coordinate is −0.25, so the reported
``[Mg/H]`` is about −0.05. The example uses synthetic data.

``abundance.csv`` contains the fitted coordinate, the reported abundance
scale and value, the conditional uncertainty, FERRE's reduced χ², the number
of pixels used, and any flags. The uncertainty is conditional on the global
parameters and continuum. ``GRID_EDGE`` indicates a solution within one-eighth
of a grid step of an edge. ``NO_UNCERTAINTY`` indicates that FERRE did not
return a usable uncertainty.

``pixels.csv`` contains the normalized observation, normalized uncertainty,
model, and weight at each wavelength. A weight of zero excludes that pixel
from the element fit.

Troubleshoot a run
------------------

* If ASPCAP cannot find the executable, pass its path with ``--ferre``.
* If the wavelengths fail validation, prepare the spectrum on the library's
  native rest-frame grid and check its air or vacuum convention.
* If the element window has fewer than two usable pixels, check the
  wavelength coverage and flux uncertainties.
* If FERRE fails, open ``ferre.log`` in the run directory reported in the
  error. To inspect the FERRE inputs after a successful run, add
  ``--keep-runs``. Use ``--run-root`` to choose an existing parent directory.
