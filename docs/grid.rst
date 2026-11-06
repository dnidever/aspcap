Working with Grid
=================

Use :class:`aspcap.Grid` to inspect a FERRE library, generate synthetic spectra,
or fit stellar labels with FERRE. Run these examples from the repository
checkout after following the installation steps in :doc:`quickstart`. They
use the bundled synthetic library and vacuum wavelengths in ångströms.

Open and inspect a library
--------------------------

.. code-block:: python

   from aspcap import Grid

   grid = Grid("examples/toy-grid.dat")
   print(grid)
   print(grid.labels)      # ('ALPHA', 'METALS')
   print(grid.ranges)      # ((-0.5, 0.5), (-0.5, 0.0))
   print(grid.wavelengths) # Native wavelength pixels

For a binary library, give the path to its ``.hdr`` file and keep the matching
``.unf`` file beside it. If FERRE is outside your PATH, pass
``executable="/path/to/ferre.x"`` when opening the grid. You can inspect the
header before installing FERRE. Model evaluation and fitting use the executable.

Evaluate synthetic spectra
--------------------------

A mapping is convenient when you know the labels by name. A positional
vector follows ``grid.labels``. Pass a sequence of vectors or mappings to
evaluate a batch with one FERRE invocation.

.. code-block:: python

   model = grid({"ALPHA": 0.20, "METALS": -0.25})
   wavelengths, flux = model["wave"], model["flux"]

   batch = grid([(0.20, -0.25), (0.00, -0.25)])
   first_flux, second_flux = batch["flux"]

   resampled = grid((0.20, -0.25), wave=[15002.5, 15004.5])
   flux_only = grid.model([15002.5, 15004.5], 0.20, -0.25)

The returned flux stays on the library's scale. Output wavelengths must use
the same units and air or vacuum convention as the library and fall within
its coverage. For repeated evaluation at the same wavelengths, pass
``outwave`` when opening the grid. ``model(wave, *pars)`` provides a
wavelength-first callback for an external optimizer.

Fit a spectrum on the native grid
---------------------------------

Supply one flux and one one-sigma uncertainty per native pixel, in the order
of ``grid.wavelengths``. Prepare the spectrum on the library's rest-frame
wavelength grid and flux scale before fitting. The constructor's ``outwave``
setting only affects model evaluation.

This example makes a synthetic observation with small residuals, then fits
both labels from the default grid midpoint.

.. code-block:: python

   observed = list(grid((0.20, -0.15))["flux"])
   observed[3] -= 0.001
   observed[4] += 0.001
   error = [0.01] * grid.npix

   fit = grid.fit(observed, error)
   print(fit.parameters)  # Near ALPHA = 0.20 and METALS = -0.15
   print(fit.errors)
   print(fit.reduced_chi2)
   best_fit_flux = fit.model

To freeze a label, give its value in the full starting vector and list only
the labels FERRE should vary. A zero weight excludes a pixel. For example,
this fit varies ``ALPHA`` only and uses four pixels around the toy feature.

.. code-block:: python

   fit = grid.fit(
       observed, error,
       start={"ALPHA": 0.0, "METALS": -0.15},
       free_labels="ALPHA",
       weights=[int(15002 <= wave <= 15005) for wave in grid.wavelengths],
   )

``fit.parameters`` contains coordinates in the library's label system.
For an element measurement that applies a fixed continuum and converts a
fitted coordinate to an abundance scale, see :doc:`python-api`. Use
``keep_runs=True`` to inspect the FERRE inputs and logs after a successful
fit. A failed FERRE run reports where it kept those files.

Fit the full APOGEE DR17 wavelength grid
----------------------------------------

Use the 5 × 5 × 5 M-giant example to fit metallicity, surface gravity, and
temperature together. Its 125 spectra cover [M/H] from −1.0 to 0.0 dex,
log g from 0.0 to 2.0, and effective temperatures from 3400 to 3800 K.
The generator reconstructs all 7,514 modeled detector pixels from the
same SDSS DR17 Synspec NLTE M-giant library as the smaller fixture below.

Install ASPCAP from the project branch to get the grid and detector
weights as package data. The installed files occupy about 4.3 MB.

.. code-block:: sh

   python -m pip install 'git+https://github.com/dnidever/aspcap.git@hge/first-abundance-fit'

To reproduce the grid from SDSS, run
``python examples/download_dr17_mgiant.py dr17-grid`` from the checkout
with NumPy installed. The script downloads the PCA header and 125
coefficient records via HTTP range requests.

The output uses the 8,575 native apStar wavelength pixels, from about
15100.8 to 16999.8 Å in vacuum. Pixels outside the three modeled detector
segments have placeholder flux of one and a weight of zero. Always pass
the provided weights to the fit, combining them with any additional mask
for your observation.

.. code-block:: python

   import json
   from importlib.resources import files

   from aspcap import Grid

   directory = files("aspcap").joinpath("data")
   grid = Grid(directory / "apogee-dr17-mgiant-full.hdr")
   weights = json.loads(
       (directory / "apogee-dr17-mgiant-full-weights.json").read_text()
   )

   truth = {"METALS": -0.32, "LOGG": 1.14, "TEFF": 3625}
   observed = list(grid(truth)["flux"])
   observed[600] += 0.0001
   observed[4200] -= 0.0001
   fit = grid.fit(observed, [0.003] * grid.npix, weights=weights)
   print(fit.parameters)  # Near -0.32, 1.14, and 3625 K

The binary data use little-endian, 32-bit floats as expected by the
standard single-precision FERRE build.

Fit a small APOGEE DR17 M-giant grid
------------------------------------

The checkout includes ``examples/apogee-dr17-mgiant-mini.dat``, a compact
ASCII FERRE library derived from the SDSS DR17 APOGEE Synspec NLTE M-giant
grid. It contains 27 spectra on a 3 × 3 × 3 grid.

* ``METALS`` ([M/H]) runs from −0.50 to 0.00 dex in steps of 0.25 dex.
* ``LOGG`` runs from 0.5 to 1.5 in steps of 0.5 dex.
* ``TEFF`` runs from 3500 to 3700 K in steps of 100 K.

Each spectrum retains 160 native vacuum pixels from about 15736.6 to
15771.2 Å. The remaining library axes are held fixed at
``LOG10VDOP = 0.38021``, ``C = 0``, ``N = 0``, and ``O Mg Si S Ca Ti = 0``.
Use this narrow subset for software tests and demonstrations.

The following example interpolates a synthetic observation between grid
nodes, adds small residuals, then fits all three labels from their grid
midpoints.

.. code-block:: python

   from aspcap import Grid

   grid = Grid("examples/apogee-dr17-mgiant-mini.dat")
   print(grid)
   truth = {"METALS": -0.32, "LOGG": 1.14, "TEFF": 3625}
   observed = list(grid(truth)["flux"])
   observed[60] += 0.0001
   observed[90] -= 0.0001

   fit = grid.fit(observed, [0.003] * grid.npix)
   print(fit.parameters)  # About -0.32, 1.14, and 3630 K

The fixture reconstructs blue-segment pixels 2739 through 2898
(zero-based) from the archived PCA coefficients and basis vectors. It
stores the resulting flux in ASCII, with the original wavelength spacing
and resolution. The source is the APO LSF-a library
``p_apssgM_200921nlte_lsfa_012_075``.
See the SDSS `DR17 APOGEE library description
<https://www.sdss4.org/dr17/irspec/apogee-libraries/>`_ and the
`source archive directory
<https://data.sdss.org/sas/dr17/apogee/spectro/speclib/synth/synspec/marcs/giantisotopes/sgM_200921nlte_lsfa/>`_.

API reference
-------------

.. autoclass:: aspcap.Grid
   :members: read, model, fit
   :special-members: __init__, __call__, __repr__

.. autoclass:: aspcap.ferre.FerreFit
