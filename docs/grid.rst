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

API reference
-------------

.. autoclass:: aspcap.Grid
   :members: read, model, fit
   :special-members: __init__, __call__, __repr__

.. autoclass:: aspcap.ferre.FerreFit
