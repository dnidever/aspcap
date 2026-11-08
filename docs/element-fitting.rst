Abundance Fitting
=================

After :doc:`global-fitting`, fit each target element in selected spectral
windows. Hold all grid dimensions fixed except the abundance dimension,
which is varied to fit the target element's abundance. Retain the supplied
continuum throughout this stage.

Use ``grid.fit_abundance(...)`` for an element measurement. It prepares the
windowed, normalized inputs, calls ``grid.fit`` with one free label, and reports
the abundance. The standalone ``aspcap.fit_abundance(grid, ...)`` remains
available and uses the same implementation.

Choosing the abundance dimension
--------------------------------

Conventional grids use an existing abundance dimension: alpha elements can use
``[alpha/M]``, carbon and nitrogen their respective dimensions, and other
elements ``[M/H]``. The exact label names come from ``grid.labels``. Changing
an alpha or metallicity dimension changes several elements in the model, so
windows must isolate features suitable for measuring the target element.

We are also developing **mini-grids**, with the six global stellar labels plus
a dedicated seventh dimension that changes the abundance of just one target
element. That element's fit holds the six global labels fixed and varies only
its dedicated abundance dimension. The present generic interface can select a
single supplied grid label, but it does not generate mini-grids, choose them
automatically, or supply a validated mini-grid analysis pipeline. A suitable
library and an explicit abundance-scale convention are required.

Windows, weights, and masks
---------------------------

``Grid.fit_abundance`` accepts inclusive wavelength windows. A usable pixel inside
any window receives weight one; a pixel outside receives zero. Overlapping
windows do not increase its weight. Measurement uncertainties also enter the
fit. A Boolean ``mask`` uses True for unusable pixels. Nonfinite flux and
nonpositive or nonfinite uncertainties are excluded.

Arbitrary element-specific weight curves are not currently exposed by
``Grid.fit_abundance``. For explicit nonnegative pixel weights and one free label,
use ``Grid.fit(..., free_labels=(label,), weights=weights)`` after preparing
the normalized arrays. This lower-level route does not perform the abundance
conversion or assemble the element-result dictionary for you.

A runnable example
------------------

Run from the repository checkout:

.. code-block:: python

   import csv
   from aspcap import Grid

   grid = Grid("examples/toy-grid.dat")
   with open("examples/toy-spectrum.csv", newline="") as handle:
       rows = list(csv.DictReader(handle))
   wave, flux, error, continuum = [
       [float(row[name]) for row in rows]
       for name in ("wavelength", "flux", "error", "continuum")
   ]
   result, pixels = grid.fit_abundance(
       wave, flux, error, continuum,
       {"ALPHA": 0.0, "METALS": -0.25},
       star_id="synthetic-star", species="Mg", label="ALPHA",
       windows=[(15002, 15005)], add_labels=("METALS",),
   )
   print(result["abundance"])
   print(result["conditional_error"])
   diagnostic_rows = list(pixels)

Supply one row per native grid pixel. Wavelengths must match the grid and be
in the stellar rest frame. Flux and error must have the same units as the
positive continuum; the function divides both by it. If arrays are already
normalized, supply a continuum of ones to avoid normalizing them again.
For the command-line equivalent and configuration files, see :doc:`quickstart`
and :doc:`inputs`.

Abundance scales
----------------

``raw_value`` is the fitted library coordinate. ``add_labels`` names distinct
fixed logarithmic labels whose global-fit values are added to it to report
``[X/H]``. In the example, an ALPHA value of 0.20 and fixed METALS of -0.25
produce a reported [Mg/H] of -0.05. This conversion is appropriate only when
those labels mean [alpha/M] and [M/H] in the supplied library.

If the free dimension already represents [X/H], leave ``add_labels`` empty.
Do not assume a mini-grid's abundance zero point or scale from its label name;
record and check the library convention. The routine does not infer it.

Outputs and interpretation
--------------------------

The result includes the species, fitted label, raw value, reported abundance,
conditional error, FERRE reduced chi-square, number of used pixels, and flags.
``GRID_EDGE`` marks a solution within one eighth of a grid step of a bound;
``NO_UNCERTAINTY`` marks an unavailable positive FERRE uncertainty. These are
limited checks rather than a complete scientific quality assessment.

Each diagnostic row contains wavelength, normalized flux, normalized error,
model flux, and window weight. Inspect the target features and blends before
accepting a measurement. Errors are conditional on the fixed global parameters
and continuum; uncertainty in those inputs and model systematics are not
propagated. Current fits use linear interpolation. Automatic iteration over an
element list and calibration of abundances remain future pipeline work.
