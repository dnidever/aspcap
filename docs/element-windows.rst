Element Windows and Weights
===========================

An element fit should emphasize the parts of the spectrum that constrain its
abundance. Windows define which wavelength intervals are considered; pixel
weights control their relative contribution. Both work alongside the measured
uncertainties and bad-pixel masks.

**The recommended format uses one file per element**, containing wavelength
and weight columns. Positive weights define the included windows and their
relative importance; zero weights exclude pixels. No separate window file is
needed. Interval-only files and the original DR17 files remain supported as
alternatives.

Windows, weights, and uncertainties
-----------------------------------

* **Windows** are inclusive wavelength intervals containing useful features of
  the target element. Multiple intervals form a union: overlaps do not count
  a pixel twice. Bounds use the grid's rest-frame wavelength units and its
  air/vacuum convention.
* **Weights** are nonnegative, dimensionless values on the fitting pixels.
  Zero excludes a pixel; larger relative values emphasize it. A smooth curve
  can downweight blends or uncertain parts of a feature.
* **Uncertainties** describe measurement noise. The contribution of a pixel to
  the fit is proportional to its window weight divided by its variance. Do not
  put inverse variance into the additional weights a second time.
* **Masks** identify unusable data, such as detector gaps and flagged pixels.
  In Spec1D and the abundance interface, True means bad. A window does not
  make a masked pixel usable.

Schematically, the objective has the form:

.. math::

   \chi^2 \propto \sum_i w_i\,
   \frac{(f_i-m_i)^2}{\sigma_i^2}.

Here flux, model, and one-sigma errors must be on compatible normalized scales.
FERRE normalizes weights internally; use its reported reduced chi-square
rather than assuming this expression specifies that statistic's normalization.

Choosing useful regions
-----------------------

Prefer features sensitive to the target abundance, with limited contamination
from other elements, uncertain lines, or sky/telluric residuals. Compare model
spectra with perturbed abundances to assess sensitivity and blends. Check the
windows across the temperature, gravity, and metallicity range being analyzed;
a useful line in one star can be weak or blended in another.

When a shared alpha or metallicity dimension represents the target abundance,
other model elements change too, making window selection especially important.
A dedicated mini-grid changes one element, but its features can still be blended.
Inspect residuals and recovery on reference spectra before adopting a window set.
The bundled toy windows demonstrate the interface; they are not a validated HGE
line-selection prescription.

DR17 windows and weights
------------------------

The package bundles the original files from the DR17 directory of
`sdss/apogee <https://github.com/sdss/apogee/tree/e134409dc14b20f69e68a0d4d34b2c1b5056a901/data/windows/dr17>`_.
They are included in installed wheels, so no download is required.

The ASPCAP method constructed pixel weights from the sensitivity of synthetic
flux to the target abundance, penalizing sensitivity to competing elements and
poor agreement with reference spectra. This is described in
`Garcia Perez et al. (2016) <https://doi.org/10.3847/0004-6256/151/6/144>`_.
For DR17, the archived files also incorporate manual window selection and
bad-region exclusions. Most windows were retained from DR16, with additional
Ce regions, as described in the
`DR17 paper <https://arxiv.org/abs/2112.02026>`_.

* ``.filt`` contains the original fractional pixel weights.
* ``.wind`` contains wavelength intervals and editing flags: 0 disables,
  1 retains, and 2 adds a region with unit weights.
* ``.mask`` contains the final weights supplied to FERRE, preserving fractional
  weights after window edits and bad-region exclusions. Use these by default.

The loader returns weights on the 7514-pixel DR17 vacuum Angstrom scale, with
three chips concatenated. It can select matching grid pixels or expand onto
the 8575-pixel apStar scale with zero weights in the gaps:

.. code-block:: python

   from aspcap import windows

   wave, weights = windows.load_dr17("Mg")
   original_wave, original_weights = windows.load_dr17("Mg", kind="filt")
   edited_regions = windows.load_dr17("Mg", kind="wind")

For an existing compatible APOGEE grid, use
``_, weights = windows.load_dr17("Mg", wavelength=grid.wavelengths)`` and
``grid.fit(spec, start=parameters, free_labels="ALPHA", weights=weights)``
if that grid represents Mg with its ``ALPHA`` label. Label names and the
conversion to [Mg/H] depend on the library. Wavelengths must match native DR17
pixels; the loader does not interpolate weights onto arbitrary spectra.

The archive includes 27 species, including unused historical placeholders.
DR17 did not attempt Ge, Rb, or Yb, and judged P, Cu, Nd, and C13 unreliable
(`DR17 abundance guidance <https://www.sdss4.org/dr17/irspec/abundances/>`_).
These are historical prescriptions, not a validation for every HGE spectrum.
Source commit, checksums, and the SDSS BSD license accompany the bundled files.

Loading your own files
----------------------

Use one two-column wavelength-weight table for new prescriptions::

   # wavelength_angstrom weight
   15749.500 0.0
   15749.718 0.4
   15749.935 1.0
   15750.153 0.7

Rows contain
wavelength and a finite nonnegative weight; zero excludes the pixel. Files use
whitespace separators and may contain ``#`` comments. Wavelength units and the
air/vacuum convention must match the model grid.

.. code-block:: python

   from aspcap import windows

   wave, weights = windows.load_weights("Mg_weights.txt")
   windows.save_weights("Mg_weights_copy.txt", wave, weights)

To fit, load against the actual native grid wavelengths. Matching may select
or reorder pixels but never interpolates or fills missing coverage:

.. code-block:: python

   _, weights = windows.load_weights(
       "Mg_weights.txt", wavelength=grid.wavelengths
   )
   fit = grid.fit(spec, start=parameters, free_labels="ALPHA", weights=weights)

For legacy one-column ``.mask`` or ``.filt`` files, provide their source scale
with ``wavelength=source_wave``, or explicitly use ``scale="dr17"`` for the
7514-pixel DR17 scale. With ``scale="dr17"``, an optional ``wavelength=`` selects
matching target pixels. File suffixes do not determine the interpretation.

Simple interval files use a separate reader, since two columns of interval
bounds cannot safely be distinguished from wavelength-weight columns:

.. code-block:: python

   intervals = windows.load_windows("Mg_intervals.txt")
   # Pass windows=intervals to grid.fit_abundance(), or make a unit-weight curve:
   wave, weights = windows.load_windows(
       "Mg_intervals.txt", wavelength=grid.wavelengths
   )

Each interval row is an inclusive lower/upper pair. Overlapping intervals form
a union. A legacy three-column ``.wind`` file is an editing recipe; apply it
explicitly to its original weights:

.. code-block:: python

   wave, base = windows.load_weights("Mg.filt", scale="dr17")
   wave, weights = windows.load_legacy_windows(
       "Mg.wind", wavelength=wave, weights=base
   )
   windows.save_weights("Mg_weights.txt", wave, weights)

This applies DR17's strict-interior edits (0 disables, 1 retains, 2 adds unit
weights). Additional bad-region exclusions in a final ``.mask`` are not encoded
in ``.wind``. To convert the final DR17 prescription directly, use
``wave, weights = windows.load_dr17("Mg")`` and then ``save_weights``.

Current abundance interface
---------------------------

``Grid.fit_abundance`` and the compatibility function ``aspcap.fit_abundance``
accept **fractional pixel weights, wavelength windows, or both**.
With windows alone, included pixels receive unit weights. With both, the
windows restrict the supplied curve; they do not replace its fractional values.
Masked pixels, nonfinite flux, and nonpositive or nonfinite errors are excluded.
At least two usable window pixels are required; that count alone does not
establish abundance sensitivity or measurement reliability.

An element configuration can be written as:

.. code-block:: python

   element = {
       "species": "Mg",
       "label": "ALPHA",
       "windows": [(15002, 15005)],
       "add_labels": ["METALS"],
   }

Use it with ``grid.fit_abundance(..., **element)`` or include it in a
:doc:`pipeline` element list. Supply the spectrum, uncertainties, fixed
continuum, global parameters, and star identifier as shown in
:doc:`element-fitting`. The pipeline preserves zero global weights as element
exclusions; positive global weight magnitudes affect only the global fit.
For a single-file element prescription:

.. code-block:: python

   import numpy as np
   from aspcap import windows

   _, weights = windows.load_weights("Mg_weights.txt", wavelength=grid.wavelengths)
   result, pixels = grid.fit_abundance(
       spec.wave, spec.flux, spec.err, np.ones(grid.npix), parameters,
       star_id="star", species="Mg", label="ALPHA",
       weights=weights, add_labels=("METALS",),
   )

Use a fixed global continuum instead of ones for unnormalized inputs. Pipeline
configurations accept the same ``weights`` array and optional ``windows``.
The number of pixels used counts positive weights, not their sum.

The backend writes these values to FERRE's ``FILTERFILE`` and supplies
one-sigma errors separately. FERRE combines them as weight divided by variance;
ASPCAP does not square the element weights or fold them into the errors.
This follows the weight-file mechanism in
`DR17 abundance runs <https://github.com/sdss/apogee/blob/e134409dc14b20f69e68a0d4d34b2c1b5056a901/python/apogee/aspcap/aspcap.py#L889-L940>`_
and
`Astra abundance runs <https://github.com/andycasey/astra/blob/5e64e310a98f189091a0d5b1aee34952a9062872/src/astra/pipelines/aspcap/abundances.py#L140-L155>`_.

Custom pixel weights with Grid.fit
----------------------------------

For a custom curve, use the general fitting method with one free abundance
label. This runnable toy example uses the native grid and an explicit curve:

.. code-block:: python

   import numpy as np
   from aspcap import Grid

   grid = Grid("examples/toy-grid.dat")
   spec = grid({"ALPHA": 0.2, "METALS": -0.25})
   spec.err = np.full(grid.npix, 0.01)
   weights = np.array([0, 0, 0.5, 1.0, 1.0, 0.5, 0, 0])
   fit = grid.fit(
       spec,
       start={"ALPHA": 0.0, "METALS": -0.25},
       free_labels="ALPHA",
       weights=weights,
   )
   mg_h = fit.parameters["ALPHA"] + fit.parameters["METALS"]
   print(mg_h)  # About -0.05 for this toy library's abundance convention

Weight arrays must have one finite, nonnegative entry per reconstructed native
pixel, concatenated in segment order. Grid.fit combines them with a Spec1D
mask. This route returns a general ``FerreFit``; abundance-scale conversion,
element result fields, and element flags must be assembled explicitly.
Record the adopted windows/curves and their grid, wavelength, and element
conventions so the measurement can be reproduced.
