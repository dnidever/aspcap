Element Windows and Weights
===========================

An element fit should emphasize the parts of the spectrum that constrain its
abundance. Windows define which wavelength intervals are considered; pixel
weights control their relative contribution. Both work alongside the measured
uncertainties and bad-pixel masks.

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

Current abundance interface
---------------------------

``Grid.fit_abundance`` and the compatibility function ``aspcap.fit_abundance``
currently use **unit weights inside the supplied windows and zero outside**.
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
Arbitrary element-specific weight curves are not yet accepted by these
abundance interfaces.

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
