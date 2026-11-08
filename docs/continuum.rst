External continuum normalization
================================

See :doc:`ferre` for FERRE's internal processing sequence, control options,
and a link to its manual.

All functions in ``aspcap.continuum`` operate outside FERRE. They make no
executable calls, change no FERRE settings, and load no Doppler model data.
Fitting a continuum and applying it are separate operations::

    from aspcap import continuum

    estimate = continuum.fit(spec, method="chebyshev", degree=4)
    normalized = estimate.apply(spec)  # leaves spec unchanged

The default is an unweighted degree-four polynomial on scaled pixel coordinates,
without clipping. Each spectral order is independent, so detector gaps are not
bridged. The Chebyshev basis represents the usual polynomial space while
avoiding large unscaled pixel powers. Wavelength coordinates are also available;
on a regular log-wavelength grid they represent a different polynomial choice
from pixel coordinates.

Role in the fitting workflow
----------------------------

Normalize before :doc:`global-fitting`, using the same continuum convention as
the synthetic library. A continuum or pseudo-continuum is an estimate rather
than a guarantee of unabsorbed flux, especially in cool, line-rich spectra.
Compare methods on representative data and inspect residuals and fitted labels.
No estimator is universally best for every stellar type or signal-to-noise.

For :doc:`element-fitting`, retain a supplied fixed continuum so that the
abundance measurement does not silently refit the normalization in each window.
If the data have already been normalized, pass a continuum of ones to
``fit_abundance``. Its division by the supplied continuum is independent of
``continuum.fit``; the abundance function does not call a continuum estimator.

Keep orders or chips separate during estimation. Preserve bad-pixel masks,
propagate the flux scaling to one-sigma errors, and inspect the estimated
continuum before applying it. Applying it does not propagate uncertainty in the
continuum estimate. Normalization does not correct radial velocity, wavelength
convention, or an :doc:`lsf` mismatch.

Available methods
-----------------

* ``chebyshev`` or ``polynomial``: direct polynomial least squares. Set
  ``weighted=True`` for inverse-variance fitting, or ``clip=(lower, upper)``
  for optional iterative robust residual rejection. Each side can be None.
* ``spline``: smoothing spline, degree three by default. ``smoothing`` is the
  permitted weighted residual sum bound; zero interpolates the retained points.
  The default noise-based bound uses errors when provided, or robust adjacent
  differences. Smoothness must be chosen for the spectrum; it is not universally
  suitable for dense molecular absorption.
* ``median``: running median of observed flux, after interpolation over bad
  samples within each order. ``window`` is an odd number of pixels.
* ``percentile``: positive-flux bin percentiles, linearly interpolated with
  constant endpoint extension. Defaults are 20 bins and the 90th percentile.
  This is inspired by Doppler's percentile envelope, not a byte-for-byte port
  of its quadratic interpolation and outlier preprocessing.
* ``model_ratio``: running median of observed/model flux ratios. Inputs must
  already be comparable in normalization, wavelength, frame, and resolution.
  The divisor is a dimensionless correction, not an independently determined
  physical continuum. Defaults are a 151-pixel window, mirror boundaries, and
  accepted ratios between 0.1 and 10. Model masks are respected.

For example::

    polynomial = continuum.fit(spec, degree=4, weighted=True,
                               clip=(2.0, 4.0))
    smooth = continuum.fit(spec, method="spline", degree=3, smoothing=0.1)
    envelope = continuum.fit(spec, method="percentile", percentile=90, bins=20)
    correction = continuum.fit(rectified_data, method="model_ratio",
                               model=rectified_model, window=151)

Clipping thresholds above are robust residual-sigma thresholds for these Python
fits. They have no relationship to FERRE's relative-error REJECTCONT setting.
Neither the generic methods nor the defaults claim exact DR17/Astra reproduction.

Diagnostics and application
---------------------------

A ``ContinuumResult`` provides ``continuum``, ``fit_mask`` (True means an original
sample was used), per-order ``flags``, ``settings``, and fitted ``parameters``.
The fit mask differs from the spectrum's bad-pixel mask. A sample excluded from
polynomial fitting by residual clipping is not automatically discarded as a
measurement when applying a valid continuum; input defects remain masked.

Insufficient support, singular fits, or invalid continua raise a descriptive
error by default. ``on_failure="mask"`` instead leaves that order's continuum
NaN and records a flag. Applying invalid scales masks them; no failure is
silently replaced with an apparently valid unity continuum.

Applying a scale C divides both flux and standard-deviation error by C. Original
bitmasks are preserved. Full continuum normalization marks flux units as
dimensionless (an empty unit string); a model-ratio correction retains the
input flux unit. The raw unit is saved with the raw data. This rescaling does not propagate uncertainty in the
estimated continuum itself, so subsequent fit errors remain conditional on it.
An explicit subsequent correction composes with the previous continuum factor.
Arrays, order extents, and wavelengths are checked before reusing a result.

Model consistency
-----------------

These functions estimate external transforms; they do not decide whether a
synthetic grid is already normalized appropriately. If data and models need
independent polynomial normalization, use consistent pixel masks and keep both
factors. A model-ratio correction assumes comparable input spectra and should
not be applied twice. Holding a final continuum fixed for abundance fitting is
a separate pipeline policy from choosing any of these estimators.

See ``examples/spectrum_continuum.py`` for a runnable multi-chip example that
needs no FERRE executable, and the spectrum page for native-grid mapping.

API
---

.. automodule:: aspcap.continuum
   :members: fit, from_array, ContinuumResult
