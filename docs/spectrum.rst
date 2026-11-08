Multi-order spectra
===================

See :doc:`observed-spectra` for apVisit and apStar products, rest-frame
velocity correction, and choosing the fitting wavelength grid.

``Spec1D`` is a lightweight spectrum container adapted from David Nidever's
Doppler design. It requires NumPy, and imports no Cannon, Payne, RV-fitting,
LSF-model, plotting, or model-download code. SciPy supplies numerical continuum
methods. Existing Doppler objects can be passed directly to those methods.

Creating orders
---------------

For one order, supply one-dimensional arrays. For multiple orders, use a NumPy
array shaped ``(pixels, orders)`` or a list of one-dimensional order arrays.
Lists represent orders, including unequal lengths::

    import numpy as np
    from aspcap import Spec1D

    lengths = [120, 100, 90]
    wave = [15000 + 200*i + np.arange(n) for i, n in enumerate(lengths)]
    flux = [np.full(n, 2.0 + i) for i, n in enumerate(lengths)]
    error = [np.full(n, 0.02) for n in lengths]
    spec = Spec1D(flux, wave=wave, err=error, instrument="example")

    print(spec.norder)           # 3
    print(spec.numpix)           # [120 100 90]
    print(spec.shape)            # (120, 3)
    print(spec[1].shape)         # (100,)

Padding is always masked. ``mask=True`` means bad, and ``bitmask`` preserves
integer flags without assigning them a survey-specific meaning. Nonfinite
flux, invalid supplied errors, and nonpositive/nonfinite wavelengths are
unusable. ``valid`` recomputes usability after direct array edits. If rectangular
input has padded wavelengths, ``numpix`` can be supplied explicitly; otherwise
trailing invalid wavelengths determine each order's extent.

``err`` is a one-sigma standard deviation, not inverse variance. It may be
omitted for unweighted estimation. ``wave_unit``, ``flux_unit``, ``wavevac``,
headers, instrument information and arbitrary ``metadata`` are retained.
``wavevac`` records a convention: assigning it does not convert coordinates.
Wavelength units, frame, and sampling must already match the intended model.

Copying and normalization
-------------------------

``copy()``, order selection, and pixel slicing return independent data. This
makes ``spec[0].flux[:] = ...`` safe for the parent spectrum. It differs from
Doppler's documented reference/view behavior. ``np.asarray(spec)`` exposes the
flux array using NumPy's normal array/copy conventions.

External continuum fitting is separate from this container::

    from aspcap import continuum

    result = continuum.fit(spec, method="chebyshev", degree=4)
    normalized = result.apply(spec)       # returns an independent copy

For Doppler-style use, ``spec.normalize(method="chebyshev", degree=4)`` updates
the object and returns it; ``inplace=False`` returns a copy. A custom
``continuum_func(spec, **options)`` can return an array or a ContinuumResult.
Options are forwarded explicitly. Changing the callable clears the cache.
``spec.cont`` lazily estimates and caches a continuum. Call
``clear_continuum()`` after changing raw measurements that have a cached fit.

The first application saves raw ``_flux`` and ``_err`` arrays. Explicit
applications of further continuum results multiply the cached continuum factors
and append their settings to ``continuum_history``. Asking ``normalize`` to
change methods on an already normalized spectrum raises an error; use the
original raw data for a different initial estimate.

Doppler interoperability
------------------------

Continuum fitting accepts a Doppler Spec1D object directly; it reads the public
flux/wave/err/mask/numpix attributes without importing Doppler. Applying the
result preserves that object's type through a deep copy by default. To obtain
the lightweight container explicitly, use ``Spec1D.from_spectrum(doppler_spec)``.

This is a focused core, not a complete replacement for Doppler. It does not
implement Doppler's readers, model preparation, LSF convolution, barycentric
corrections, or automatic air/vacuum conversions. An optional LSF object is
opaque metadata. Its order-dependent numerical interpretation is not handled
here. Doppler's original MIT notice is included as ``LICENSE-DOPPLER``.

Native grid mapping
-------------------

An apStar-style spectrum has one full wavelength array including detector gaps.
Split it using explicit, inclusive wavelength regions::

    # These indices describe our bundled 8575-pixel APOGEE example.
    # native_wave, native_flux and native_error are the original 1-D arrays.
    regions = [(native_wave[a], native_wave[a+n-1])
               for a, n in [(246, 3028), (3585, 2495), (6344, 1991)]]
    chips = Spec1D.from_native(native_wave, native_flux, err=native_error,
                              regions=regions)
    result = continuum.fit(chips, method="chebyshev", degree=4)
    native = result.apply(chips).to_native()

``to_native`` preserves original pixel positions without resampling. Uncovered
gaps have mask=True, NaN flux, and infinite supplied errors. ``flatten`` only
concatenates actual order samples; it does not insert gaps. The saved mapping
survives order selection and flattening, so a mapped selection can still be
restored with other positions masked.

For ``Grid.fit``, pass ``weights=native.valid.astype(float)``. For the existing
abundance function, pass ``mask=native.mask``. The latter still requires a
finite positive continuum array; for already normalized observations on the
grid's model scale, supply an array of ones. These choices do not change
FERRE's internal normalization settings or automatically normalize a model grid.

API
---

.. autoclass:: aspcap.spectrum.Spec1D
   :members:
