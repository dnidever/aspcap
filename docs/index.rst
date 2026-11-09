ASPCAP for HGE
==============

ASPCAP (the APOGEE Stellar Parameters and Chemical Abundances Pipeline) compares observed spectra with synthetic spectral grids using FERRE. The workflow combines continuum normalization with two fitting stages: a global stellar-parameter fit and individual-element abundance fits.

The ASPCAP workflow
-------------------

1. **Prepare and normalize the spectrum.** Keep wavelengths, flux, uncertainties, and bad-pixel masks together, preserving separate orders or detector chips where appropriate. Estimate the continuum and divide both flux and uncertainty by it. Observations and models must use compatible normalization and wavelength conventions.

2. **Fit global stellar parameters.** Fit the usable spectrum jointly across a multidimensional synthetic grid. The planned HGE global fit has six free labels: effective temperature, surface gravity (log g), overall metallicity [M/H], alpha enhancement [alpha/M], carbon [C/M], and nitrogen [N/M]. Microturbulence is prescribed by a relation with log g rather than fitted as an independent parameter. These labels describe effects across the spectrum, including molecular absorption. The relation and its coefficients must be specified consistently with the synthetic grid; this overview does not prescribe a particular calibration.

3. **Fit individual elemental abundances.** Start from the global solution and hold all grid dimensions fixed except the abundance dimension, which is varied to fit the target element's abundance. Fit selected element-sensitive wavelength windows, with pixel weights emphasizing useful features and excluding unsuitable pixels. Alpha elements can use [alpha/M], carbon and nitrogen their respective coordinates, and other elements [M/H]. The coordinate is a proxy for the target element within its windows: changing an alpha or metallicity coordinate changes multiple elements in the synthetic grid. Convert the fitted coordinate to [X/H] according to the grid's abundance convention and the fixed global labels.

   We are also developing an alternative called mini-grids. Each mini-grid adds a dedicated seventh dimension to the six global stellar labels, varying the abundance of a single target element. During that element's abundance fit, the six global labels are held at their global-fit values and only the dedicated element-abundance dimension is varied. This provides an element-specific fitting coordinate instead of using the shared alpha or metallicity coordinate as a proxy.

See the `DR17 ASPCAP description <https://www.sdss4.org/dr17/irspec/aspcap/>`_ for the survey workflow and its grid-dependent choices.

What this package currently provides
------------------------------------

This HGE package provides Spec1D objects, external continuum estimators, a Grid interface for FERRE interpolation and fits with multiple free labels, and fit_abundance plus a command-line interface for individual-element fits. A global fit can be performed with Grid.fit using a suitable supplied grid. Pipeline now coordinates normalization, global fitting, and the configured element fits for one prepared spectrum. Grids, starting parameters, and continuum settings are supplied explicitly; automatic input reading and grid selection are not yet provided.

The current fit_abundance interface accepts supplied global parameters and a fixed continuum. It uses unit weights inside the selected windows and zero outside, excludes bad pixels, and uses measurement uncertainties in the fit. Arbitrary element-specific weight curves are not yet exposed by this interface; Grid.fit accepts explicit pixel weights. External normalization is separate from FERRE's internal continuum controls.

See :doc:`ferre` for the fitting engine, its internal processing options,
and the FERRE manual.

Use :doc:`pipeline` to run normalization and both fitting stages together.

The main topic pages are :doc:`grid`, :doc:`lsf`, :doc:`continuum`,
:doc:`initial-estimates`, :doc:`global-fitting`, and :doc:`element-fitting`.

Start with :doc:`spectrum` for Spec1D objects, multiple orders, uncertainties,
masks, and APOGEE chip mapping. See :doc:`continuum` for polynomial, spline,
median, percentile, and model-ratio estimators, options, and worked examples.
Spectrum objects and external continuum estimation require no FERRE executable
or Doppler model files. See :doc:`grid` for interpolation and global fitting,
and :doc:`element-fitting` for the individual-element workflow.

.. toctree::
   :maxdepth: 2
   :caption: User guide

   quickstart
   installation
   observed-spectra
   ferre
   grid
   lsf
   continuum
   initial-estimates
   global-fitting
   element-fitting
   pipeline
   spectrum
   inputs
   python-api
