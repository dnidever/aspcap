ASPCAP for HGE
==============

ASPCAP (the APOGEE Stellar Parameters and Chemical Abundances Pipeline) compares observed spectra with synthetic spectral grids using FERRE. The workflow combines continuum normalization with two fitting stages: a global stellar-parameter fit and individual-element abundance fits.

The ASPCAP workflow
-------------------

1. **Prepare and normalize the spectrum.** Keep wavelengths, flux, uncertainties, and bad-pixel masks together, preserving separate orders or detector chips where appropriate. Estimate the continuum and divide both flux and uncertainty by it. Observations and models must use compatible normalization and wavelength conventions.

2. **Fit global stellar parameters.** Fit the usable spectrum jointly across a multidimensional synthetic grid. A typical seven-dimensional giant-star grid varies effective temperature, surface gravity, microturbulence, overall metallicity [M/H], alpha enhancement [alpha/M], carbon [C/M], and nitrogen [N/M]. These labels describe effects across the spectrum, including molecular absorption. The free labels depend on the grid and stellar population; not every fit needs seven dimensions.

3. **Fit individual elemental abundances.** Start from the global solution and fix every grid label except the coordinate used to represent the target abundance. Fit selected element-sensitive wavelength windows, with pixel weights emphasizing useful features and excluding unsuitable pixels. Alpha elements can use [alpha/M], carbon and nitrogen their respective coordinates, and other elements [M/H]. The coordinate is a proxy for the target element within its windows: changing an alpha or metallicity coordinate changes multiple elements in the synthetic grid. Convert the fitted coordinate to [X/H] according to the grid's abundance convention and the fixed global labels.

See the `SDSS ASPCAP description <https://www.sdss4.org/dr17/irspec/aspcap/>`_ for the survey workflow and its grid-dependent choices.

What this package currently provides
------------------------------------

This HGE package provides Spec1D objects, external continuum estimators, a Grid interface for FERRE interpolation and fits with multiple free labels, and fit_abundance plus a command-line interface for individual-element fits. A global fit can be performed with Grid.fit using a suitable supplied grid. An automatic end-to-end pipeline that selects grids, normalizes spectra, runs the global fit, and loops over elements is not yet provided.

The current fit_abundance interface accepts supplied global parameters and a fixed continuum. It uses unit weights inside the selected windows and zero outside, excludes bad pixels, and uses measurement uncertainties in the fit. Arbitrary element-specific weight curves are not yet exposed by this interface; Grid.fit accepts explicit pixel weights. External normalization is separate from FERRE's internal continuum controls.

Start with :doc:`spectrum` for Spec1D objects, multiple orders, uncertainties,
masks, and APOGEE chip mapping. See :doc:`continuum` for polynomial, spline,
median, percentile, and model-ratio estimators, options, and worked examples.
Spectrum objects and external continuum estimation require no FERRE executable
or Doppler model files. See :doc:`grid` for interpolation and global fitting,
and :doc:`quickstart` for the individual-element workflow.

.. toctree::
   :maxdepth: 2
   :caption: User guide

   quickstart
   inputs
   grid
   spectrum
   continuum
   python-api
