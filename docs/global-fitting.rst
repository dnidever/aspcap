Global Fitting
==============

The global stage fits the usable spectrum jointly to determine the stellar
labels that influence many spectral features. It supplies the parameter solution
used by the subsequent :doc:`element-fitting` stage.

See :doc:`initial-estimates` for starting values, multiple-start examples,
and a comparison with QASPCAP, DR17 ASPCAP, and Astra.

Labels and microturbulence
--------------------------

The planned HGE global analysis has six free labels:

* Effective temperature.
* Surface gravity, log g.
* Overall metallicity, [M/H].
* Alpha enhancement, [alpha/M].
* Carbon, [C/M].
* Nitrogen, [N/M].

Microturbulence is prescribed by a relation with log g instead of being an
independent fitted label. No particular relation or coefficients are implemented
by the generic ``Grid.fit`` interface. A six-dimensional library can encode
this prescription when its spectra are calculated. A library with a separate
microturbulence axis would require a fitting implementation that ties that axis
to the trial log g throughout optimization; fixing it at the starting value is
not equivalent. ``Grid.fit`` now accepts explicit affine ties; see :doc:`grid`. The relation
must be expressed in the actual stored coordinate, including any logarithm.

``Grid.fit`` is a general fitting tool: it varies all labels in the supplied
grid by default, or only those named in ``free_labels``. It does not automatically
identify the six HGE labels or apply physical relations between them.

Prepare the spectrum and grid
-----------------------------

Use a library covering the expected stellar parameters and having an appropriate
:doc:`lsf`. Place the observation on its native rest-frame wavelength pixels,
using the same units and air/vacuum convention. Normalize observations and models
consistently; see :doc:`continuum`. Supply the normalized one-sigma errors along
with the flux. Exclude gaps and bad pixels using zero weights.

``Grid.fit`` accepts either a ``Spec1D`` object or separate flux/error arrays.
With Spec1D, its errors and mask are used automatically and its wavelengths
are checked against the grid. Saved native mappings restore detector gaps;
otherwise orders concatenate without padding. See :doc:`spectrum`.
Neither the fit nor the mapping applies a radial-velocity correction or an LSF
convolution. The fit does not estimate a continuum.

A runnable example
------------------

From the repository checkout:

.. code-block:: python

   import numpy as np
   from aspcap import Grid

   grid = Grid("examples/toy-grid.dat")
   truth = {"ALPHA": 0.20, "METALS": -0.15}
   spec = grid(truth)
   spec.err = np.full(grid.npix, 0.01)
   fit = grid.fit(
       spec,
       start={"ALPHA": 0.0, "METALS": -0.25},
       free_labels=("ALPHA", "METALS"),
       weights=[1.0] * grid.npix,
   )
   print(fit.parameters)
   print(fit.errors)
   print(fit.reduced_chi2)

This small two-label example demonstrates the interface; it does not represent
a six-label stellar analysis or an independent interpolation-accuracy test.
The :doc:`grid` page also provides an installed three-dimensional APOGEE example.
For a real global fit, supply a scientifically suitable grid and starting values
keyed by its exact, case-sensitive ``grid.labels``. Include starting values for
both free and fixed labels. Omitted ``free_labels`` means all labels are free.

Results and diagnostics
-----------------------

The returned ``FerreFit`` contains ``parameters``, ``errors``,
``reduced_chi2``, and the native-grid ``model``. Inspect residuals across chips
and important spectral features, the location of the solution within grid
bounds, and sensitivity to starting values. Errors are conditional on any
fixed labels and the supplied continuum and do not include all systematic
uncertainties. An unavailable positive FERRE uncertainty is returned as None.

Use ``keep_runs=True`` to retain FERRE inputs and logs after success.
Failed runs retain their files and report their location. Current wrapper
calls use linear grid interpolation and disable FERRE's internal continuum
adjustment. Use :doc:`pipeline` to coordinate normalization, this global fit, and element
measurements from Python. A global-fitting command is not yet provided.
