FERRE: the fitting engine
=========================

FERRE is the Fortran model-interpolation and optimization engine used by
ASPCAP. Given observations, uncertainties, weights, and a numerical library of
synthetic spectra, it searches for labels whose predicted spectrum minimizes
the fitting statistic. ASPCAP supplies the surrounding spectral objects,
preparation, configuration, and interpretation of stellar and elemental results.

The `FERRE manual (PDF)
<https://github.com/dnidever/ferre/blob/master/docs/ferre.pdf>`_ describes its
control files, formats, and options. See also the
`FERRE repository <https://github.com/dnidever/ferre>`_. The summary below
covers capabilities of the executable; not all are exposed by ASPCAP's current
Python interface. Consult the manual and source corresponding to your installed
FERRE version when configuring direct runs.

Rest-frame wavelengths are required
-----------------------------------

For the ASPCAP workflow, remove the stellar radial-velocity Doppler shift
before fitting: the observed spectrum must be on a **stellar rest-frame
wavelength scale**, matching the model's wavelength units and air/vacuum
convention. A barycentric correction alone does not remove the star's radial
velocity. Record the correction applied to the wavelength array.

Rest frame does not mean identical wavelength sampling. FERRE can resample a
model onto a different observed wavelength array, but resampling itself does
not estimate or remove a Doppler shift. The current ASPCAP wrapper does not
fit radial velocity or shift wavelengths for you.

What happens to a trial model
-----------------------------

In the inspected FERRE implementation, model evaluation proceeds through:

1. Interpolation in the stellar-label grid at the trial parameter values.
2. PCA reconstruction, when using a supported compressed library and evaluating
   in reconstructed flux space.
3. Optional internal LSF convolution of the model spectrum.
4. Optional wavelength resampling onto the observed pixels.
5. Optional boxcar filtering and continuum treatment.
6. Comparison with the processed observation, using errors and pixel weights.

The optimizer repeats this evaluation as it varies the free labels. A model
can also be evaluated without fitting labels. Fixed labels remain at their
supplied values. The exact processing depends on the control settings and
library format; PCA-space fitting is a separate configuration.

Interpolation between grid nodes
--------------------------------

``INTER`` selects interpolation in stellar-label space:

.. list-table:: Grid interpolation
   :header-rows: 1
   :widths: 20 80

   * - Setting
     - Method
   * - ``INTER <= 0``
     - Nearest neighbor.
   * - ``INTER=1``
     - Linear interpolation.
   * - ``INTER=2``
     - Quadratic Bézier interpolation.
   * - ``INTER=3``
     - Cubic Bézier interpolation.
   * - ``INTER=4``
     - Cubic splines.

Higher-order methods use larger neighboring stencils and need enough valid
nodes in each relevant dimension. They can improve accuracy for smooth spectra,
but should be tested against independently calculated spectra. Cubic is not
a guarantee of the smallest error for every grid. ``MONO`` provides a
monotonicity option for cubic Bézier interpolation; see the manual for its use.
These settings differ from wavelength interpolation below. See :doc:`grid`
for grid preparation and the current wrapper's format restrictions.

LSF convolution inside FERRE
----------------------------

FERRE can convolve the model internally during each trial evaluation. Supply
a sufficiently high-resolution library and LSF information in ``LSFFILE``;
set ``LSF`` to choose the representation. This avoids requiring a different
pre-convolved library for every observation.

Modes 1–4 use one profile specification for all objects. Modes 11–14 allow
an individual LSF for each spectrum. Both wavelength-independent and
wavelength-dependent sampled profiles are supported, as are Gaussian profiles
specified by FWHM. The LSF sampling must match the synthetic library, and
Gaussian widths are in library-pixel units. An observation's LSF object must be
converted to FERRE's required input representation; it is not a generic Python
object passed to the executable. See :doc:`lsf` for mode details and the
alternative of pre-convolved APOGEE LSF groups.

Wavelength resampling
---------------------

``WINTER`` controls whether the observation and model share their sampling:

* ``WINTER=0``: no internal wavelength resampling; pixels must already match.
* ``WINTER=1``: resample observations onto the model wavelengths.
* ``WINTER=2``: resample models onto the observed wavelengths.

Supply the observed wavelengths through ``WFILE`` when required. Resampling
the noiseless model is generally preferable to interpolating noisy observations
and their uncertainties. The observed wavelengths must lie within the usable
model coverage. Resampling does not compensate for an LSF mismatch.

In the inspected executable, ``TWINTER=0`` selects linear wavelength
interpolation and ``TWINTER=1`` selects cubic-spline wavelength interpolation.
This is independent of ``INTER``. The manual's resampling discussion describes
linear interpolation; check your executable's supported controls before using
additional modes.

Internal continuum normalization
--------------------------------

This section describes **normalization inside FERRE**, separate from the
external estimators in :doc:`continuum`.

For positive ``CONT`` with ``OBSCONT`` enabled, FERRE estimates a continuum for
the observation and divides both its flux and one-sigma errors by that
continuum. Each trial model is divided by its own continuum estimated with the
same algorithm. The two continuum curves are computed separately; the intent
is to apply a consistent transformation to data and models, not to divide them
by one shared curve. Processing respects the library's spectral segments.

.. list-table:: Internal continuum choices
   :header-rows: 1
   :widths: 20 80

   * - Setting
     - Behavior
   * - ``CONT=0``
     - Disable internal continuum normalization.
   * - ``CONT=1``
     - Fit a polynomial in pixel coordinates; ``NCONT`` is its degree
       (zero gives a constant).
   * - ``CONT=2``
     - Split each spectral segment into ``NCONT`` pieces and divide by
       the mean in each piece.
   * - ``CONT=3``
     - Divide by a running mean. In the inspected implementation the
       filter width is ``NCONT + 1`` pixels.
   * - ``CONT=-1, -2, -3``
     - Apply the corresponding estimator to the model/observation ratio,
       then divide the model by the estimated ratio correction.

``OBSCONT=0`` disables the observation-normalization step while leaving the
selected model treatment active. Negative ``CONT`` modes correct the model
relative to the observation rather than independently flattening both. Their
interpretation must be considered together with ``OBSCONT``; do not assume
that every sign/setting combination has the same normalization policy.

``REJECTCONT`` is a relative observational-error threshold for selecting pixels
in the polynomial continuum fit. It is not iterative sigma clipping of
absorption lines. The inspected polynomial implementation uses a single,
unweighted polynomial fit to the retained samples. The manual describes
``NCONT`` as the running-mean width, whereas the implementation uses
``NCONT + 1``; use the source for the exact behavior of your build. The implemented
``CONT=-2`` mode uses piecewise means of the ratio, not a spline.

Use ``SFFILE`` for the processed observed spectra and ``OFFILE`` for fitted
model spectra when checking a direct FERRE run. Compare arrays processed on
compatible scales, rather than a normalized model against the untouched input.
Continuum choices can affect abundances and should be validated on representative
spectra, especially in cool stars with dense molecular absorption.

Other useful controls
---------------------

* ``INIT``, ``INDINI``, and ``NRUNS`` control initial guesses and multiple
  searches. See :doc:`initial-estimates`.
* ``NOV`` and ``INDV`` select the labels that can vary. The other labels are
  fixed at their supplied values.
* ``ALGOR`` selects the search algorithm, including supported local, global,
  and sampling approaches. ``ERRBAR`` selects uncertainty estimation.
  Settings depend on the analysis and executable version.
* ``FILTERFILE`` supplies pixel weights; uncertainties also enter the fitting
  statistic. Masking and appropriate errors remain necessary after normalization.
* ``F_FORMAT``, ``F_ACCESS``, and PCA controls govern library representation
  and memory/disk access. ``NTHREADS`` controls OpenMP parallelism over objects.

These are pointers to manual options, not a recommended single control file
for every analysis. Weights, processing choices, and free-label selections
must be consistent when comparing fit statistics.

What ASPCAP currently exposes
-----------------------------

The current ``Grid`` interface returns model Spec1D objects by default
(``as_dict=True`` requests a dictionary) and fits selected labels from Spec1D
observations or explicit flux/error arrays, with errors and weights. Its FERRE calls request ``INTER=1``, ``WINTER=0``,
``CONT=0``, ``OBSCONT=0``, ``INIT=0``, and ``NRUNS=1``. They do not currently
write the LSF inputs or expose the internal continuum and resampling controls.

Thus prepare observations on the native rest-frame grid and compatible flux
scale before using ``Grid.fit``. Optional output resampling in ``Grid.__call__``
or ``Grid.model`` is performed by the Python wrapper after model evaluation;
it is distinct from FERRE's internal ``WINTER=2`` route. ``fit_abundance``
divides the observation and errors by a supplied fixed continuum before calling
FERRE, with internal normalization disabled. A stored ``Spec1D.lsf`` is metadata
until a future interface explicitly converts it into FERRE inputs.

Implementation references
-------------------------

The processing description was checked against FERRE's model evaluation,
observation handling, continuum estimator, wavelength resampler, and control
variables. The links below pin the inspected repository revision.

* `share.f90 <https://github.com/dnidever/ferre/blob/4016d9e558d0bff6c1ec41f65bb1ebd04533efec/src/share.f90>`_
* `flx.f90 <https://github.com/dnidever/ferre/blob/4016d9e558d0bff6c1ec41f65bb1ebd04533efec/src/flx.f90>`_
* `ferre.f90 <https://github.com/dnidever/ferre/blob/4016d9e558d0bff6c1ec41f65bb1ebd04533efec/src/ferre.f90>`_
* `continuum.f90 <https://github.com/dnidever/ferre/blob/4016d9e558d0bff6c1ec41f65bb1ebd04533efec/src/continuum.f90>`_
* `wresample.f90 <https://github.com/dnidever/ferre/blob/4016d9e558d0bff6c1ec41f65bb1ebd04533efec/src/wresample.f90>`_
