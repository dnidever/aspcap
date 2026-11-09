Run the Pipeline
================

``Pipeline`` coordinates spectrum validation, external normalization, global
fitting, and the configured element fits. It reuses the existing ``Spec1D``,
continuum, and ``Grid`` interfaces, preserving the input and intermediate results.
Each single-spectrum call processes one prepared spectrum with explicitly supplied
grids, starting parameters, and element prescriptions.

A complete example
------------------

From the repository checkout, this example creates a synthetic observation
and fits it through both stages:

.. code-block:: python

   import numpy as np
   from aspcap import Grid, Pipeline

   grid = Grid("examples/toy-grid.dat")
   spectrum = grid({"ALPHA": 0.2, "METALS": -0.25})
   spectrum.flux *= 2.0
   spectrum.err = np.full(grid.npix, 0.02)

   pipeline = Pipeline(
       grid=grid,
       elements=[{
           "species": "Mg",
           "label": "ALPHA",
           "windows": [(15002, 15005)],
           "add_labels": ["METALS"],
       }],
       continuum=np.full(grid.npix, 2.0),
   )
   result = pipeline.fit(
       spectrum,
       initial={"ALPHA": 0.0, "METALS": -0.25},
       star_id="synthetic-star",
   )
   print(result.success)
   print(result.parameters)  # ALPHA about 0.2, METALS about -0.25
   print(result.abundances)  # Mg about -0.05 dex

The same example is available as ``examples/pipeline.py``. This toy grid has
two labels; it demonstrates orchestration rather than a complete HGE giant-star
library or a test of interpolation accuracy.

Inputs and configuration
------------------------

* **Spectrum:** a ``Spec1D`` with errors, native rest-frame wavelengths in
  Angstroms, a matching air/vacuum convention, and compatible model flux scale
  and LSF. See :doc:`observed-spectra` and :doc:`lsf`. Saved native mappings
  restore detector gaps; otherwise orders concatenate without padding.
* **Continuum:** a mapping of options for :doc:`continuum`, a supplied array,
  a list of per-order arrays, or a ``ContinuumResult``. Array matrices use
  the spectrum's ``(pixels, orders)`` layout. Explicit ``None`` means skip
  normalization because the spectrum is already on the library's scale.
  No estimator is selected automatically.
* **Initial parameters:** a complete mapping or vector in ``grid.labels`` order,
  supplied to the constructor or to ``fit``. This version requires estimates;
  it does not silently fall back to the midpoint. See :doc:`initial-estimates`.
* **Elements:** a list of dictionaries or ``ElementConfig`` objects. A mapping
  such as ``{"Mg": {"label": "ALPHA", "windows": ..., "add_labels": ...}}``
  is also accepted. An empty element list runs only the global stage.
* **Free labels and weights:** ``free_labels`` controls the global fit;
  ``weights`` is a native-pixel vector. The default varies all supplied grid
  labels. Microturbulence relations must be encoded consistently in the grid;
  explicit ``ties`` can impose a validated relation; none is chosen automatically.

For an estimated continuum, replace the supplied array with, for example,
``continuum={"method": "chebyshev", "degree": 4, "clip": (3, 3)}``.
Choose and validate the estimator for your real spectra; the polynomial
example is not a universal normalization prescription.
Per-spectrum ``initial``, ``continuum``, and ``weights`` override constructor
settings. ``continuum=None`` explicitly disables a configured estimator.

Normalization and masks
-----------------------

The continuum stage divides both flux and errors, independently by order.
Global fitting uses the resulting spectrum. Abundance fitting then receives
that already normalized spectrum with a continuum of ones, so the correction
is not applied twice and remains fixed throughout the element stage.
The original and applied continuum are retained for inspection.

Bad-pixel masks and invalid flux/errors are excluded by the shared Grid input
validation. Zero global weights also exclude pixels from element fits.
Positive global weight magnitudes apply only to the global fit; element fits
use their own fractional weights, restricted by any configured windows. Global weights
do not replace the spectrum mask used during continuum estimation. Set the
Spec1D mask for pixels that should also be excluded from that stage.

Inspecting results and failures
-------------------------------

``PipelineResult`` retains:

* ``original_spectrum`` and ``prepared_spectrum``: independent input snapshots.
* ``spectrum`` and ``continuum``: normalized spectrum and its ``ContinuumResult``;
  the latter is None when normalization was skipped.
* ``initial_parameters`` and ``global_fit``: supplied starts and the ``FerreFit``.
  ``parameters`` is a convenience copy of the fitted global labels, and
  ``residuals`` gives observed-minus-model values at used native pixels.
* ``element_fits``: per-element results, reusable pixel rows, flags, and errors.
  ``abundances`` maps attempted species to [X/H], or None on an element failure.
* ``flags``, ``errors``, and ``config``: processing diagnostics and the settings
  used in the run. Inspect scientific residuals and component flags as well.

``success`` is True when the global fit and all requested element fits complete.
An individual element failure preserves the global fit and other abundances,
with ``PARTIAL_ABUNDANCES`` and an ``element:<species>`` error entry.
A failure in preparation, normalization, initial estimates, or global fitting
stops subsequent stages while retaining products already computed. The failed
stage appears in ``errors`` and ``flags``. Set ``raise_on_error=True`` to re-raise
expected data/FERRE errors for debugging. Invalid input types and malformed
constructor configuration raise immediately.

Errors on measured labels and abundances remain conditional on the fixed inputs
and continuum; the pipeline does not add systematic or continuum uncertainties.
Use ``keep_runs=True`` and optionally ``run_root`` on the constructor to retain
successful FERRE inputs and logs. Failed FERRE runs retain their own files.

Rerun stages independently
--------------------------

Each stage is callable. Continuing the complete example:

.. code-block:: python

   prepared = pipeline.prepare(spectrum)
   normalized, estimate = pipeline.normalize(prepared)
   initial = pipeline.estimate_initial(normalized, initial=result.initial_parameters)
   global_fit = pipeline.fit_global(normalized, initial)
   element_fits = pipeline.fit_elements(
       normalized, global_fit.parameters, star_id="synthetic-star"
   )

To repeat just Mg, call ``fit_elements(result.spectrum, result.parameters,
species="Mg", star_id=result.config["star_id"])``. Keep the normalized spectrum and its global solution
together; using different continuum choices changes the measurement.

Supplied mini-grids
-------------------

An ``ElementConfig`` can specify a separate ``grid`` and an ``initial`` value
for its free abundance label. Fixed labels are copied by exact name from the
global solution. If the dedicated abundance label is absent from that solution,
its explicit initial value is required. The supplied grid must use matching
native pixels, wavelength convention, and compatible LSF/normalization. Set
``add_labels`` according to its abundance scale; see :doc:`element-fitting`.
No mini-grid generation, selection, or resampling is automatic.

Current scope
-------------

This interface accepts Spec1D inputs through Python. Use read_apstar for FITS
input and the batch methods for lists of spectra. :doc:`pipeline-recipes`
provides ordered stages, explicit multiple starts, and a command-line runner.
Automatic row selection, initial-estimate prediction, LSF convolution
configuration, and arbitrary observed-spectrum resampling remain future work.

API reference
-------------

.. autoclass:: aspcap.Pipeline
   :members: prepare, normalize, estimate_initial, fit_global, fit_elements, fit, fit_many, run_many

.. autoclass:: aspcap.PipelineResult
   :members: parameters, abundances, residuals, success

.. autoclass:: aspcap.ElementConfig

.. autoclass:: aspcap.ElementFitResult
   :members: success, flags

Element prescription files
--------------------------

The recommended storage format is one wavelength-weight file per element;
see :doc:`element-windows`. Load the file against the element grid, then use
an element configuration such as:

.. code-block:: python

   from aspcap import windows

   _, mg_weights = windows.load_weights("Mg_weights.txt", wavelength=grid.wavelengths)
   element = dict(species="Mg", label="ALPHA", weights=mg_weights,
                  add_labels=("METALS",))
   pipeline = Pipeline(grid, [element], initial=parameters)

Optional ``windows`` restrict the curve. Interval-only prescriptions retain
unit weights and can be loaded with ``windows.load_windows``. These are
alternatives; users do not need to maintain both files.

Parameter ties
--------------

Set ``ties=`` and ``tie_mode=`` on ``Pipeline`` for the global fit. Set them
separately on an ``ElementConfig`` or element dictionary for abundance fits.
Both use the :doc:`grid` tie syntax; the default mode is relative. Ties are
recorded in the result configuration. A metallicity-compensation example is
shown in :doc:`element-fitting`. Untied labels retain the usual fixed/free
behavior.

Text configuration
------------------

Use :doc:`pipeline-recipes` to describe ordered stages in a JSON file,
including repeated fits, multiple starts, and element weights/ties.
The basic recipe runs on bundled test data; a DR17-style template provides
a starting point for production grids.

Many spectra
------------

``Pipeline.fit_many`` runs complete pipelines; ``Pipeline.run_many`` applies a
single named stage to a list of Spec1D objects. ``PipelineRecipe.fit_many``
executes recipes for a list, optionally stopping after a selected step.
See :doc:`batch` for Python and command-line examples and failure handling.
