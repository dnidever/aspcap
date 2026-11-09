Batch Processing
================

Use the same methods for one spectrum or many: ``fit(spec)`` returns a
single result, while ``fit([spec1, spec2])`` returns ordered BatchItems.
The component pipeline methods follow the same convention. The first batch
implementation processes inputs **serially**, using independent FERRE runs.
It preserves order and isolates failures; it does not yet combine stars into
one FERRE execution or launch parallel workers.

Python lists
------------

For an existing Pipeline and a list of Spec1D objects:

.. code-block:: python

   items = pipeline.fit([spec1, spec2], star_ids=["star1", "star2"])
   for item in items:
       if item.success:
           print(item.star_id, item.value.parameters)
       else:
           print(item.star_id, item.error or item.value.errors)

Each BatchItem holds its input index, star_id, normal single-spectrum return
value, and any exception message. Recorded pipeline failures remain inside
``item.value.errors``; ``item.success`` is false for either kind of failure.
Inputs are not modified. Identities default to spectrum metadata, filenames,
or generated names; duplicate identities retain distinct indices.

Apply a single stage, without implicitly running earlier stages:

.. code-block:: python

   prepared = pipeline.prepare([spec1, spec2])
   normalized = pipeline.normalize([spec1, spec2])
   global_fits = pipeline.fit_global(
       normalized_specs, initial=[initial1, initial2],
   )
   elements = pipeline.fit_elements(
       normalized_specs, parameters=[solution1, solution2],
   )

The methods prepare, normalize, estimate_initial, fit_global, fit_elements,
and fit all accept one Spec1D or a list/tuple. A single object returns the
normal single-spectrum value; a list returns a list of BatchItems. Each value has that stage's normal return type:
normalize returns a (spectrum, continuum) pair; fit_global returns a FerreFit;
fit_elements returns an element-result dictionary. Keep only successful
items when constructing inputs for a later standalone stage. A parameter
sequence must contain one complete mapping per input, in the same order.
A single mapping can instead be shared by every spectrum. Other options,
including continuum settings and pixel weights, are shared across the batch.
Full ``fit(..., initial=[initial1, initial2])`` also supports separate starts;
omit initial to use the pipeline's configured shared starting values.

Recipe batches
--------------

.. code-block:: python

   from aspcap import PipelineRecipe

   recipe = PipelineRecipe.read("pipeline.json")
   items = recipe.fit([spec1, spec2])
   prefix = recipe.fit([spec1, spec2], stop_after=4)

``stop_after`` is a one-based inclusive step number. It runs all prerequisites
through that step, rather than running the selected step in isolation.
A successful prefix ending before any fit has ``parameters=None``. Invalid
inputs or failed stages do not prevent other spectra from being processed.
Shared recipe settings can be overridden with an ``initial`` sequence and
``star_ids`` sequence, both matching the input count.

Command line
------------

The recipe runner accepts a single apStar FITS or prepared CSV, comma-separated
paths, or a text file containing one path per line:

.. code-block:: console

   aspcap run pipeline.json "star1.fits,star2.fits" --output results.json
   aspcap run pipeline.json @spectra.txt --output results.json
   aspcap run pipeline.json --input-list spectra.txt --output results.json
   aspcap run pipeline.json @spectra.txt --stop-after 4 --output coarse-results.json

In a list file, blank lines and lines beginning with ``#`` are ignored. Relative
paths resolve beside the list file. Comma-separated paths resolve from the
working directory; use a list file for filenames containing commas. Recipe
paths remain relative to the recipe file. Input order and duplicates are kept.
FITS and CSV inputs may be mixed.

``--fits-row 0`` selects the default combined apStar row. For individual visits,
supply the actual zero-based FITS row index using ``--fits-row``. This selection
applies to all FITS files in that invocation. ``--bad-bits 255`` specifies pixel
flags to exclude; hexadecimal integers are accepted too. See
:doc:`observed-spectra` for the reader's pixel selection and mask policy.
CSV input columns are wavelength, flux, error, and optional mask. Errors are
one-sigma uncertainties, not inverse variances. CSV pixels must already match
the recipe grid; FITS pixels can be selected from the apStar wavelength scale.

A batch output JSON contains success, count, and a results list, with source
path, input index, identity, stage history, parameters, abundances, and errors
for each input. Single-file invocations retain the single-result structure.
Unreadable inputs are recorded and later inputs continue. The command exits
nonzero if any input fails, after saving all results. A malformed recipe or
unavailable grid is a shared setup error and stops before the batch starts.
Fitting does not apply survey calibrations or propagate resampling covariance.

Reference
---------

.. autoclass:: aspcap.BatchItem
   :members: success

.. autofunction:: aspcap.read_apstar

.. autofunction:: aspcap.read_spectrum
