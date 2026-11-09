Text Pipeline Recipes
=====================

A recipe is an ASCII JSON file listing the pipeline's steps in execution order.
It can be reviewed, edited, and shared without writing Python code. Relative
paths resolve beside the recipe file, independent of the working directory.
Only the supported stages run; arbitrary Python and shell commands are not
executed. The file has a version number so future changes can be checked.

Basic test
----------

From a repository checkout, run:

.. code-block:: console

   aspcap run examples/pipeline-basic.json examples/pipeline-spectrum.csv --output result.json

The spectrum is a prepared, normalized toy observation. The recipe explicitly
skips normalization, fits the two toy labels, and measures Mg. It should recover
ALPHA near 0.2, METALS near -0.25, and [Mg/H] near -0.05. The output JSON records
parameters, abundances, errors, the recipe, and a stage/candidate history.
A failed run returns a nonzero exit status and saves partial results.

The examples also ship with installed packages:

.. code-block:: python

   from importlib.resources import files
   from aspcap import PipelineRecipe

   recipe = PipelineRecipe.example("basic")
   example_dir = files("aspcap").joinpath("data/recipes")
   # pipeline-spectrum.csv and toy-grid.dat are beside the installed recipe.

For a prepared Spec1D object:

.. code-block:: python

   recipe = PipelineRecipe.read("my-pipeline.json")
   result = recipe.fit(spec, star_id="my-star")
   print(result.parameters, result.abundances, result.errors)

``initial=parameters`` on ``fit`` overrides the recipe's initial values for a
particular spectrum. The input object is copied. The command-line reader accepts
a flat CSV with ``wavelength,flux,error`` and optional ``mask`` (1/true means bad).
Use ``--wavelength-frame air`` for air wavelengths; vacuum is the default.
The runner accepts prepared rest-frame spectra matching native grid pixels;
it does not read apStar FITS or perform velocity correction.

File structure
--------------

This is a complete small recipe, assuming the toy grid sits beside it:

.. code-block:: json

   {
     "version": 1,
     "grids": {"toy": "toy-grid.dat"},
     "steps": [
       {"type": "prepare", "grid": "toy"},
       {"type": "normalize", "method": "skip"},
       {"type": "initial", "parameters": {"ALPHA": 0.0, "METALS": -0.25}},
       {"type": "global_fit"},
       {"type": "abundances", "elements": [
         {"species": "Mg", "label": "ALPHA",
          "windows": [[15002, 15005]], "add_labels": ["METALS"]}
       ]}
     ]
   }

JSON uses double quotes, has no trailing commas, and does not support comments.
Optional ``description`` fields document the whole recipe or individual steps;
``name`` labels a step in the output history. An optional top-level ``ferre``
sets the executable; ``--ferre`` or the Python executable argument overrides it.

Available steps
---------------

* **prepare:** select a named ``grid`` and validate the input. Required first.
* **normalize:** choose ``skip``, ``constant`` (with ``value``), or an external
  continuum method: chebyshev, polynomial, spline, median, percentile, or
  model_ratio. Put estimator settings in ``options``. model_ratio requires a
  preceding global fit and a flat native-pixel spectrum. Repeated normalization
  applies a correction to the latest spectrum, dividing flux and errors.
* **initial:** supply a complete ``parameters`` mapping in exact grid labels.
  There is no implicit stellar-parameter estimator.
* **global_fit:** fit from the current parameters, or try an explicit ``starts``
  list of complete parameter mappings. Optional ``grids`` names candidate
  libraries; their labels, pixels, and wavelength frames must match. Select
  the successful candidate with lowest reduced chi-square. This is an
  unpenalized comparison, not historical survey grid-selection logic.
  Optional settings are ``free_labels``, ``weights``, ``ties``, and ``tie_mode``.
* **abundances:** fit an ``elements`` list using the latest global solution and
  fixed current continuum. Elements accept label, windows, weights, add_labels,
  ties, tie_mode, initial, and an optional mini-grid path. ``windows_file`` reads
  a simple two-column interval file. Each element is recorded independently.

Normalize and initialize explicitly before the first fit. After a normalization
correction, run another global fit before measuring abundances. Unknown stages
and options are rejected rather than ignored. Structural validation happens
when reading the recipe; library/pixel checks occur during execution.

A weight source can be a native-pixel numeric array, ``{"file": "Mg.txt"}``
for the recommended wavelength-weight table, or
``{"legacy_file": "Mg.mask"}`` for a one-column file on the fitting scale.
``{"dr17": "Mg"}`` loads bundled DR17 element weights, and
``{"dr17": "global"}`` loads its global mask. See :doc:`element-windows`.
Global zero weights exclude pixels from subsequent element fits; positive
magnitudes apply only to that global fit. Continuum fitting uses the spectrum's
bad-pixel mask, not the global weight curve. Ties are explicit and use the
syntax in :doc:`grid`; global and element ties are separate settings.

DR17-style template
-------------------

``examples/pipeline-dr17-style.json`` provides:

1. Native-grid input preparation.
2. External polynomial continuum estimation.
3. Explicit initial stellar parameters.
4. A coarse multi-start fit using the archived DR17 global mask.
5. Model-guided continuum correction.
6. A refined global fit.
7. Fixed-continuum fits for all 24 species attempted in DR17, using archived
   fractional weights and relative C/N/alpha compensation for METALS fits.

**Edit this template before use.** Supply pre-convolved production grid paths,
replace the illustrative label names and starting parameters, and ensure the
starts lie within the library. In particular, ALPHA and LOG10VDOP may have
different names in your library. Grid families/LSF groups can be supplied as
explicit candidates; the runner does not choose fiber groups automatically.
P, Cu, Nd, and C13 were judged unreliable in DR17; their presence represents
historical attempted fits, not recommended measurements.

This is a DR17-style fitting workflow, not an exact reproduction of DR17.
It uses current wrapper defaults (linear interpolation and one FERRE start per
candidate) and external continuum estimators. The historical pipeline used
additional FERRE continuum/optimizer settings, coarse-grid initialization,
selection rules, line-by-line checks, survey quality flags, and calibrations.
Those are not implemented by this recipe runner. Missing survey stages are
not silently described as completed. See the
`DR17 pipeline source <https://github.com/sdss/apogee/blob/e134409dc14b20f69e68a0d4d34b2c1b5056a901/pro/aspcap/doaspcap.pro>`_
and :doc:`initial-estimates` for historical context.

Reference
---------

.. autoclass:: aspcap.PipelineRecipe
   :members: read, example, fit

.. autoclass:: aspcap.RecipeResult
   :members: success, parameters, abundances
