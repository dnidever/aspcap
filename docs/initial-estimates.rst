Initial estimates for global fitting
====================================

An initial estimate supplies the starting value for each grid label before
:doc:`global-fitting`. A good estimate can reduce wasted searches and help
avoid a poor local solution. It is a starting point, not a prior or a fixed
measurement: labels listed as free may move during the fit.

Choosing a starting point and choosing a suitable :doc:`grid` are related but
separate tasks. An accurate temperature estimate cannot rescue a library that
does not cover the star, its composition, or its instrumental profile.

Current ASPCAP interface
------------------------

``Grid.fit(flux, error, start=...)`` accepts a complete mapping keyed by the
exact, case-sensitive ``grid.labels``, or a vector in that order. Supply values
for both free and fixed labels. Labels omitted from ``free_labels`` remain at
the values in ``start``. If ``free_labels`` is omitted, all grid labels vary.

If ``start`` is omitted, the wrapper uses the midpoint of each label's range.
It rejects nonfinite or out-of-range starting values rather than silently
clipping them. The wrapper requests one FERRE run with the supplied starting
parameters (``INIT=0``, ``NRUNS=1``); it does not automatically estimate labels,
try several starts, select among libraries, or run a coarse preliminary stage.
Midpoints are convenient for demonstrations, but are not a general scientific
initialization strategy.

.. code-block:: python

   from aspcap import Grid

   grid = Grid("examples/toy-grid.dat")
   flux = grid({"ALPHA": 0.20, "METALS": -0.15})["flux"]
   error = [0.01] * grid.npix
   fit = grid.fit(
       flux, error,
       start={"ALPHA": 0.0, "METALS": -0.25},
       free_labels=("ALPHA", "METALS"),
   )
   print(fit.parameters)

Run this example from the repository checkout with FERRE installed. It
illustrates the starting-value interface using two toy labels, not a physical
six-label initial estimator.

How the earlier implementations handled initialization
------------------------------------------------------

These are source-level comparisons of the three implementations discussed for
this project. Their configuration and available upstream estimates differ;
they should not be read as a single prescription shared by every survey release.
The links pin the inspected revisions so the comparison remains reproducible.

QASPCAP: FERRE-generated multiple starts
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The inspected ``qaspcap.pro`` writes ``INDINI`` and ``NRUNS`` to configure
FERRE's multiple-start initialization. It begins with two starts per dimension,
then changes selected leading dimensions to one start as dimensionality grows.
``NRUNS`` is the product of these counts. The dimension order comes from the
library; the counts should not be interpreted as fixed choices for particular
physical labels. This path delegates starting-point generation to FERRE rather
than obtaining a per-star atmospheric estimate from a learned predictor.

Source: `QASPCAP qaspcap.pro, initialization block
<https://github.com/callendeprieto/qaspcap/blob/636af1374158acce27cb02dfcf3966efd08eb1f1/qaspcap.pro#L399-L405>`_.

SDSS APOGEE ASPCAP: configurable coarse-to-fine fitting
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The inspected survey IDL driver reads library-specific settings including
``INITPAR``, ``coarse``, ``init``, and ``indini``. For subsequent non-coarse
classes, it selects a previous coarse result through ``aspcap_bestclass`` and
maps its fitted parameters into the next library by label. The FERRE input
writer exposes ``INIT`` and ``INDINI`` and derives the run count from the
initialization counts when supplied. Thus the survey orchestration can use
coarse fits to initialize later fits, with behavior controlled by its library
configuration.

The Python ``ferre.writenml`` helper also exposes ``init``, ``indini``, and
``nruns``. That helper writes FERRE controls; it is not itself a complete
initial-estimate pipeline. The production orchestration is important context
when inspecting only the Python directory originally linked for this project.

Sources: `SDSS doaspcap.pro
<https://github.com/sdss/apogee/blob/e134409dc14b20f69e68a0d4d34b2c1b5056a901/pro/aspcap/doaspcap.pro#L188-L214>`_,
`coarse-result transfer
<https://github.com/sdss/apogee/blob/e134409dc14b20f69e68a0d4d34b2c1b5056a901/pro/aspcap/doaspcap.pro#L333-L344>`_,
`IDL FERRE input writer
<https://github.com/sdss/apogee/blob/e134409dc14b20f69e68a0d4d34b2c1b5056a901/pro/aspcap/writeferre.pro>`_, and
`Python FERRE input writer
<https://github.com/sdss/apogee/blob/e134409dc14b20f69e68a0d4d34b2c1b5056a901/python/apogee/aspcap/ferre.py#L38-L79>`_.

Astra: upstream estimates, coarse search, then refinement
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The inspected default estimate generator tries ApogeeNet results and then
quality-selected Gaia XP estimates from Zhang et al. (2023) for remaining
spectra. It supplies temperature, gravity, and a metallicity estimate, sets
initial alpha, carbon, and nitrogen ratios to zero, and estimates microturbulence
from gravity. An additional specialized path reads existing DR17 parameters.

The coarse-stage planner accepts a configurable estimate generator, identifies
suitable libraries, and clips guesses to their allowed interiors. Its fallback
branches can adjust a finite guess toward a nearby grid or try grid centers
when supplied estimates are invalid. These branches operate on spectra yielded
by the generator; they do not guarantee coverage of spectra it never yields.

The refinement stage selects coarse results using a penalized reduced
chi-square that accounts for quality flags, computes a model-guided continuum,
and passes the selected fitted labels into the stellar-parameter stage.
Initialization is therefore coupled to library selection and normalization,
rather than just a single starting vector.

Sources: `Astra initial.py
<https://github.com/andycasey/astra/blob/5e64e310a98f189091a0d5b1aee34952a9062872/src/astra/pipelines/aspcap/initial.py>`_,
`coarse.py
<https://github.com/andycasey/astra/blob/5e64e310a98f189091a0d5b1aee34952a9062872/src/astra/pipelines/aspcap/coarse.py>`_, and
`stellar_parameters.py
<https://github.com/andycasey/astra/blob/5e64e310a98f189091a0d5b1aee34952a9062872/src/astra/pipelines/aspcap/stellar_parameters.py>`_.

Proposed HGE strategy
---------------------

The following is a design direction, not an implemented automatic estimator:

1. Accept user-supplied estimates with their provenance. Optionally obtain
   temperature, gravity, and metallicity from an upstream analysis, or from
   photometry and other data through an explicitly chosen calibration.
2. Choose compatible libraries and check their label coverage. Use explicit
   defaults for unconstrained abundance ratios, recording that they are defaults.
   A zero ratio can be a starting guess without becoming an abundance prior.
3. When estimates are absent or uncertain, try a small set of physically
   distinct starts or a coarse spectral search rather than relying on one
   midpoint. Include alternative libraries near overlapping grid boundaries.
4. Refine the best acceptable candidates using consistent masks, uncertainties,
   normalization, and instrumental broadening. Compare their residuals and
   solution stability, not only whether FERRE returned successfully.
5. Save the original estimate, its source, any bounds adjustment, all attempted
   starts, the selected library, and convergence or boundary diagnostics.

For the planned HGE analysis, microturbulence follows a relation with log g.
An initial gravity can provide an initial microturbulence value, but the global
fit must keep the relation consistent as gravity changes. See
:doc:`global-fitting` for the current parameter-tying limitation. The historical
implementations' microturbulence handling is not automatically adopted here.

Trying several starts with the existing interface
-------------------------------------------------

You can currently make separate ``Grid.fit`` calls. Continuing the toy example:

.. code-block:: python

   starts = [
       {"ALPHA": -0.2, "METALS": -0.4},
       {"ALPHA": 0.0, "METALS": -0.25},
       {"ALPHA": 0.2, "METALS": -0.1},
   ]
   candidates = [grid.fit(flux, error, start=start) for start in starts]
   usable = [
       fit for fit in candidates
       if fit.reduced_chi2 is not None
       and 0 <= fit.reduced_chi2 < float("inf")
   ]
   best = min(usable, key=lambda fit: fit.reduced_chi2)
   print(best.parameters)

This minimal example assumes every call succeeds and at least one statistic
is usable. Production orchestration should retain failures and handle the case
of no acceptable candidate. Compare statistics directly only when the spectra,
pixel selection, errors, weights, and free-label setup are comparable. Different
libraries or continuum choices may require additional selection checks.
Starting-point agreement is a useful diagnostic, not proof that the model or
reported abundances are accurate.
