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
   flux = grid({"ALPHA": 0.20, "METALS": -0.15}).flux
   error = [0.01] * grid.npix
   fit = grid.fit(
       flux, error,
       start={"ALPHA": 0.0, "METALS": -0.25},
       free_labels=("ALPHA", "METALS"),
   )
   print(fit.parameters)

Run this example from the repository checkout. It
illustrates the starting-value interface using two toy labels, not a physical
six-label initial estimator.

How the earlier implementations handled initialization
------------------------------------------------------

The inspected versions use three different approaches:

* **QASPCAP — multiple starts.** Configures FERRE's ``INDINI`` and ``NRUNS``
  to search from several starting points.
  `QASPCAP source <https://github.com/callendeprieto/qaspcap/blob/636af1374158acce27cb02dfcf3966efd08eb1f1/qaspcap.pro#L399-L405>`_.

* **DR17 ASPCAP — coarse-to-fine fitting.** Runs preliminary fits, selects a
  preferred coarse solution, and uses its parameters to initialize the next
  library fit. Library configuration controls the initialization.
  `DR17 ASPCAP source <https://github.com/sdss/apogee/blob/e134409dc14b20f69e68a0d4d34b2c1b5056a901/pro/aspcap/doaspcap.pro#L333-L344>`_.

* **Astra — upstream estimates plus refinement.** Starts from ApogeeNet or
  quality-selected Gaia XP estimates, tries suitable grids in a coarse stage,
  then initializes refinement from the preferred coarse result.
  `Estimates <https://github.com/andycasey/astra/blob/5e64e310a98f189091a0d5b1aee34952a9062872/src/astra/pipelines/aspcap/initial.py>`_,
  `coarse search <https://github.com/andycasey/astra/blob/5e64e310a98f189091a0d5b1aee34952a9062872/src/astra/pipelines/aspcap/coarse.py>`_,
  `refinement <https://github.com/andycasey/astra/blob/5e64e310a98f189091a0d5b1aee34952a9062872/src/astra/pipelines/aspcap/stellar_parameters.py>`_.

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
