Synthetic spectral grids
========================

See :doc:`ferre` for FERRE's internal processing sequence, control options,
and a link to its manual.

Use :class:`aspcap.Grid` to inspect a FERRE library, generate synthetic spectra,
or fit stellar labels with FERRE. Run these examples from the repository
checkout after following the installation steps in :doc:`quickstart`. They
use the bundled synthetic library and vacuum wavelengths in ångströms.

A grid is a library of synthetic spectra indexed by physical labels. Each node
stores a spectrum at one combination of label values. FERRE interpolates
between nodes to evaluate trial models during a fit. The labels and their
abundance conventions belong to the library; the wrapper does not infer them
from a stellar type or element name.

Historical grid families and the HGE plan
--------------------------------------------------

**Temperature and gravity groups.** Earlier ASPCAP analyses maintained several
related libraries covering different stellar parameter ranges: for example, DR17 ASPCAP used separate GK-giant
and M-giant grids, alongside grids for dwarfs and hotter stars. GK and M refer
to stellar temperature regimes, rather than different observed wavelength
bands; their spectra cover the APOGEE spectral region. Splitting these
parameter ranges helped keep individual libraries manageable in memory.
The families overlapped in temperature and other stellar parameters, so a
spectrum could be fitted with more than one candidate grid. ASPCAP then had
to compare those solutions, apply its grid-selection and quality criteria,
and choose the preferred solution and library for that spectrum.

**LSF broadening.** Another distinction was instrumental broadening. Versions of a stellar grid
were pre-convolved with different representative LSFs, allowing the analysis
to select a library appropriate to the observation. These LSF versions were
separate from the division into stellar temperature and gravity ranges.
See :doc:`lsf` for the release-specific LSF groups and the alternative of
internal, on-the-fly convolution in FERRE.

**The HGE plan.** We aim to use **one global stellar grid covering giants** across
the required temperature and gravity range, rather than separate GK-giant
and M-giant global libraries. Its six free labels are temperature, log g,
metallicity, alpha enhancement, carbon, and nitrogen; microturbulence follows
a relation with log g. The final coverage and spacing still need to be chosen
and validated. A single global grid still needs an explicit LSF treatment:
either pre-convolved versions or a suitable high-resolution grid with FERRE
performing convolution for each observation.

**Production grid and test fixtures.** This is the planned production library,
not a claim that the bundled small
APOGEE fixtures already provide that coverage. Element-specific mini-grids,
described in :doc:`element-fitting`, are an additional abundance-fitting
resource and do not change the goal of one global giant-star grid.

For the historical library families, see the
`DR17 spectral-library documentation
<https://www.sdss4.org/dr17/irspec/apogee-libraries/>`_.

Choosing and preparing a grid
-----------------------------

Choose label coverage and spacing appropriate for the stars and the precision
you need. Check temperature, gravity, metallicity, abundance conventions,
microturbulence prescription, wavelength coverage, and instrumental broadening.
See :doc:`global-fitting`, :doc:`element-fitting`, and :doc:`lsf` for how these
choices enter each fitting stage. Observations and models must use compatible
normalization; see :doc:`continuum`.

For the planned six-label HGE global fit, microturbulence follows a relation
with log g. This prescription must be encoded in the library or explicitly
enforced by a future parameter-tying fitting interface. ``Grid`` does not
supply that relation. Element-specific mini-grids add a seventh, dedicated
abundance dimension; their intended use is described in :doc:`element-fitting`.

Supported file formats and limits
---------------------------------

The current ASPCAP reader supports a regular Cartesian grid with constant
positive steps in each label dimension, a single wavelength segment, a positive
header resolution, and uncompressed spectra. Use either an ASCII library
containing its header and spectra, or a ``.hdr`` header beside a matching
``.unf`` binary file. The binary file must match the FERRE build's precision
and byte order; the bundled fixture uses little-endian 32-bit floats.

Required header information includes the dimension count, ordered labels,
``LLIMITS``, ``STEPS``, ``N_P``, ``NPIX``, ``WAVE``, ``LOGW``, ``VACUUM``, and
``RESOLUTION``. Axis values are lower bound plus integer index times step.
``LOGW=0`` gives linear wavelengths and ``LOGW=1`` base-ten logarithmic
wavelengths. ``VACUUM`` declares the wavelength convention; it does not convert
the observation. Reading the header does not validate every stored spectrum.

PCA-compressed libraries, multiple spectral segments, and sparse/nonuniform
coordinate-table grids are not exposed by this ASPCAP reader. The separate
FERRE package's support for additional formats does not automatically extend
this wrapper. The small packaged grids are demonstration fixtures, not a
complete production library for all six HGE labels.

Interpolation and accuracy
--------------------------

Current ASPCAP model evaluations and fits request linear interpolation from
FERRE. The interface does not yet expose quadratic or cubic interpolation
controls. ``Grid.model`` resamples the evaluated flux when requested; this is
separate from interpolation between stellar-label nodes and does not apply
an LSF convolution.

To assess a proposed label spacing, compare interpolated flux with independently
calculated spectra between nodes. Reproducing spectra generated by the same
interpolator tests software consistency, not interpolation accuracy. Start with
one-label tests at several stellar parameters, then test interacting labels
and representative multidimensional patches. Examine spectral residuals and
biases in recovered parameters and abundances. Missing spectra and excluded
regions require explicit format and validity support; zero-filled nodes should
not be assumed to create exclusion boundaries in this wrapper.

Open and inspect a library
--------------------------

.. code-block:: python

   from aspcap import Grid

   grid = Grid("examples/toy-grid.dat")
   print(grid)
   print(grid.labels)      # ('ALPHA', 'METALS')
   print(grid.ranges)      # ((-0.5, 0.5), (-0.5, 0.0))
   print(grid.wavelengths) # Native wavelength pixels

For a binary library, give the path to its ``.hdr`` file and keep the matching
``.unf`` file beside it. If FERRE is outside your PATH, pass
``executable="/path/to/ferre"`` when opening the grid. You can inspect the
header before installing FERRE. Model evaluation and fitting use the executable.

Evaluate synthetic spectra
--------------------------

A mapping is convenient when you know the labels by name. A positional
vector follows ``grid.labels``. Pass a sequence of vectors or mappings to
evaluate a batch with one FERRE invocation.

.. code-block:: python

   model = grid({"ALPHA": 0.20, "METALS": -0.25})
   wavelengths, flux = model["wave"], model["flux"]

   batch = grid([(0.20, -0.25), (0.00, -0.25)])
   first_flux, second_flux = batch["flux"]

   resampled = grid((0.20, -0.25), wave=[15002.5, 15004.5])
   flux_only = grid.model([15002.5, 15004.5], 0.20, -0.25)

The returned flux stays on the library's scale. Output wavelengths must use
the same units and air or vacuum convention as the library and fall within
its coverage. For repeated evaluation at the same wavelengths, pass
``outwave`` when opening the grid. ``model(wave, *pars)`` provides a
wavelength-first callback for an external optimizer.

Fit a spectrum on the native grid
---------------------------------

Supply one flux and one one-sigma uncertainty per native pixel, in the order
of ``grid.wavelengths``. Prepare the spectrum on the library's rest-frame
wavelength grid and flux scale before fitting. The constructor's ``outwave``
setting only affects model evaluation.

This example makes a synthetic observation with small residuals, then fits
both labels from the default grid midpoint.

.. code-block:: python

   observed = list(grid((0.20, -0.15))["flux"])
   observed[3] -= 0.001
   observed[4] += 0.001
   error = [0.01] * grid.npix

   fit = grid.fit(observed, error)
   print(fit.parameters)  # Near ALPHA = 0.20 and METALS = -0.15
   print(fit.errors)
   print(fit.reduced_chi2)
   best_fit_flux = fit.model

To freeze a label, give its value in the full starting vector and list only
the labels FERRE should vary. A zero weight excludes a pixel. For example,
this fit varies ``ALPHA`` only and uses four pixels around the toy feature.

.. code-block:: python

   fit = grid.fit(
       observed, error,
       start={"ALPHA": 0.0, "METALS": -0.15},
       free_labels="ALPHA",
       weights=[int(15002 <= wave <= 15005) for wave in grid.wavelengths],
   )

``fit.parameters`` contains coordinates in the library's label system.
For an element measurement that applies a fixed continuum and converts a
fitted coordinate to an abundance scale, see :doc:`python-api`. Use
``keep_runs=True`` to inspect the FERRE inputs and logs after a successful
fit. A failed FERRE run reports where it kept those files.

Fit the full APOGEE DR17 wavelength grid
----------------------------------------

Use the 5 × 5 × 5 M-giant example to fit metallicity, surface gravity, and
temperature together. Its 125 spectra cover [M/H] from −1.0 to 0.0 dex,
log g from 0.0 to 2.0, and effective temperatures from 3400 to 3800 K.
The generator reconstructs all 7,514 modeled detector pixels from the
same SDSS DR17 Synspec NLTE M-giant library as the smaller fixture below.

Install ASPCAP from PyPI to get the grid and detector
weights as package data. The installed files occupy about 4.3 MB.

.. code-block:: sh

   python -m pip install --upgrade aspcap

To reproduce the grid from SDSS, run
``python examples/download_dr17_mgiant.py dr17-grid`` from the checkout
with NumPy installed. The script downloads the PCA header and 125
coefficient records via HTTP range requests.

The output uses the 8,575 native apStar wavelength pixels, from about
15100.8 to 16999.8 Å in vacuum. Pixels outside the three modeled detector
segments have placeholder flux of one and a weight of zero. Always pass
the provided weights to the fit, combining them with any additional mask
for your observation.

.. code-block:: python

   import json
   from importlib.resources import files

   from aspcap import Grid

   directory = files("aspcap").joinpath("data")
   grid = Grid(directory / "apogee-dr17-mgiant-full.hdr")
   weights = json.loads(
       (directory / "apogee-dr17-mgiant-full-weights.json").read_text()
   )

   truth = {"METALS": -0.32, "LOGG": 1.14, "TEFF": 3625}
   observed = list(grid(truth)["flux"])
   observed[600] += 0.0001
   observed[4200] -= 0.0001
   fit = grid.fit(observed, [0.003] * grid.npix, weights=weights)
   print(fit.parameters)  # Near -0.32, 1.14, and 3625 K

The binary data use little-endian, 32-bit floats as expected by the
standard single-precision FERRE build.

Fit a small APOGEE DR17 M-giant grid
------------------------------------

The checkout includes ``examples/apogee-dr17-mgiant-mini.dat``, a compact
ASCII FERRE library derived from the SDSS DR17 APOGEE Synspec NLTE M-giant
grid. It contains 27 spectra on a 3 × 3 × 3 grid.

* ``METALS`` ([M/H]) runs from −0.50 to 0.00 dex in steps of 0.25 dex.
* ``LOGG`` runs from 0.5 to 1.5 in steps of 0.5 dex.
* ``TEFF`` runs from 3500 to 3700 K in steps of 100 K.

Each spectrum retains 160 native vacuum pixels from about 15736.6 to
15771.2 Å. The remaining library axes are held fixed at
``LOG10VDOP = 0.38021``, ``C = 0``, ``N = 0``, and ``O Mg Si S Ca Ti = 0``.
Use this narrow subset for software tests and demonstrations.

The following example interpolates a synthetic observation between grid
nodes, adds small residuals, then fits all three labels from their grid
midpoints.

.. code-block:: python

   from aspcap import Grid

   grid = Grid("examples/apogee-dr17-mgiant-mini.dat")
   print(grid)
   truth = {"METALS": -0.32, "LOGG": 1.14, "TEFF": 3625}
   observed = list(grid(truth)["flux"])
   observed[60] += 0.0001
   observed[90] -= 0.0001

   fit = grid.fit(observed, [0.003] * grid.npix)
   print(fit.parameters)  # About -0.32, 1.14, and 3630 K

The fixture reconstructs blue-segment pixels 2739 through 2898
(zero-based) from the archived PCA coefficients and basis vectors. It
stores the resulting flux in ASCII, with the original wavelength spacing
and resolution. The source is the APO LSF-a library
``p_apssgM_200921nlte_lsfa_012_075``.
See the SDSS `DR17 APOGEE library description
<https://www.sdss4.org/dr17/irspec/apogee-libraries/>`_ and the
`source archive directory
<https://data.sdss.org/sas/dr17/apogee/spectro/speclib/synth/synspec/marcs/giantisotopes/sgM_200921nlte_lsfa/>`_.

API reference
-------------

.. autoclass:: aspcap.Grid
   :members: read, model, fit
   :special-members: __init__, __call__, __repr__

.. autoclass:: aspcap.ferre.FerreFit
