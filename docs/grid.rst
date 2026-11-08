Synthetic spectral grids
========================

See :doc:`ferre` for FERRE's internal processing sequence, control options,
and a link to its manual.

Use :class:`aspcap.Grid` to inspect a FERRE library, generate synthetic spectra,
or fit stellar labels with FERRE. Run these examples from the repository
checkout. They
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

For the historical library families, see the
`DR17 spectral-library documentation
<https://www.sdss4.org/dr17/irspec/apogee-libraries/>`_.

Compression with PCA
--------------------

PCA stores a shared set of basis spectra and a short coefficient vector for
each grid node, rather than every wavelength pixel. FERRE interpolates the
coefficients and reconstructs the model spectrum as needed. DR17 libraries
achieved nearly tenfold compression; the saving depends on how many components
are retained.

Discarding components introduces reconstruction error, which should be checked
in abundance-sensitive features. ASPCAP passes PCA libraries directly to FERRE
for interpolation and reconstruction; ``grid.npix`` and returned flux arrays
refer to the reconstructed spectrum, not the coefficient vector.

See the `DR17 spectral-library documentation
<https://www.sdss4.org/dr17/irspec/apogee-libraries/>`_ and :doc:`ferre`.

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

Starting with ASPCAP 0.2.2, ``Grid`` reads header geometry and passes the original library paths directly
to FERRE. It accepts PCA-compressed and uncompressed libraries, multiple spectral
segments, transposed storage, regular axes, and sparse/nonuniform ``NODEFILE``
coordinate tables. FERRE remains responsible for reading the spectral payload,
reconstructing PCA models, and checking the combinations supported by the
installed executable.

Use an ASCII library, a ``.hdr`` with its ``.unf`` payload, or a formatted
``.hdr``/``.fmt`` pair. Keep the original PCA basis records in the header and
coordinate sidecars beside their library. You may also pass a sequence of paths
to supply separate FERRE ``SYNTHFILE`` modules with matching label axes.
No Python conversion or decompression of the stored spectra is performed.

.. code-block:: python

   grid = Grid("/path/to/compressed-library.hdr")
   print(grid.compressed, grid.npix)
   print(grid.segments)
   model = grid(parameters)  # reconstructed flux in native segment order

   # Separate modules sharing the same stellar labels and axes:
   grid = Grid(["/path/to/blue.hdr", "/path/to/red.hdr"])

   # Explicit FERRE storage/access settings when needed:
   grid = Grid("/path/to/library.hdr", f_format=1, f_access=1)

``f_format`` is FERRE's ``F_FORMAT`` (0 ASCII/formatted, 1 binary), and
``f_access`` is ``F_ACCESS`` (0 memory, 1 direct disk access). They are inferred
for conventional suffixes; a header with only a matching ``.fmt`` defaults to
formatted direct access. Mixed modules use common format/access controls,
as required by FERRE. Binary precision and byte order must match the executable.

``grid.wavelengths`` concatenates reconstructed pixels in file/segment order.
``grid.segments`` exposes each segment's start and exclusive stop indices,
wavelengths, resolution, and air/vacuum convention. ``grid.stored_npix`` records
the stored width for each supplied module, which can be a PCA coefficient count.
Zero-resolution photometric modules retain FERRE's zero-wavelength convention.
Model-output resampling is restricted to real spectral segments and does not
interpolate across detector gaps or photometric entries.

For regular axes, the header gives ``LLIMITS``, ``STEPS``, and ``N_P``.
For sparse grids, ``NODEFILE`` supplies coordinates; ``grid.axis_values`` exposes
the actual sorted values and ``grid.nspectra`` counts stored rows. ``grid.step``
is the mean spacing for an uneven axis, not a claim of uniform sampling.
Missing regions and invalid nodes follow the installed FERRE engine's rules.
Label bounds alone do not guarantee that a sparse interpolation stencil exists.

``LOGW=0``, 1, and 2 describe linear, base-ten logarithmic, and natural-log
wavelengths respectively. Metadata inspection does not validate the entire
payload or guarantee that a particular FERRE build supports every storage,
interpolation, and access combination. Engine errors retain the run files and
log for diagnosis. The small packaged grids remain examples rather than a
complete six-label production library.

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
biases in recovered parameters and abundances. For sparse libraries,
FERRE validates missing nodes and excluded regions; inspect invalid-stencil
failures rather than treating them as usable spectra.

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
``executable="/path/to/ferre"`` when opening the grid. Header inspection reads metadata only. Model evaluation and fitting use FERRE.

Check parameter coverage
------------------------

Use ``grid.contains(parameters)`` to check one point without running FERRE:

.. code-block:: python

   grid.contains({"ALPHA": 0.20, "METALS": -0.25})  # True
   grid.contains({"ALPHA": 0.80, "METALS": -0.25})  # False
   grid.contains((0.20, -0.25))  # positional values follow grid.labels

The check includes label bounds and, for sparse grids, the stored neighbors
needed for linear interpolation. A point can lie inside the outer bounds but
still return False because an interpolation corner is missing. At an exact
stored node, only neighbors with nonzero interpolation weights are required.
Malformed labels or parameter vectors raise ValueError; nonfinite coordinates
return False.

To also evaluate the model and detect invalid spectra recognized by FERRE,
use ``grid.contains(parameters, check_model=True)``. This requires the executable
and can take longer. The default geometric check does not read flux arrays,
so it cannot identify a zero-spectrum exclusion by itself. Coverage does not
establish scientific accuracy or wavelength/LSF compatibility.

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

The returned flux stays on the library's scale; PCA arrays are reconstructed
by FERRE before return. Output wavelengths must use
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

The grid and detector weights are bundled as package data and occupy about
4.3 MB.

To reproduce the grid from SDSS, run
``python examples/download_dr17_mgiant.py dr17-grid`` from the checkout.
The script downloads the PCA header and 125 coefficient records via HTTP
range requests.

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
   :members: read, contains, model, fit
   :special-members: __init__, __call__, __repr__

.. autoclass:: aspcap.ferre.FerreFit
