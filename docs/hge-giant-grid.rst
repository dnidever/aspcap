Reduced HGE Giant Grid
======================

The derived HGE grid combines the DR17 Synspec NLTE **M-giant and GK-giant
APOGEE-N LSF-a** libraries into one six-label library. It covers 3000–5500 K
and ends at log g=3. Its free labels are C, N, alpha (named
``O Mg Si S Ca Ti`` in these libraries), METALS, LOGG, and TEFF.

Temperature–gravity selection
-----------------------------

.. list-table:: Retained nodes
   :header-rows: 1

   * - Temperature nodes (K)
     - Minimum log g
     - Maximum log g
   * - 3000, 3100
     - -0.5
     - 0.5
   * - 3200, 3300
     - -0.5
     - 1.0
   * - 3400, 3500
     - -0.5
     - 1.5
   * - 3600, 3700
     - -0.5
     - 2.0
   * - 3800, 3900
     - -0.5
     - 2.5
   * - 4000
     - -0.5
     - 3.0
   * - 4250
     - 0.0
     - 3.0
   * - 4500
     - 0.5
     - 3.0
   * - 4750
     - 1.0
     - 3.0
   * - 5000
     - 1.5
     - 3.0
   * - 5250
     - 2.0
     - 3.0
   * - 5500
     - 2.5
     - 3.0

Gravity spacing is 0.5 dex. Temperature spacing is 100 K through 4000 K and
250 K above that. The NODEFILE retains 85 temperature–gravity combinations
and omits the excluded corners. All original abundance-axis nodes are retained,
giving 673,200 spectra. These are proposed HGE coverage choices, not a calibrated
stellar-selection boundary. Interpolation requires all contributing corners;
higher-order interpolation requires a wider valid stencil.

Microturbulence
---------------

At each stored gravity node, the builder interpolates the original source's
logarithmic microturbulence coordinate at:

.. math::

   v_{\rm micro}\;[\mathrm{km\,s^{-1}}]
   =2.33797955-0.43433033\,\log g_{\rm raw}.

This new empirical line was fitted to the medians of 0.25-dex gravity bins
using 166,518 quality-selected, deduplicated DR17 RGB stars. **LOGG is the raw
FERRE coordinate**, corresponding to catalog LOGG_SPEC, rather than calibrated
gravity. The velocity dimension is absent from the new coordinate table;
no fitting-time microturbulence tie is needed. Temperature and metallicity
subgroups have residual differences, and individual stars scatter about the
relation. The sampled gravity range reached approximately -0.437 to 3.0;
the -0.5 grid node slightly extends it.

PCA representation and use
--------------------------

Both source PCA bases are retained with a mean-offset component: 151 components
per block, over 12 blocks. The shared representation reproduces the source
models without additional PCA truncation. It is a union of basis vectors,
not an orthogonal PCA basis. Set **PCACHI=0 and PCAPROJECT=0**; these are the
ASPCAP wrapper defaults. The file stores all 7514 original modeled wavelength
pixels over the three chips. Observed apStar pixels can be selected with
:func:`aspcap.read_apstar`.

The file bundle is approximately 4.6 GiB. Use direct access:

.. code-block:: python

   from aspcap import Grid

   grid = Grid("hge-grid/hge-giants-dr17-3000-5500-apo-a.hdr", f_access=1)
   spec = grid({"C": 0, "N": 0, "O Mg Si S Ca Ti": 0,
                "METALS": -0.5, "LOGG": 1.5, "TEFF": 4500})

Keep the .hdr, .nodes, and .unf files together. The large binary is not bundled
with ASPCAP or hosted by the package's public catalog. A reproducible builder
and approved node list are in ``examples/hge-giant-grid`` in the repository;
see its README for preflight, small-grid checks, full construction, and resume
commands. Building requires network access; using the completed library does not.

The M/GK source switch is at 4000–4250 K. Source overlap comparisons and
checks of reconstruction, interpolation, exclusions, and synthetic recovery
are recorded with the build. These spot checks do not replace interpolation
and abundance validation over the full label range. Source atmosphere-hole
interpolation and source PCA approximation remain inherited from DR17.

The noiseless six-label recovery check uses tighter FERRE convergence controls
(STOPCR=SIMP=1e-9, SCOPE=0.05). Default stopping settings left differences of
about 0.046 dex in C, 0.041 dex in gravity, and 3 K in temperature in this test;
the stricter settings recovered the input labels to better than 0.0002 dex.
These validation controls do not change the package fitting defaults.
