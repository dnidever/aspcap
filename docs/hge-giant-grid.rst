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

The adopted library uses **105 shared PCA components per block**, over 12
blocks. The basis was trained only within the approved HGE region, using
537,370 spectra and holding out 135,830 spectra by abundance combination,
including the abundance corners. Both M- and GK-giant spectra share this basis.

Additional reconstruction RMS is 0.0095% of the continuum at the median and
0.077% at the 99th percentile on held-out spectra. The worst held-out RMS is
0.352%; the largest pixel error checked on 6,355 spectra is 2.08%. In the
sampled C/N/Mg/Fe window tests with measurable abundance signals, maximum
abundance shifts were about 0.006 dex. Weak C/N signals at extreme abundances
can produce larger shifts. These tests do not validate all elements or fits.

The coefficient file is **3.16 GiB**, down from 4.54 GiB for the original
151-component union of both source bases. The grid retains all 7514 original
modeled wavelength pixels over the three chips. Observed apStar pixels can
be selected with :func:`aspcap.read_apstar`. Keep PCACHI=0 and PCAPROJECT=0
(the wrapper defaults); PCA-space fitting has not been validated here.
See :doc:`pca-compression` for reusable training and compression code.

Use direct disk access to avoid loading the full coefficient array:

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

Noiseless synthetic validation uses tighter FERRE convergence controls
(STOPCR=SIMP=1e-9, SCOPE=0.05). These validation controls do not change the
package fitting defaults. The original union-grid checks and the final
compressed-grid checks are recorded separately.
