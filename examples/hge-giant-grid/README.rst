Reduced six-label HGE giant grid
===============================

This builder combines the DR17 Synspec NLTE APOGEE-N LSF-a M- and GK-giant
libraries. It retains the approved 85 Teff/logg nodes from 3000 to 5500 K,
ending at logg=3. All native C, N, alpha, and METALS nodes are retained.

Microturbulence is prescribed by the new fit to raw DR17 LOGG_SPEC::

    vmicro [km/s] = 2.3379795498263993 - 0.43433033061467907 * LOGG

The fit used 166,518 quality-selected, deduplicated RGB stars and gave equal
weight to medians in 0.25-dex raw-gravity bins. It is an empirical prescription,
not an assertion that every star has that velocity. The observed calibration
sample reaches about logg=-0.437; the -0.5 node slightly extends that range.

From a source checkout::

    python examples/hge-giant-grid/build_grid.py ./hge-grid --preflight
    python examples/hge-giant-grid/build_grid.py ./hge-grid --mini --ferre /path/to/ferre
    python examples/hge-giant-grid/build_grid.py ./hge-grid --workers 8

The full build requires network access, several GB of downloads, and roughly
4.6 GiB for its output. It fetches validated HTTP ranges instead of saving
both complete production source grids. Headers are cached; output chunks have
checksums and configuration fingerprints. Repeat the same command to resume.
Use --cache and --proposal to override the header cache and coordinate proposal.
The mini check requires an installed FERRE executable; the full build itself
only manipulates grid coefficients. Set OPENBLAS_NUM_THREADS=1 if appropriate
when using several download workers.

The result is a .hdr / .nodes / .unf bundle. Coordinate rows match spectrum
records, with no record-number column. Excluded combinations have no records.
Retain all three files together. The library contains 673,200 spectra and 7514
modeled wavelength pixels across all three original detector segments.

The source PCA bases differ. The output retains their union, plus their mean
spectrum difference, giving 151 components in each of 12 blocks. This preserves
both source reconstructions without new PCA truncation. This union basis is
not an orthogonal PCA basis: use PCACHI=0 and PCAPROJECT=0, as ASPCAP does by
default. Microturbulence interpolation is linear in the source LOG10VDOP axis.

Use direct file access to avoid loading the full coefficient array::

    from aspcap import Grid
    grid = Grid("hge-grid/hge-giants-dr17-3000-5500-apo-a.hdr", f_access=1)
    spectrum = grid({"C": 0, "N": 0, "O Mg Si S Ca Ti": 0,
                     "METALS": -0.5, "LOGG": 1.5, "TEFF": 4500})

This is a new sparse grid, not a historical DR17 analysis product. Keep the
source bibliographic attribution and generated provenance with redistributed
copies. The large binary grid is not included in the Python distribution.

Validation of the completed build::

    python examples/hge-giant-grid/validate_grid.py ./hge-grid --ferre /path/to/ferre

This checks every output chunk checksum, compares sampled nodes with the
original source records, checks off-node interpolation, confirms excluded
corners, and runs synthetic fits including all six labels. Source comparisons
need network access (and the cached headers). Keep the generated validation.json,
build-metadata.json, build-state.json, and final file checksums with the grid.

Source library description and references:
https://www.sdss4.org/dr17/irspec/apogee-libraries/

The noiseless six-label recovery check uses tighter FERRE convergence controls
(STOPCR=SIMP=1e-9, SCOPE=0.05). Default stopping settings left differences of
about 0.046 dex in C, 0.041 dex in gravity, and 3 K in temperature in this test;
the stricter settings recovered the input labels to better than 0.0002 dex.
These validation controls do not change the package fitting defaults.
