# Changelog

## 0.7.0 — 2026-11-09

- Add `Grid.from_download` and the pinned `hge-giants-test-v1` catalog entry, with anonymous Google Drive downloads, progress, safe resume, publisher checksum verification and offline cache reuse. The loader defaults to direct disk access.
- Add reusable blockwise `PCABasis` training and recompression, including existing nonorthogonal PCA representations and streaming batches.
- Expand grid handling for sparse/nonuniform and PCA libraries, Spec1D inputs and outputs, element weights and parameter ties.
- Add multi-order Spec1D data, external continuum normalization, pipeline recipes and batch processing.
- Bundle DR17 element-window weights, custom weight/window readers and a catalog of production libraries.
- Expand documentation for observed spectra, grids, LSF, FERRE, fitting, pipelines, PCA and installation. Include reproducible HGE grid builders and validation examples.

The multi-gigabyte HGE test grid is downloaded separately; it is not bundled in the Python wheel. Its source is DR17 Synspec NLTE APOGEE-N LSF-a, with a prescribed raw-gravity microturbulence relation and 105 shared PCA components per block. It remains an initial test library.
