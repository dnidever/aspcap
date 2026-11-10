PCA Compression
===============

:class:`aspcap.PCABasis` trains a shared, blockwise PCA representation and
converts spectra to and from coefficients. Each wavelength block stores a mean
spectrum and a small set of basis vectors. A model is reconstructed as the mean
plus a weighted sum of those vectors. Fewer components save space but lose
spectral information; validate both flux residuals and abundance fits.

Training from spectra
---------------------

Select the physical parameter region before training and reserve independent
validation spectra. For example, hold out entire abundance combinations rather
than randomly splitting nearly identical neighboring spectra.

.. code-block:: python

   from aspcap import PCABasis

   # training_flux: (number of spectra, number of wavelength pixels).
   # Use blocks appropriate to your wavelength layout; do not cross chip gaps.
   basis = PCABasis.fit(training_flux, ncomponents=30,
                        block_lengths=[500, 500, 500], batch_size=512)
   basis.save("shared-pca.npz")

   coefficients = basis.project(validation_flux)
   reconstructed = basis.reconstruct(coefficients)
   residual = reconstructed - validation_flux

The component count is per block. Training uses equal spectrum and pixel
weights in physical flux, with no automatic normalization, masking, clipping,
or region selection. Input values must be finite and share the same wavelength
sampling and flux scale. Large arrays may be memory mapped. For custom data
loading, ``PCABasis.from_batches`` accepts an iterable of batches.

Recompressing existing PCA models
---------------------------------

Existing basis vectors need not be orthogonal. Training within their complete
spectral subspace avoids expanding every spectrum into wavelength pixels:

.. code-block:: python

   basis = PCABasis.fit_coefficients(
       training_coefficients, old_mean, old_components,
       ncomponents=105, block_lengths=block_lengths, batch_size=1024,
   )
   new_coefficients = basis.project_coefficients(
       old_coefficients, old_mean, old_components,
   )

Coefficient arrays have shape ``(nspectra, nblocks, ncomponents)``; component
arrays have shape ``(ncomponents, npix)``. Mean spectra have shape ``(npix,)``.
The source and target must have identical wavelength pixels and compression
blocks. Recompression inherits errors in the original models.

The repository example ``examples/hge-giant-grid/recompress_ferre.py`` trains
and rewrites the spectrum-major, little-endian float32 FERRE PCA layout used by
the DR17/HGE grid. It preserves coordinate order and wavelength metadata and
refuses to overwrite an existing destination. Other binary layouts require an
appropriate reader; the Python PCA class is independent of file format.

.. code-block:: console

   python examples/hge-giant-grid/recompress_ferre.py train source.hdr shared-pca.npz --components 105
   python examples/hge-giant-grid/recompress_ferre.py compress source.hdr shared-pca.npz compressed-grid

For the original approved HGE grid, ``make_training_rows.py`` reproduces the
adopted abundance-held-out split::

   python examples/hge-giant-grid/make_training_rows.py source.nodes training-records.npy

``--rows training-records.npy`` restricts training to an explicit array of
zero-based source record indices. Omitting it trains on every source row.
The generic example does not infer the HGE selection or a validation split.
Keep the saved basis, training selection, provenance, and validation reports.
The derived HGE library's adopted 105-component basis is described in
:doc:`hge-giant-grid`.

Disk access
-----------

The coefficient file size is approximately
``nspectra * nblocks * ncomponents * 4`` bytes for float32 storage. This is
a disk size; it need not all reside in RAM. ``Grid(header, f_access=1)`` uses
FERRE direct disk access. PCA training and rewriting also process batches.

API
---

.. autoclass:: aspcap.PCABasis
   :members:
