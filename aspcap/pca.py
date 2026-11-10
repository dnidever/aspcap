"""Train and apply shared, blockwise PCA bases for synthetic spectra.

Training accepts batches, so a library need not fit in memory. Existing PCA
coefficients can be recompressed exactly within their original spectral
subspace without expanding every training spectrum into wavelength pixels.
"""

from pathlib import Path
import numpy as np


class PCABasis:
    """A mean spectrum and orthonormal PCA vectors in wavelength blocks.

    Parameters
    ----------
    mean : array-like, shape (npix,)
        Mean physical flux, including any original library mean.
    components : array-like, shape (ncomponents, npix)
        Basis vectors. Each block is independently orthonormal. The same
        component count is used in every block, as in FERRE's PCA format.
    block_lengths : sequence of int
        Number of wavelength pixels in each contiguous compression block.
        Blocks must not span detector gaps.
    """

    def __init__(self, mean, components, block_lengths):
        self.mean = np.array(mean, dtype=float, copy=True)
        self.components = np.array(components, dtype=float, copy=True)
        raw = np.asarray(block_lengths)
        self.block_lengths = np.asarray(block_lengths, dtype=int)
        if (
            raw.ndim != 1
            or len(raw) == 0
            or not np.array_equal(raw, self.block_lengths)
            or np.any(self.block_lengths <= 0)
        ):
            raise ValueError("block_lengths must be positive integers")
        if (
            self.mean.ndim != 1
            or self.components.ndim != 2
            or self.components.shape[1] != len(self.mean)
            or sum(self.block_lengths) != len(self.mean)
            or self.components.shape[0] < 1
            or not np.isfinite(self.mean).all()
            or not np.isfinite(self.components).all()
        ):
            raise ValueError("Finite mean, components and block geometry must agree")
        self.ncomponents = self.components.shape[0]
        self.nblocks = len(self.block_lengths)
        for sl in self.slices:
            vectors = self.components[:, sl]
            if not np.allclose(
                vectors @ vectors.T, np.eye(self.ncomponents), atol=1e-8, rtol=1e-8
            ):
                raise ValueError("Components must be orthonormal within each block")

    @property
    def slices(self):
        """Wavelength slices in stored order (separate from detector segments)."""
        edges = np.r_[0, np.cumsum(self.block_lengths)]
        return tuple(slice(a, b) for a, b in zip(edges[:-1], edges[1:]))

    @classmethod
    def fit(cls, spectra, ncomponents, block_lengths, *, batch_size=1024):
        """Train on a 2D physical-flux array, including memory-mapped arrays.

        Every row has equal weight; no continuum scaling, uncertainty weights,
        masking or automatic row selection is applied. Select the intended
        physical grid region and reserve validation rows before training.
        """
        return cls.from_batches(
            cls._batches(spectra, batch_size), ncomponents, block_lengths
        )

    @classmethod
    def fit_coefficients(
        cls,
        coefficients,
        source_mean,
        source_components,
        ncomponents,
        block_lengths,
        *,
        batch_size=1024
    ):
        """Train from existing, possibly nonorthogonal PCA representations.

        coefficients has shape (nspectra, nblocks, nsource_components).
        source_mean has shape (npix,) and source_components has shape
        (nsource_components, npix). Physical flux in each block is
        source_mean + coefficients @ source_components. A QR factorization
        reduces training to this subspace; results are equivalent to PCA of
        expanded physical flux, up to floating-point roundoff.
        """
        return cls.from_batches(
            cls._batches(coefficients, batch_size),
            ncomponents,
            block_lengths,
            source_mean=source_mean,
            source_components=source_components,
        )

    @staticmethod
    def _batches(array, batch_size):
        if int(batch_size) != batch_size or batch_size < 1:
            raise ValueError("batch_size must be a positive integer")
        for start in range(0, len(array), int(batch_size)):
            yield array[start : start + int(batch_size)]

    @classmethod
    def from_batches(
        cls,
        batches,
        ncomponents,
        block_lengths,
        *,
        source_mean=None,
        source_components=None
    ):
        """Train from an iterable of flux or source-coefficient batches.

        Supply both source arguments for coefficient batches; otherwise each
        batch must be (nrows, npix) physical flux. This entry point permits
        arbitrary selection and streaming from multiple compatible libraries.
        Only selected training rows are consumed. Memory use scales with the
        batch size and per-block covariance, rather than total library size.
        The covariance is accumulated in float64 with centered batch updates.
        """
        lengths = np.asarray(block_lengths, dtype=int)
        if (
            lengths.ndim != 1
            or not np.array_equal(lengths, block_lengths)
            or np.any(lengths <= 0)
            or len(lengths) == 0
            or int(ncomponents) != ncomponents
            or ncomponents < 1
        ):
            raise ValueError(
                "Positive integer block lengths and component count required"
            )
        rank = int(ncomponents)
        npix = int(sum(lengths))
        compressed = source_components is not None
        if compressed != (source_mean is not None):
            raise ValueError("Supply both source_mean and source_components")
        mean0 = (
            np.zeros(npix) if not compressed else np.asarray(source_mean, dtype=float)
        )
        if mean0.shape != (npix,) or not np.isfinite(mean0).all():
            raise ValueError("Invalid source mean")
        if compressed:
            B = np.asarray(source_components, dtype=float)
            if B.ndim != 2 or B.shape[1] != npix or not np.isfinite(B).all():
                raise ValueError("Invalid source components")
        edges = np.r_[0, np.cumsum(lengths)]
        slices = [slice(a, b) for a, b in zip(edges[:-1], edges[1:])]
        qr = [
            (
                np.linalg.qr(B[:, sl].T, mode="reduced")
                if compressed
                else (np.eye(sl.stop - sl.start), None)
            )
            for sl in slices
        ]
        if any(rank > Q.shape[1] for Q, _ in qr):
            raise ValueError("Component count exceeds a source block's subspace")
        means = [np.zeros(Q.shape[1]) for Q, _ in qr]
        moments = [np.zeros((len(mu), len(mu))) for mu in means]
        count = 0
        for batch in batches:
            batch = np.asarray(batch, dtype=float)
            expected = (len(lengths), B.shape[0]) if compressed else (npix,)
            if (
                batch.ndim != len(expected) + 1
                or batch.shape[1:] != expected
                or not np.isfinite(batch).all()
            ):
                raise ValueError(
                    "Training batches must have matching shapes and finite values"
                )
            n = len(batch)
            if n == 0:
                continue
            for k, (sl, (Q, R)) in enumerate(zip(slices, qr)):
                values = batch[:, k] @ R.T if compressed else batch[:, sl]
                mu = values.mean(axis=0)
                centered = values - mu
                delta = mu - means[k]
                moments[k] += centered.T @ centered + np.outer(delta, delta) * (
                    count * n / (count + n)
                )
                means[k] += delta * n / (count + n)
            count += n
        if count < rank + 1:
            raise ValueError("Need at least ncomponents + 1 training spectra")
        result_mean = mean0.copy()
        components = np.empty((rank, npix))
        for k, (sl, (Q, _)) in enumerate(zip(slices, qr)):
            covariance = moments[k] / count
            _, vectors = np.linalg.eigh((covariance + covariance.T) / 2)
            components[:, sl] = vectors[:, -rank:][:, ::-1].T @ Q.T
            result_mean[sl] += means[k] @ Q.T
        return cls(result_mean, components, lengths)

    def project(self, spectra):
        """Return coefficients shaped (..., nblocks, ncomponents) from flux."""
        flux = np.asarray(spectra, dtype=float)
        if (
            flux.ndim < 1
            or flux.shape[-1] != len(self.mean)
            or not np.isfinite(flux).all()
        ):
            raise ValueError("Flux must be finite and match the wavelength pixel count")
        return np.stack(
            [
                (flux[..., sl] - self.mean[sl]) @ self.components[:, sl].T
                for sl in self.slices
            ],
            axis=-2,
        )

    def reconstruct(self, coefficients):
        """Reconstruct physical flux; leading single/batch dimensions persist."""
        c = np.asarray(coefficients, dtype=float)
        if (
            c.ndim < 2
            or c.shape[-2:] != (self.nblocks, self.ncomponents)
            or not np.isfinite(c).all()
        ):
            raise ValueError("Coefficients must be finite and match this basis")
        return np.concatenate(
            [
                self.mean[sl] + c[..., k, :] @ self.components[:, sl]
                for k, sl in enumerate(self.slices)
            ],
            axis=-1,
        )

    def project_coefficients(self, coefficients, source_mean, source_components):
        """Convert old PCA coefficients without expanding spectra to pixels.

        Source and destination must use exactly the same physical wavelength
        pixels, units, normalization and block boundaries. No resampling,
        convolution, continuum correction or error weighting is performed.
        """
        c = np.asarray(coefficients, dtype=float)
        B = np.asarray(source_components, dtype=float)
        mean = np.asarray(source_mean, dtype=float)
        if (
            B.ndim != 2
            or B.shape[1] != len(self.mean)
            or mean.shape != self.mean.shape
            or c.ndim < 2
            or c.shape[-2:] != (self.nblocks, B.shape[0])
            or not all(np.isfinite(x).all() for x in (c, B, mean))
        ):
            raise ValueError(
                "Source representation does not match the destination geometry"
            )
        return np.stack(
            [
                c[..., k, :] @ (B[:, sl] @ self.components[:, sl].T)
                + (mean[sl] - self.mean[sl]) @ self.components[:, sl].T
                for k, sl in enumerate(self.slices)
            ],
            axis=-2,
        )

    def save(self, filename):
        """Save a portable NumPy archive containing only numerical arrays."""
        with Path(filename).open("wb") as stream:
            np.savez_compressed(
                stream,
                mean=self.mean,
                components=self.components,
                block_lengths=self.block_lengths,
            )

    @classmethod
    def load(cls, filename):
        """Load and validate a saved basis without enabling pickle loading."""
        with np.load(filename, allow_pickle=False) as archive:
            return cls(archive["mean"], archive["components"], archive["block_lengths"])
