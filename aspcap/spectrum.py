"""Lightweight multi-order spectra, adapted from Doppler's Spec1D design.

Arrays follow Doppler's convention: one-dimensional for a single order and
``(maximum_pixels, orders)`` otherwise. ``mask=True`` always means unusable.
The implementation is independent of Doppler's fitting, LSF and model loaders.
See LICENSE-DOPPLER for the original source notice.
"""

from copy import deepcopy
import numpy as np


def _orders(value, name):
    if isinstance(value, (list, tuple)) and value and np.ndim(value[0]) > 0:
        arrays = [np.asarray(v) for v in value]
    else:
        array = np.asarray(value)
        if array.ndim == 1:
            arrays = [array]
        elif array.ndim == 2 and array.shape[1] > 0:
            arrays = [array[:, i] for i in range(array.shape[1])]
        else:
            raise ValueError(f"{name} must be one-dimensional or (pixels, orders)")
    if any(a.ndim != 1 for a in arrays):
        raise ValueError(f"{name} must contain one-dimensional orders")
    return arrays


def _pack(arrays, shape, fill, dtype):
    out = np.full(shape, fill, dtype=dtype)
    for i, array in enumerate(arrays):
        if out.ndim == 1:
            out[: len(array)] = array
        else:
            out[: len(array), i] = array
    return out


class Spec1D:
    """Store one or more spectral orders without fitting-model dependencies.

    Parameters
    ----------
    flux : array-like or sequence of arrays
        One order, a ``(pixels, orders)`` array, or unequal-length order arrays.
    err, wave, mask, bitmask : array-like, optional
        Matching uncertainties, wavelengths, bad-pixel flags and integer flags.
        Errors are one-sigma standard deviations. True mask values are bad.
        Nonfinite flux, invalid wavelengths, invalid supplied errors and padded
        samples are automatically masked. No uncertainty is invented if absent.
    numpix : sequence of int, optional
        Valid extent of each order in rectangular input. If omitted, trailing
        invalid wavelengths are treated as padding; interior defects stay masked.
    wavevac : bool or None, optional
        Whether wavelengths are in vacuum. This is metadata; assigning it does
        not convert the coordinates. None denotes an unknown or mixed frame.
        Units/frame conversion is explicitly external.
    wave_unit, flux_unit : object, optional
        Unit metadata; numerical arrays are not converted automatically.
    head, instrument, filename, metadata, lsf : optional
        Observation metadata. LSF information is optional opaque metadata; this
        class does not construct, convolve or interpret an LSF model.

    Notes
    -----
    Order selection and ``copy`` return independent arrays. Public arrays are
    mutable; call ``clear_continuum`` after editing data with a cached continuum.
    This intentionally implements the data interface, not every Doppler method.
    """

    def __init__(
        self,
        flux,
        err=None,
        wave=None,
        mask=None,
        bitmask=None,
        *,
        numpix=None,
        wavevac=True,
        wave_unit="Angstrom",
        flux_unit=None,
        head=None,
        instrument=None,
        filename=None,
        metadata=None,
        lsf=None,
    ):
        fs = _orders(flux, "flux")
        lengths = np.array([len(a) for a in fs], dtype=int)
        shape = (int(max(lengths)),) if len(fs) == 1 else (int(max(lengths)), len(fs))

        def matching(value, name, fill, dtype):
            if value is None:
                return None
            arrays = _orders(value, name)
            if len(arrays) != len(fs) or any(
                len(a) != n for a, n in zip(arrays, lengths)
            ):
                raise ValueError(f"{name} must match the flux order lengths")
            return _pack(arrays, shape, fill, dtype)

        self.flux = _pack(fs, shape, np.nan, float)
        self.err = matching(err, "err", np.inf, float)
        self.wave = matching(wave, "wave", 0.0, float)
        self.mask = matching(mask, "mask", True, bool)
        if self.mask is None:
            self.mask = _pack([np.zeros(n, bool) for n in lengths], shape, True, bool)
        if bitmask is not None and any(
            a.dtype.kind not in "iu" for a in _orders(bitmask, "bitmask")
        ):
            raise ValueError("bitmask must contain integers")
        bitdtype = (
            np.result_type(*[a.dtype for a in _orders(bitmask, "bitmask")])
            if bitmask is not None
            else np.dtype(np.int64)
        )
        if bitdtype.kind not in "iu":
            raise ValueError(
                "Mixed bitmask types cannot be represented without precision loss"
            )
        self.bitmask = matching(bitmask, "bitmask", 0, bitdtype)
        if numpix is None:
            self.numpix = lengths.copy()
            if self.wave is not None:
                for i in range(len(fs)):
                    w = self._column(self.wave, i)
                    good = np.flatnonzero(np.isfinite(w) & (w > 0))
                    self.numpix[i] = good[-1] + 1 if len(good) else 0
        else:
            counts = np.atleast_1d(np.asarray(numpix))
            if (
                counts.shape != lengths.shape
                or not np.all(np.isfinite(counts))
                or np.any(counts != np.floor(counts))
            ):
                raise ValueError("numpix must contain one integer per order")
            if np.any(counts < 0) or np.any(counts > lengths):
                raise ValueError("numpix lies outside the supplied order lengths")
            self.numpix = counts.astype(int)
        self.mask |= ~np.isfinite(self.flux)
        if self.wave is not None:
            self.mask |= ~np.isfinite(self.wave) | (self.wave <= 0)
        if self.err is not None:
            self.mask |= ~np.isfinite(self.err) | (self.err <= 0)
        for i, count in enumerate(self.numpix):
            self._column(self.mask, i)[count:] = True
        self.wavevac = None if wavevac is None else bool(wavevac)
        self.wave_unit, self.flux_unit = wave_unit, flux_unit
        self.head, self.instrument, self.filename = deepcopy(head), instrument, filename
        self.metadata, self.lsf = deepcopy(metadata or {}), deepcopy(lsf)
        self.normalized = False
        self._cont = None
        self._continuum_func = None
        self.continuum_history = []
        self._native_wave = None
        self._native_indices = None

    @staticmethod
    def _column(array, order):
        return array if array.ndim == 1 else array[:, order]

    @property
    def shape(self):
        return self.flux.shape

    @property
    def ndim(self):
        return self.flux.ndim

    @property
    def npix(self):
        return self.shape[0]

    @property
    def norder(self):
        return 1 if self.ndim == 1 else self.shape[1]

    @property
    def size(self):
        """Number of samples within order extents, including masked samples."""
        return int(sum(self.numpix))

    @property
    def valid(self):
        """Current usable pixels, including edits to data and order extents."""
        good = ~self.mask & np.isfinite(self.flux)
        if self.wave is not None:
            good &= np.isfinite(self.wave) & (self.wave > 0)
        if self.err is not None:
            good &= np.isfinite(self.err) & (self.err > 0)
        for i, count in enumerate(self.numpix):
            self._column(good, i)[count:] = False
        return good

    @property
    def snr(self):
        if self.err is None:
            return None
        good = self.valid
        return (
            float(np.median(self.flux[good] / self.err[good]))
            if np.any(good)
            else np.nan
        )

    def __len__(self):
        return self.norder

    def __iter__(self):
        return (self[i] for i in range(self.norder))

    def __array__(self, dtype=None, copy=None):
        """Expose flux to NumPy, following its optional copy contract."""
        if copy is False and dtype is not None and np.dtype(dtype) != self.flux.dtype:
            raise ValueError("Cannot change dtype without copying")
        array = np.asarray(self.flux, dtype=dtype)
        return array.copy() if copy else array

    def __repr__(self):
        return f"Spec1D(norder={self.norder}, numpix={self.numpix.tolist()}, normalized={self.normalized})"

    def copy(self):
        """Deep-copy data, masks, continuum state and metadata."""
        return deepcopy(self)

    def __getitem__(self, index):
        """Copy an order, order slice, or ``spectrum[pixel_slice, order]``.

        A slice of a single-order spectrum selects pixels. All results are
        independent copies, unlike a NumPy view. Native pixel mappings survive.
        """
        pixels = slice(None)
        if isinstance(index, tuple):
            pixels, index = index
            if not isinstance(pixels, slice) or not isinstance(
                index, (int, np.integer)
            ):
                raise IndexError("Use a pixel slice and an integer order")
        elif isinstance(index, slice) and self.norder == 1:
            pixels, index = index, 0
        selected = np.arange(self.norder)[index]
        selected = np.atleast_1d(selected)
        if not len(selected):
            raise IndexError("No orders selected")
        arrays = {}
        for name in ("flux", "wave", "err", "mask", "bitmask"):
            value = getattr(self, name)
            arrays[name] = (
                None
                if value is None
                else [
                    self._column(value, i)[: self.numpix[i]][pixels].copy()
                    for i in selected
                ]
            )
        out = Spec1D(
            **arrays,
            wavevac=self.wavevac,
            wave_unit=self.wave_unit,
            flux_unit=self.flux_unit,
            head=self.head,
            instrument=self.instrument,
            filename=self.filename,
            metadata=self.metadata,
            lsf=self.lsf,
        )
        for name in ("_flux", "_err"):
            value = getattr(self, name, None)
            if value is not None:
                setattr(
                    out,
                    name,
                    _pack(
                        [
                            self._column(value, i)[: self.numpix[i]][pixels]
                            for i in selected
                        ],
                        out.shape,
                        np.nan,
                        float,
                    ),
                )
        out.normalized = self.normalized
        if hasattr(self, "_flux_unit"):
            out._flux_unit = deepcopy(self._flux_unit)
        out._continuum_func = self._continuum_func
        out.continuum_history = deepcopy(self.continuum_history)
        if self._cont is not None:
            out.cont = _pack(
                [
                    self._column(self._cont, i)[: self.numpix[i]][pixels]
                    for i in selected
                ],
                out.shape,
                np.nan,
                float,
            )
        if self._native_indices is not None:
            out._native_wave = self._native_wave.copy()
            out._native_indices = [
                self._native_indices[i][pixels].copy() for i in selected
            ]
        return out

    @property
    def continuum_func(self):
        """Optional callable ``function(spec, **options)``; changing it clears the cache."""
        return self._continuum_func

    @continuum_func.setter
    def continuum_func(self, value):
        if value is not None and not callable(value):
            raise TypeError("continuum_func must be callable or None")
        self._continuum_func = value
        self.clear_continuum()

    def clear_continuum(self):
        """Invalidate a cached estimate after changes to raw data or settings."""
        self._cont = None

    @property
    def cont(self):
        """Cached external continuum, computed lazily when first requested."""
        if self._cont is None:
            from .continuum import fit

            value = self.continuum_func(self) if self.continuum_func else fit(self)
            self.cont = value.continuum if hasattr(value, "continuum") else value
        return self._cont

    @cont.setter
    def cont(self, value):
        if value is None:
            self._cont = None
            return
        array = np.asarray(value, dtype=float)
        if array.shape != self.shape:
            raise ValueError("Continuum must match the flux shape")
        self._cont = array.copy()

    def normalize(self, *, method=None, inplace=True, **options):
        """Apply external normalization; forward options to the estimator.

        Like Doppler, the default updates this object. ``inplace=False`` returns
        a copy. Prefer ``continuum.fit(spec).apply(spec)`` to retain diagnostics.
        Repeating normalization without options is a no-op; changing methods on
        an already normalized spectrum requires the original measurements.
        """
        from .continuum import fit, from_array

        if self.normalized:
            if method is not None or options:
                raise ValueError(
                    "Already normalized; fit a new continuum to the original spectrum"
                )
            return self if inplace else self.copy()
        if method is not None:
            result = fit(self, method=method, **options)
        elif self.continuum_func:
            value = self.continuum_func(self, **options)
            result = value if hasattr(value, "apply") else from_array(self, value)
        elif options or self._cont is None:
            result = fit(self, **options)
        else:
            result = from_array(self, self._cont)
        return result.apply(self, inplace=inplace)

    @classmethod
    def from_spectrum(cls, spectrum):
        """Copy a Doppler-compatible object without importing Doppler.

        Required attribute is flux; wave/err/mask/numpix and metadata are
        optional (wavelength-coordinate methods require wave). Wavelengths and fluxes must already be numerical arrays in the
        desired units/frame. No model initialization or downloads occur.
        """
        if isinstance(spectrum, cls):
            return spectrum.copy()
        kw = {
            name: getattr(spectrum, name, None)
            for name in (
                "err",
                "wave",
                "mask",
                "bitmask",
                "numpix",
                "head",
                "instrument",
                "filename",
                "lsf",
            )
        }
        out = cls(
            spectrum.flux,
            **kw,
            wavevac=getattr(spectrum, "wavevac", True),
            wave_unit=getattr(spectrum, "wave_unit", "Angstrom"),
            flux_unit=getattr(spectrum, "flux_unit", None),
            metadata=getattr(spectrum, "metadata", None),
        )
        out.normalized = bool(getattr(spectrum, "normalized", False))
        cached = getattr(spectrum, "_cont", None)
        if cached is not None:
            out.cont = cached
        for name in ("_flux", "_err", "_flux_unit", "continuum_history"):
            if hasattr(spectrum, name):
                setattr(out, name, deepcopy(getattr(spectrum, name)))
        callback = getattr(spectrum, "continuum_func", None)
        if callable(callback):
            # Preserve the cache as well as the callable; the public setter
            # intentionally invalidates a cache when a user changes methods.
            out._continuum_func = callback
        return out

    @classmethod
    def from_native(
        cls, wave, flux, *, regions, err=None, mask=None, bitmask=None, **metadata
    ):
        """Split a native 1-D grid into orders, preserving original pixel indices.

        ``regions`` is a sequence of inclusive wavelength bounds, one per chip.
        Regions may not overlap. Gaps are retained in the saved native mapping
        and restored as masked samples by ``to_native``; no resampling occurs.
        """
        original = cls(flux, wave=wave, err=err, mask=mask, bitmask=bitmask, **metadata)
        if original.norder != 1 or original.wave is None:
            raise ValueError("from_native requires one-dimensional wave and flux")
        indices, used = [], set()
        for lo, hi in regions:
            if not np.isfinite([lo, hi]).all() or lo > hi:
                raise ValueError("Invalid region bounds")
            selected = np.flatnonzero((original.wave >= lo) & (original.wave <= hi))
            if not len(selected) or used.intersection(selected.tolist()):
                raise ValueError("Regions must be nonempty and non-overlapping")
            used.update(selected.tolist())
            indices.append(selected)
        if not indices:
            raise ValueError("Provide at least one region")
        kw = {
            name: (
                None
                if getattr(original, name) is None
                else [getattr(original, name)[ind] for ind in indices]
            )
            for name in ("flux", "wave", "err", "mask", "bitmask")
        }
        out = cls(**kw, **metadata)
        out._native_wave = original.wave.copy()
        out._native_indices = [ind.copy() for ind in indices]
        return out

    def to_native(self):
        """Restore saved native positions, marking uncovered samples as bad.

        Returns a new single-order spectrum. Without a mapping, use ``flatten``
        for concatenation instead. Gaps get NaN flux, infinite supplied errors,
        and mask=True. Masks/bitmasks are preserved without inventing flag bits.
        """
        if self._native_indices is None:
            raise ValueError("No native pixel mapping; construct with from_native")
        n = len(self._native_wave)
        flux, mask = np.full(n, np.nan), np.ones(n, bool)
        err = None if self.err is None else np.full(n, np.inf)
        bitmask = None if self.bitmask is None else np.zeros(n, self.bitmask.dtype)
        cont = None if self._cont is None else np.full(n, np.nan)
        for i, ind in enumerate(self._native_indices):
            count = self.numpix[i]
            flux[ind] = self._column(self.flux, i)[:count]
            mask[ind] = self._column(self.mask, i)[:count]
            if err is not None:
                err[ind] = self._column(self.err, i)[:count]
            if bitmask is not None:
                bitmask[ind] = self._column(self.bitmask, i)[:count]
            if cont is not None:
                cont[ind] = self._column(self._cont, i)[:count]
        out = Spec1D(
            flux,
            wave=self._native_wave,
            err=err,
            mask=mask,
            bitmask=bitmask,
            wavevac=self.wavevac,
            wave_unit=self.wave_unit,
            flux_unit=self.flux_unit,
            metadata=self.metadata,
            head=self.head,
            instrument=self.instrument,
            filename=self.filename,
            lsf=self.lsf,
        )
        out.normalized = self.normalized
        if hasattr(self, "_flux_unit"):
            out._flux_unit = deepcopy(self._flux_unit)
        out.cont = cont
        for name in ("_flux", "_err"):
            value = getattr(self, name, None)
            if value is not None:
                saved = np.full(n, np.nan)
                for i, ind in enumerate(self._native_indices):
                    saved[ind] = self._column(value, i)[: self.numpix[i]]
                setattr(out, name, saved)
        out.continuum_history = deepcopy(self.continuum_history)
        return out

    def flatten(self):
        """Concatenate order samples; this does not insert detector gaps."""
        kw = {
            name: (
                None
                if getattr(self, name) is None
                else np.concatenate(
                    [
                        self._column(getattr(self, name), i)[:n]
                        for i, n in enumerate(self.numpix)
                    ]
                )
            )
            for name in ("flux", "wave", "err", "mask", "bitmask")
        }
        out = Spec1D(
            **kw,
            wavevac=self.wavevac,
            wave_unit=self.wave_unit,
            flux_unit=self.flux_unit,
            metadata=self.metadata,
            head=self.head,
            instrument=self.instrument,
            filename=self.filename,
            lsf=self.lsf,
        )
        out.normalized = self.normalized
        if hasattr(self, "_flux_unit"):
            out._flux_unit = deepcopy(self._flux_unit)
        for name in ("_flux", "_err"):
            value = getattr(self, name, None)
            if value is not None:
                setattr(
                    out,
                    name,
                    np.concatenate(
                        [self._column(value, i)[:n] for i, n in enumerate(self.numpix)]
                    ),
                )
        if self._cont is not None:
            out.cont = np.concatenate(
                [self._column(self._cont, i)[:n] for i, n in enumerate(self.numpix)]
            )
        out.continuum_history = deepcopy(self.continuum_history)
        if self._native_indices is not None:
            out._native_wave = self._native_wave.copy()
            out._native_indices = [np.concatenate(self._native_indices)]
        return out
