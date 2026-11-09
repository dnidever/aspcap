"""Load, interpolate, and fit a synthetic spectral library with FERRE."""

from __future__ import annotations
from bisect import bisect_left
from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from itertools import product
from math import isfinite
from pathlib import Path
import shutil
from os import PathLike
import subprocess
import tempfile

import numpy as np

from .spectrum import Spec1D

from ._grid_headers import read_library, Segment


class Grid:
    """Inspect, interpolate, and fit a synthetic spectral library with FERRE.

    The header is read when the grid is opened. FERRE is invoked only when
    evaluating models or fitting spectra. Libraries may contain PCA coefficients, multiple spectral segments,
    transposed records, or sparse coordinate tables. FERRE handles their
    payloads; this class exposes native output geometry and label bounds.

    Attributes
    ----------
    path : pathlib.Path
        Resolved path to the ASCII library or binary header.
    labels : tuple of str
        Case-sensitive labels in FERRE axis order.
    lower, step, counts : tuple
        Lower limits, grid spacings, and point counts in label order.
    ranges : tuple of tuple of float
        Inclusive physical bounds, one (lower, upper) pair per label.
    wavelengths : tuple of float
        Native wavelength pixels in the library's wavelength units.
        Logarithmically sampled headers are converted to physical wavelengths.
    vacuum, binary : bool
        Common wavelength convention (None for mixed segments), and whether
        F_FORMAT selects binary data.
    resolution : float
        Header resolution, or a tuple for segments with differing resolutions.
    nlabels, nspectra, npix : int
        Number of labels, stored spectra (including sparse grids), and native
        reconstructed output pixels.
    outwave : tuple of float or None
        Default output sampling for model evaluation.
    segments : tuple of Segment
        Native segment slices, wavelengths, frame, and resolution.
    axis_values : tuple of tuple of float
        Physical label coordinates, including nonuniform sparse axes.
    compressed, sparse : bool
        Whether a supplied library has PCA coefficients or a NODEFILE.
    paths : tuple of pathlib.Path
        Original library paths, passed directly to FERRE.
    executable : str or path-like
        FERRE executable used for evaluation and fitting.
    """

    def __init__(
        self, path, outwave=None, *, executable="ferre", f_format=None, f_access=None
    ):
        """Open one FERRE library or an ordered sequence of spectral modules.

        Parameters
        ----------
        path : path-like or sequence of path-like
            ASCII library, .hdr header, or .unf/.fmt payload path. A sequence
            supplies SYNTHFILE(1), SYNTHFILE(2), ... to FERRE. Keep headers,
            PCA basis records, and coordinate sidecars in their original layout.
        outwave : sequence of float, optional
            Default model-output wavelengths. Each value must lie in an actual
            spectral segment; resampling never bridges gaps between segments.
        executable : path-like, optional
            FERRE executable; resolved when evaluating or fitting.
        f_format : int, optional
            FERRE F_FORMAT (0 ASCII/formatted, 1 unformatted binary). Inferred
            from .hdr/.unf/.fmt when omitted. Applies to every supplied module.
        f_access : int, optional
            FERRE F_ACCESS (0 memory, 1 direct disk access). Defaults to memory,
            except .fmt pairs use formatted direct access. FERRE validates the
            combinations supported by the installed engine.

        Notes
        -----
        PCA reconstruction, transposed payloads, and library validity are handled
        by FERRE. Python reads headers and NODEFILE coordinates for metadata;
        it does not decode the spectral payload. Output arrays always contain
        reconstructed native flux pixels, not PCA coefficients.
        """
        supplied = (path,) if isinstance(path, (str, PathLike)) else tuple(path)
        if not supplied:
            raise ValueError("Supply at least one FERRE library")
        libraries = [read_library(item) for item in supplied]
        first = libraries[0]
        if any(
            item["labels"] != first["labels"] or item["axes"] != first["axes"]
            for item in libraries[1:]
        ):
            raise ValueError("Spectral modules must have matching labels and axes")
        self.paths = tuple(item["path"] for item in libraries)
        self.path = self.paths[0]
        self.labels, self.axis_values = first["labels"], first["axes"]
        self.lower = tuple(axis[0] for axis in self.axis_values)
        self.counts = tuple(len(axis) for axis in self.axis_values)
        self.step = tuple(
            (axis[-1] - axis[0]) / (len(axis) - 1) if len(axis) > 1 else 0.0
            for axis in self.axis_values
        )
        self.ranges = tuple((axis[0], axis[-1]) for axis in self.axis_values)
        self.nspectra = first["nspectra"]
        self.sparse = any(item["sparse"] for item in libraries)
        self._node_tables = tuple(
            item["coordinates"] for item in libraries if item["coordinates"] is not None
        )
        self.compressed = any(item["compressed"] for item in libraries)
        self.headers = tuple(item["header"] for item in libraries)
        self.stored_npix = tuple(item["stored_npix"] for item in libraries)
        self.wavelengths = tuple(w for item in libraries for w in item["wavelengths"])
        self.segments = []
        offset = 0
        for item in libraries:
            for segment in item["segments"]:
                self.segments.append(
                    Segment(
                        segment.start + offset,
                        segment.stop + offset,
                        segment.wavelengths,
                        segment.vacuum,
                        segment.resolution,
                    )
                )
            offset += len(item["wavelengths"])
        self.segments = tuple(self.segments)
        frames = {segment.vacuum for segment in self.segments}
        self.vacuum = frames.pop() if len(frames) == 1 else None
        resolutions = tuple(segment.resolution for segment in self.segments)
        self.resolution = resolutions[0] if len(set(resolutions)) == 1 else resolutions
        formats = {
            int(item.suffix.lower() in (".hdr", ".unf")) for item in map(Path, supplied)
        }
        # Header-only formatted direct-access pairs are also accepted.
        formatted = all(
            p.suffix.lower() == ".hdr"
            and not p.with_suffix(".unf").exists()
            and p.with_suffix(".fmt").exists()
            for p in self.paths
        )
        if f_format is None:
            if len(formats) != 1:
                raise ValueError("Set f_format explicitly for mixed path suffixes")
            f_format = 0 if formatted else formats.pop()
        self.f_format = int(f_format)
        self.f_access = int(
            f_access
            if f_access is not None
            else (formatted or any(Path(p).suffix.lower() == ".fmt" for p in supplied))
        )
        self.binary = self.f_format == 1
        self.executable = executable
        self.nlabels, self.npix = len(self.labels), len(self.wavelengths)
        self.outwave = None if outwave is None else self._output_wave(outwave)

    @classmethod
    def read(cls, path, **options):
        """Open a library using the same path and options as the constructor.

        Pass constructor options such as outwave, executable, f_format, or
        f_access as keywords. Returns an instance
        of the class on which this method is called.
        """
        return cls(path, **options)

    def __repr__(self):
        """Show library size, label ranges, and native wavelength coverage."""
        summary = (
            f"Grid({str(self.path)!r}, {self.nlabels} labels, "
            f"{self.nspectra} spectra, {self.npix} pixels)"
        )
        axes = [
            f"  {name}: {lo:g} to {hi:g} (step {step:g}, {count} points)"
            for name, (lo, hi), step, count in zip(
                self.labels, self.ranges, self.step, self.counts
            )
        ]
        frame = (
            "mixed frames"
            if self.vacuum is None
            else ("vacuum" if self.vacuum else "air")
        )
        wave = f"  Wavelength: {self.wavelengths[0]:g} to {self.wavelengths[-1]:g} ({frame})"
        return "\n".join([summary, *axes, wave])

    def _parameters(self, pars, *, check_bounds=True):
        if isinstance(pars, Mapping):
            if set(pars) != set(self.labels):
                raise ValueError(f"Parameters must match grid labels: {self.labels}")
            pars = [pars[label] for label in self.labels]
        try:
            pars = tuple(float(value) for value in pars)
        except (TypeError, ValueError) as error:
            raise ValueError("Supply one numeric value per grid label") from error
        if len(pars) != self.nlabels:
            raise ValueError(
                f"Expected {self.nlabels} parameters in order {self.labels}"
            )
        for label, value, (low, high) in zip(self.labels, pars, self.ranges):
            if check_bounds and not low <= value <= high:
                raise ValueError(f"Value for {label!r} is outside [{low:g}, {high:g}]")
        return pars

    def contains(self, pars, *, check_model=False):
        """Return whether one parameter point is covered by this grid.

        Parameters
        ----------
        pars : mapping or sequence of float
            Complete parameter mapping using exact grid labels, or one vector
            in ``grid.labels`` order. This method accepts a single point.
        check_model : bool, optional
            Also evaluate the point with FERRE to check model validity. Defaults
            to False: geometry is checked without reading spectra or launching
            the executable. Enable this to detect invalid stored spectra such
            as zero-spectrum exclusions recognized by FERRE in sparse grids.

        Returns
        -------
        bool
            False for nonfinite coordinates, points outside the inclusive label
            bounds, or a missing linear-interpolation corner in a sparse grid.
            Exact sparse nodes need only the neighbors with nonzero weights.
            With check_model=True, a model flagged invalid by FERRE also returns
            False. True is a coverage check, not a scientific quality assessment.

        Raises
        ------
        ValueError
            Labels, vector length, or numeric inputs are malformed, or the
            executable cannot be located when check_model=True.
        RuntimeError
            Model checking fails for an execution or input-file error other
            than FERRE returning an invalid model.

        Notes
        -----
        The geometric check uses the current wrapper's linear interpolation.
        Sparse coordinates are indexed once when the grid is opened. It does
        not scan spectral payloads, validate a FERRE build's format support,
        or test wavelength/LSF compatibility. With check_model=True, FERRE
        evaluates the point using the current model-evaluation settings.
        """
        point = self._parameters(pars, check_bounds=False)
        if any(
            not isfinite(value) or not lo <= value <= hi
            for value, (lo, hi) in zip(point, self.ranges)
        ):
            return False
        if self._node_tables:
            neighbors = []
            for value, axis in zip(point, self.axis_values):
                right = bisect_left(axis, value)
                if right < len(axis) and axis[right] == value:
                    neighbors.append((axis[right],))
                else:
                    neighbors.append((axis[right - 1], axis[right]))
            if any(
                any(corner not in table for table in self._node_tables)
                for corner in product(*neighbors)
            ):
                return False
        if check_model:
            try:
                self(point, as_dict=True)
            except RuntimeError as error:
                if "FERRE returned an invalid model" in str(error):
                    return False
                raise
        return True

    def _sample(self, value):
        for segment in self.segments:
            wave = segment.wavelengths
            if segment.resolution == 0:
                continue
            if wave[0] <= value <= wave[-1]:
                if len(wave) == 1:
                    return segment.start, segment.start, 0.0
                right = max(1, bisect_left(wave, value))
                left = right - 1
                return (
                    segment.start + left,
                    segment.start + right,
                    (value - wave[left]) / (wave[right] - wave[left]),
                )
        raise ValueError(
            "Output wavelengths must fall within a spectral segment, not a gap"
        )

    def _output_wave(self, wave):
        wave = tuple(float(value) for value in wave)
        if not wave or not all(map(isfinite, wave)):
            raise ValueError("Supply finite output wavelengths")
        for value in wave:
            self._sample(value)
        return wave

    def __call__(self, pars, wave=None, *, as_dict=False):
        """Interpolate one model or a batch of models in a single FERRE run.

        Parameters
        ----------
        pars : mapping or sequence
            One complete mapping from labels to values, one vector in
            grid-label order, or a nonempty batch of mappings or vectors.
            Values must be finite and within their label bounds.
        wave : sequence of float, optional
            Output wavelengths within the native range, in the library's
            units and frame. Overrides outwave for this call. If neither is
            supplied, use the native pixels.
        as_dict : bool, optional
            Return the legacy wave/flux dictionary instead of Spec1D objects.
            Defaults to False.

        Returns
        -------
        Spec1D, list of Spec1D, or dict
            A single model returns Spec1D; a batch returns a list of independent
            objects, including a one-row batch. Wavelengths and model flux are
            one-dimensional in native segment order unless output sampling was
            requested. No errors are invented. With as_dict=True, wave is a
            tuple and flux is one tuple or a tuple of spectra for a batch.
            All flux remains on the library's scale.

        Notes
        -----
        FERRE interpolates linearly in the label dimensions. Any wavelength
        resampling is linear and preserves the requested wavelength order.

        Raises
        ------
        ValueError
            The parameters or output wavelengths are invalid, or the
            executable cannot be located.
        RuntimeError
            FERRE failed or returned invalid model output. The run directory
            is retained and included in the error message.
        """
        if isinstance(pars, Mapping):
            batch, rows = False, (pars,)
        else:
            rows = tuple(pars)
            if not rows:
                raise ValueError("Supply at least one set of parameters")
            batch = isinstance(rows[0], Iterable) and not isinstance(
                rows[0], (str, bytes)
            )
            if not batch:
                rows = (rows,)
        rows = tuple(self._parameters(row) for row in rows)
        output_wave = self.outwave if wave is None else self._output_wave(wave)
        if output_wave is None:
            output_wave = self.wavelengths
        models = _execute(self, rows, (), executable=self.executable)
        if output_wave != self.wavelengths:
            samples = [self._sample(value) for value in output_wave]
            models = tuple(
                tuple(model[i] + t * (model[j] - model[i]) for i, j, t in samples)
                for model in models
            )
        if as_dict:
            return {"wave": output_wave, "flux": models if batch else models[0]}
        spectra = []
        for parameters, values in zip(rows, models):
            spectrum = Spec1D(
                values,
                wave=output_wave,
                wavevac=self.vacuum,
                numpix=[len(output_wave)],
                metadata={
                    "parameters": dict(zip(self.labels, parameters)),
                    "grid_paths": [str(path) for path in self.paths],
                    "pca_reconstructed": self.compressed,
                },
            )
            spectra.append(spectrum)
        return spectra if batch else spectra[0]

    def model(self, wave, *pars):
        """Evaluate one model and return its flux tuple.

        Parameters
        ----------
        wave : sequence of float or None
            Output wavelengths as in __call__. None uses the grid default.
        *pars : float
            One positional value per grid label, in the order of labels.

        Returns
        -------
        tuple of float
            Model flux at the requested wavelengths. The wavelength-first
            signature is suitable for optimizers such as scipy.optimize.curve_fit.
        """
        return self(pars, wave=wave, as_dict=True)["flux"]

    def _spectrum_fit_inputs(self, flux, error=None, weights=None):
        """Validate/map a spectrum and combine masks with native weights."""
        spectrum = (
            flux.to_native() if flux._native_indices is not None else flux.flatten()
        )
        if spectrum.wave is None:
            raise ValueError("Spec1D fitting requires wavelengths")
        if (
            self.vacuum is None
            or spectrum.wavevac is None
            or spectrum.wavevac != self.vacuum
        ):
            raise ValueError(
                "Spec1D wavelength frame must match the grid's common air/vacuum frame"
            )
        unit = str(spectrum.wave_unit).strip().lower()
        if unit not in ("angstrom", "angstroms", "aa", "å"):
            raise ValueError(
                "Spec1D wavelengths must use Angstrom units; convert explicitly"
            )
        native = np.asarray(self.wavelengths)
        tolerance = np.maximum(1e-6, np.abs(native) * 1e-10)
        if len(spectrum.wave) != self.npix or not np.all(
            np.isfinite(spectrum.wave) & (np.abs(spectrum.wave - native) <= tolerance)
        ):
            raise ValueError(
                "Spec1D wavelengths must match the native grid; no resampling is applied"
            )
        error = spectrum.err if error is None else np.asarray(error, dtype=float)
        if error is None:
            raise ValueError(
                "Spec1D has no uncertainties; supply spec.err or an explicit error vector"
            )
        error = np.asarray(error, dtype=float)
        if error.shape != (self.npix,):
            raise ValueError(
                "Spec1D fitting errors must have one value per native pixel"
            )
        valid = (
            ~spectrum.mask
            & np.isfinite(spectrum.flux)
            & np.isfinite(error)
            & (error > 0)
        )
        if weights is None:
            weights = np.ones(self.npix)
        else:
            weights = np.asarray(weights, dtype=float)
            if weights.shape != (self.npix,) or not np.all(
                np.isfinite(weights) & (weights >= 0)
            ):
                raise ValueError(
                    "Weights must be finite and nonnegative with one entry per native pixel"
                )
        weights = np.where(valid, weights, 0.0)
        flux = spectrum.flux
        return flux, error, weights

    def fit(
        self,
        flux,
        error=None,
        start=None,
        *,
        free_labels=None,
        weights=None,
        executable=None,
        keep_runs=False,
        run_root=None,
    ):
        """Fit a prepared spectrum on the native grid with FERRE.

        Parameters
        ----------
        flux : Spec1D or sequence of float
            Prepared observation on the library's flux scale. Spec1D inputs
            supply wavelengths, errors, and bad-pixel masks. Saved native pixel
            mappings are restored; otherwise orders concatenate in order, with
            padding removed. Wavelengths must match the native grid. This
            method applies no resampling, velocity, LSF, or continuum correction.
        error : sequence of float, optional
            One-sigma uncertainties in native pixel order. Required for array
            input; defaults to spec.err for Spec1D input. An explicit native
            error vector overrides the object's errors, retaining its mask.
            Invalid Spec1D flux/errors are masked automatically. Array inputs
            require finite flux and positive finite errors on weighted pixels.
        start : mapping or sequence of float, optional
            Full starting parameters keyed by label or in grid-label order.
            Defaults to the midpoint of each axis. Values for frozen labels
            are taken from this vector.
        free_labels : str or sequence of str, optional
            One or more distinct labels to vary. Defaults to every grid label.
        weights : sequence of float, optional
            Nonnegative finite weights with one entry per native pixel. Zero
            excludes a pixel. Defaults to unit weights. More pixels must have
            positive weight than there are free labels.
        executable : str or path-like, optional
            Override this grid's executable for the fit.
        keep_runs : bool, optional
            Keep the FERRE files after success. Failed runs retain their files.
        run_root : str or path-like, optional
            Existing parent directory for the run files. Defaults to the
            system temporary directory.

        Returns
        -------
        FerreFit
            Fitted parameters and conditional errors keyed by label, FERRE's
            reduced chi-square, and model flux on the native grid.

        Raises
        ------
        ValueError
            The inputs are invalid or the executable cannot be located.
        RuntimeError
            FERRE failed or returned invalid fit output. The run directory
            is retained and included in the error message.
        """
        if isinstance(flux, Spec1D):
            flux, error, weights = self._spectrum_fit_inputs(flux, error, weights)
        elif error is None:
            raise ValueError("Supply one-sigma errors with array flux input")
        if start is None:
            start = tuple((lo + hi) / 2 for lo, hi in self.ranges)
        start = self._parameters(start)
        if free_labels is None:
            free_labels = self.labels
        elif isinstance(free_labels, str):
            free_labels = (free_labels,)
        else:
            free_labels = tuple(free_labels)
        if not free_labels or len(set(free_labels)) != len(free_labels):
            raise ValueError("Supply one or more distinct free grid labels")
        axes = tuple(self.labels.index(label) for label in free_labels)
        if weights is None:
            weights = (1,) * self.npix
        if any(len(values) != self.npix for values in (flux, error, weights)):
            raise ValueError("Flux, error, and weights must match the grid pixel count")
        if not all(isfinite(w) and w >= 0 for w in weights):
            raise ValueError("Weights must be finite and nonnegative")
        if sum(w > 0 for w in weights) <= len(axes):
            raise ValueError("The fit requires more weighted pixels than free labels")
        if any(
            w > 0 and not (isfinite(f) and isfinite(e) and e > 0)
            for f, e, w in zip(flux, error, weights)
        ):
            raise ValueError(
                "Weighted pixels require finite flux and positive finite error"
            )
        flux = tuple(f if w > 0 else 1.0 for f, w in zip(flux, weights))
        error = tuple(e if w > 0 else 1e10 for e, w in zip(error, weights))
        return _execute(
            self,
            (start,),
            axes,
            flux=flux,
            error=error,
            weights=weights,
            executable=self.executable if executable is None else executable,
            keep_runs=keep_runs,
            run_root=run_root,
        )

    def fit_abundance(
        self,
        wavelength,
        flux,
        error,
        continuum,
        parameters,
        *,
        star_id,
        species,
        label,
        windows=None,
        weights=None,
        add_labels=(),
        executable=None,
        keep_runs=False,
        run_root=None,
        mask=None,
    ):
        """Fit one element with a fixed continuum and other labels held fixed.

        Parameters
        ----------
        wavelength, flux, error, continuum : sequence of float
            One value per reconstructed native pixel. Wavelengths must match
            this library. Flux and one-sigma errors are divided by the supplied
            finite positive continuum; use ones for already normalized arrays.
        parameters : mapping
            Full global parameter solution, including a starting value for the
            one free abundance label. Other labels remain fixed.
        star_id, species : str
            Spectrum identifier and target element name recorded in the result.
        label : str
            Exact grid label of the abundance dimension to vary.
        windows : sequence of (float, float), optional
            Inclusive intervals restricting the fit. Supply windows or weights.
        weights : sequence of float, optional
            Finite nonnegative values per native pixel, separate from inverse
            variance. Defaults to unit weights inside windows. With both,
            windows restrict the supplied weights. Bad pixels receive zero.
        add_labels : sequence of str, optional
            Distinct fixed logarithmic labels added to the fitted coordinate
            to report [X/H]. Empty when the free label already represents [X/H].
        executable : path-like, optional
            Override this grid's FERRE executable for this fit.
        keep_runs : bool, optional
            Retain successful FERRE run files; defaults to False.
        run_root : path-like, optional
            Existing parent directory for FERRE run files.
        mask : sequence of bool, optional
            One value per native pixel; True excludes a bad pixel.

        Returns
        -------
        result : dict
            Element abundance, raw coordinate, conditional uncertainty,
            reduced chi-square, used-pixel count, and basic flags.
        pixels : iterator
            Rows of wavelength, normalized flux, normalized error, model,
            and window weight.

        Notes
        -----
        Preparation and reporting are shared with the standalone
        :func:`aspcap.fit_abundance`; the numerical fit calls ``self.fit`` with
        one free label. This method does not refit the continuum or select an
        element's grid dimension, windows, or abundance convention for you.
        Errors from input validation and FERRE propagate as in ``Grid.fit``.
        """
        from .abundances import fit_abundance

        return fit_abundance(
            self,
            wavelength,
            flux,
            error,
            continuum,
            parameters,
            star_id=star_id,
            species=species,
            label=label,
            windows=windows,
            weights=weights,
            add_labels=add_labels,
            executable=executable,
            keep_runs=keep_runs,
            run_root=run_root,
            mask=mask,
        )


@dataclass(frozen=True)
class FerreFit:
    """Result of a native-grid FERRE fit.

    Parameters and errors are keyed by grid label. Errors are conditional
    on any frozen labels and are None when FERRE reports no usable positive
    uncertainty. reduced_chi2 is FERRE's reported statistic, and model is the
    best-fit flux tuple in native pixel order.
    """

    parameters: dict[str, float]
    errors: dict[str, float | None]
    reduced_chi2: float
    model: tuple[float, ...]


def run_ferre(
    grid,
    flux,
    error,
    start,
    free_label,
    weights,
    *,
    executable="ferre",
    keep_runs=False,
    run_root=None,
):
    """Run the abundance workflow's single-label fit on the native grid."""
    return grid.fit(
        flux,
        error,
        start,
        free_labels=(free_label,),
        weights=weights,
        executable=executable,
        keep_runs=keep_runs,
        run_root=run_root,
    )


def _execute(
    grid,
    rows,
    axes,
    *,
    flux=None,
    error=None,
    weights=None,
    executable,
    keep_runs=False,
    run_root=None,
):
    ndim = grid.nlabels
    resolved = shutil.which(str(executable))
    if resolved is None:
        raise ValueError(f"FERRE executable not found or not executable: {executable}")
    executable = Path(resolved).resolve()
    directory = Path(tempfile.mkdtemp(prefix="aspcap-", dir=run_root))
    try:
        # Original absolute paths preserve NODEFILE resolution and PCA sidecars.
        # No payload is copied, renamed, decoded, or decompressed in Python.
        synthfiles = "\n".join(
            f" SYNTHFILE({i})='{str(path).replace(chr(39), chr(39) * 2)}'"
            for i, path in enumerate(grid.paths, 1)
        )

        def write_row(name, values):
            row = " ".join(f"{x:.12e}" for x in values)
            (directory / name).write_text(row + "\n", encoding="ascii")

        # FERRE uses a half-open normalized interval [0, 1) internally.
        # Nudge exact upper endpoints inwards far enough for its single precision.
        inputs = [
            " ".join(
                f"{(hi - 2e-6 * (hi - lo) if value == hi and not grid.sparse else value):.12e}"
                for value, (lo, hi) in zip(row, grid.ranges)
            )
            for row in rows
        ]
        (directory / "start.ipf").write_text(
            "".join(f"star{i} {row}\n" for i, row in enumerate(inputs)),
            encoding="ascii",
        )
        fit_files = ""
        if axes:
            write_row("flux.frd", flux)
            write_row("error.err", error)
            write_row("weights.dat", weights)
            indices = " ".join(str(axis + 1) for axis in axes)
            fit_files = f""" INDV={indices}
 FFILE='flux.frd'
 ERFILE='error.err'
 FILTERFILE='weights.dat'
 OPFILE='fit.opf'
"""
        (directory / "input.nml").write_text(
            f"""&LISTA
 NDIM={ndim}
 NOV={len(axes)}
 NOBJ={len(rows)}
{synthfiles}
 F_FORMAT={grid.f_format}
 F_ACCESS={grid.f_access}
 PCACHI=0
 PCAPROJECT=0
 INTER=1
 WINTER=0
 ALGOR=1
 INIT=0
 NRUNS=1
 NTHREADS=1
 CONT=0
 NCONT=0
 OBSCONT=0
 ERRBAR=1
 COVPRINT=1
 PFILE='start.ipf'
{fit_files} OFFILE='model.mdl'
/
""",
            encoding="ascii",
        )
        with (directory / "ferre.log").open("w") as log:
            process = subprocess.run(
                [str(executable)], cwd=directory, stdout=log, stderr=subprocess.STDOUT
            )
        # FERRE 5.1.3 reports a normal completion with STOP 1.
        if process.returncode and not (
            process.returncode == 1
            and "quit: Run ended"
            in (directory / "ferre.log").read_text(errors="replace")
        ):
            process.check_returncode()
        models = tuple(
            tuple(float(x) for x in line.split())
            for line in (directory / "model.mdl").read_text().splitlines()
        )
        if len(models) != len(rows) or any(
            len(model) != grid.npix
            or not all(map(isfinite, model))
            or all(value == -1000 for value in model)
            for model in models
        ):
            raise ValueError("FERRE returned an invalid model")
        result = models
        if axes:
            fields = (directory / "fit.opf").read_text().split()
            if len(fields) < 2 * ndim + 4 or fields[0] != "star0":
                raise ValueError("Unexpected FERRE parameter output")
            params = tuple(float(x) for x in fields[1 : ndim + 1])
            sigma = tuple(float(x) for x in fields[ndim + 1 : 2 * ndim + 1])
            chi2 = 10 ** float(fields[2 * ndim + 3])
            if not all(map(isfinite, (*params, chi2))) or any(
                x == -999.999 for x in params
            ):
                raise ValueError("FERRE returned an invalid solution")
            for i, label in enumerate(grid.labels):
                lo, hi = grid.ranges[i]
                tolerance = max(5e-5 * abs(params[i]), 3e-6 * (hi - lo), 1e-8)
                if i not in axes and abs(params[i] - rows[0][i]) > tolerance:
                    raise ValueError(f"FERRE changed frozen label {label!r}")
            result = FerreFit(
                dict(zip(grid.labels, params)),
                {
                    label: x if isfinite(x) and x > 0 else None
                    for label, x in zip(grid.labels, sigma)
                },
                chi2,
                models[0],
            )
    except (OSError, ValueError, subprocess.SubprocessError, OverflowError) as error:
        raise RuntimeError(
            f"FERRE fit failed. Run files kept in {directory}. {error}"
        ) from error
    if not keep_runs:
        shutil.rmtree(directory)
    return result
