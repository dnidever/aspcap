"""Load, interpolate, and fit a synthetic spectral library with FERRE."""

from __future__ import annotations
from bisect import bisect_left
from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from math import isfinite
from pathlib import Path
import shutil
from os import PathLike
import subprocess
import tempfile

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

    def _parameters(self, pars):
        if isinstance(pars, Mapping):
            if set(pars) != set(self.labels):
                raise ValueError(f"Parameters must match grid labels: {self.labels}")
            pars = [pars[label] for label in self.labels]
        pars = tuple(float(value) for value in pars)
        if len(pars) != self.nlabels:
            raise ValueError(
                f"Expected {self.nlabels} parameters in order {self.labels}"
            )
        for label, value, (low, high) in zip(self.labels, pars, self.ranges):
            if not low <= value <= high:
                raise ValueError(f"Value for {label!r} is outside [{low:g}, {high:g}]")
        return pars

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

    def __call__(self, pars, wave=None):
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

        Returns
        -------
        dict
            The wave entry is a wavelength tuple. The flux entry is one flux
            tuple for a single model or a tuple of flux tuples in batch order,
            including for a one-row batch. Flux is on the library's scale.

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
        return {"wave": output_wave, "flux": models if batch else models[0]}

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
        return self(pars, wave=wave)["flux"]

    def fit(
        self,
        flux,
        error,
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
        flux, error : sequence of float
            Observed flux and one-sigma uncertainties, in native pixel order
            and on the library's flux scale. Each sequence must have npix
            entries. Weighted pixels require finite flux and positive finite
            errors. The fit does not apply wavelength or continuum corrections.
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
