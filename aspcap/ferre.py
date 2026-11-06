"""Load, interpolate, and fit a synthetic spectral library with FERRE."""


from __future__ import annotations
from bisect import bisect_left
from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from math import isfinite, prod
from pathlib import Path
import re
import shutil
import subprocess
import tempfile


class Grid:
    """Inspect, interpolate, and fit a synthetic spectral library with FERRE.

    The header is read when the grid is opened. FERRE is invoked only when
    evaluating models or fitting spectra. Supported libraries have a single
    segment, positive spectral resolution, and uncompressed spectra.
    Evaluation and fitting require at least two points on every label axis.

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
        Whether wavelengths are in vacuum and whether library data are binary.
    resolution : float
        Spectral resolution recorded in the header.
    nlabels, nspectra, npix : int
        Number of axes, product of the axis sizes, and native pixel count.
    outwave : tuple of float or None
        Default output sampling for model evaluation.
    executable : str or path-like
        FERRE executable used for evaluation and fitting.
    """

    def __init__(self, path, outwave=None, *, executable="ferre.x"):
        """Open a library and validate its metadata.

        Parameters
        ----------
        path : str or path-like
            ASCII library path, or a binary library's .hdr path. The binary
            data must be in a sibling .unf file with the same stem.
        outwave : sequence of float, optional
            Default output wavelengths for evaluation. Values must lie
            within the native wavelength range and use its units and frame.
        executable : str or path-like, optional
            Executable name on PATH or a path to a FERRE executable. Defaults
            to ferre.x. It is resolved when a run starts.

        Raises
        ------
        OSError
            The library header cannot be opened.
        ValueError
            The header or wavelength selection is invalid or unsupported,
            or the binary data file is missing.
        """
        path = Path(path).resolve()
        binary = path.suffix == ".hdr"
        if binary and not path.with_suffix(".unf").is_file():
            raise ValueError(f"Missing binary grid data: {path.with_suffix('.unf')}")
        with path.open("rb") as file:
            lines = []
            for raw in file:
                line = raw.decode("ascii").split("!", 1)[0].strip()
                if line == "/":
                    break
                lines.append(line)
            else:
                raise ValueError(f"Unterminated FERRE header: {path}")
            if file.peek(256).lstrip().startswith(b"&"):
                raise ValueError("Only single-segment FERRE libraries are supported")
        header = " ".join(lines)
        fields = re.split(r"(\b[A-Za-z_][\w]*(?:\(\s*\d+\s*\))?)\s*=", header)
        values = {
            re.sub(r"\s+", "", key).upper(): value.strip().rstrip(",")
            for key, value in zip(fields[1::2], fields[2::2])
        }

        def numbers(key):
            return tuple(float(x.replace("D", "E").replace("d", "e"))
                         for x in values[key].replace(",", " ").split())

        try:
            ndim = int(values.get("N_OF_DIM", values.get("NDIM", "0")))
            labels = tuple(values[f"LABEL({i})"].strip("'\042")
                           for i in range(1, ndim + 1))
            lower, step = numbers("LLIMITS"), numbers("STEPS")
            counts = tuple(int(n) for n in numbers("N_P"))
            start, spacing = numbers("WAVE")
            npix = int(values["NPIX"])
            logw = int(values["LOGW"])
            vacuum = int(values["VACUUM"])
            resolution, = numbers("RESOLUTION")
        except (KeyError, ValueError) as error:
            raise ValueError(f"Invalid or incomplete FERRE header in {path}: {error}") from error
        if (int(values.get("MULTI", "0")) != 0
                or ("NPCA" in values and any(numbers("NPCA")))
                or not isfinite(resolution) or resolution <= 0
                or logw not in (0, 1) or vacuum not in (0, 1)
                or not ndim or len(set(labels)) != ndim
                or not (len(lower) == len(step) == len(counts) == ndim)
                or not all(isfinite(x) for x in (*lower, *step, start, spacing))
                or any(n < 1 for n in counts) or any(s <= 0 for s in step)
                or npix < 2 or spacing <= 0):
            raise ValueError(f"Unsupported FERRE header in {path}")
        wave = tuple(start + i * spacing for i in range(npix))
        if logw:
            wave = tuple(10**x for x in wave)
        self.path = path
        self.labels = labels
        self.lower = lower
        self.step = step
        self.counts = counts
        self.wavelengths = wave
        self.vacuum = bool(vacuum)
        self.binary = binary
        self.resolution = resolution
        self.executable = executable
        self.nlabels = len(labels)
        self.npix = len(wave)
        self.nspectra = prod(counts)
        self.ranges = tuple((lo, lo + dx * (n - 1))
                            for lo, dx, n in zip(lower, step, counts))
        self.outwave = None if outwave is None else self._output_wave(outwave)

    @classmethod
    def read(cls, path, **options):
        """Open a library using the same path and options as the constructor.

        Pass outwave or executable as keyword options. Returns an instance
        of the class on which this method is called.
        """
        return cls(path, **options)

    def __repr__(self):
        """Show library size, label ranges, and native wavelength coverage."""
        summary = (f"Grid({str(self.path)!r}, {self.nlabels} labels, "
                   f"{self.nspectra} spectra, {self.npix} pixels)")
        axes = [f"  {name}: {lo:g} to {hi:g} (step {step:g}, {count} points)"
                for name, (lo, hi), step, count
                in zip(self.labels, self.ranges, self.step, self.counts)]
        frame = "vacuum" if self.vacuum else "air"
        wave = f"  Wavelength: {self.wavelengths[0]:g} to {self.wavelengths[-1]:g} ({frame})"
        return "\n".join([summary, *axes, wave])

    def _parameters(self, pars):
        if isinstance(pars, Mapping):
            if set(pars) != set(self.labels):
                raise ValueError(f"Parameters must match grid labels: {self.labels}")
            pars = [pars[label] for label in self.labels]
        pars = tuple(float(value) for value in pars)
        if len(pars) != self.nlabels:
            raise ValueError(f"Expected {self.nlabels} parameters in order {self.labels}")
        for label, value, (low, high) in zip(self.labels, pars, self.ranges):
            if not low <= value <= high:
                raise ValueError(f"Value for {label!r} is outside [{low:g}, {high:g}]")
        return pars

    def _output_wave(self, wave):
        wave = tuple(float(value) for value in wave)
        if not wave or not all(self.wavelengths[0] <= value <= self.wavelengths[-1]
                               for value in wave):
            raise ValueError("Output wavelengths must fall within the grid")
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
            batch = isinstance(rows[0], Iterable) and not isinstance(rows[0], (str, bytes))
            if not batch:
                rows = (rows,)
        rows = tuple(self._parameters(row) for row in rows)
        output_wave = self.outwave if wave is None else self._output_wave(wave)
        if output_wave is None:
            output_wave = self.wavelengths
        models = _execute(self, rows, (), executable=self.executable)
        if output_wave != self.wavelengths:
            samples = []
            for value in output_wave:
                right = max(1, bisect_left(self.wavelengths, value))
                left = right - 1
                fraction = ((value - self.wavelengths[left]) /
                            (self.wavelengths[right] - self.wavelengths[left]))
                samples.append((left, right, fraction))
            models = tuple(tuple(model[i] + t * (model[j] - model[i])
                                 for i, j, t in samples) for model in models)
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

    def fit(self, flux, error, start=None, *, free_labels=None, weights=None,
            executable=None, keep_runs=False, run_root=None):
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
        if any(w > 0 and not (isfinite(f) and isfinite(e) and e > 0)
               for f, e, w in zip(flux, error, weights)):
            raise ValueError("Weighted pixels require finite flux and positive finite error")
        flux = tuple(f if w > 0 else 1.0 for f, w in zip(flux, weights))
        error = tuple(e if w > 0 else 1e10 for e, w in zip(error, weights))
        return _execute(self, (start,), axes, flux=flux, error=error, weights=weights,
                        executable=self.executable if executable is None else executable,
                        keep_runs=keep_runs, run_root=run_root)


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


def run_ferre(grid, flux, error, start, free_label, weights, *,
              executable="ferre.x", keep_runs=False, run_root=None):
    """Run the abundance workflow's single-label fit on the native grid."""
    return grid.fit(flux, error, start, free_labels=(free_label,), weights=weights,
                    executable=executable, keep_runs=keep_runs, run_root=run_root)


def _execute(grid, rows, axes, *, flux=None, error=None, weights=None,
             executable, keep_runs=False, run_root=None):
    ndim = grid.nlabels
    if any(count < 2 for count in grid.counts):
        raise ValueError("FERRE requires at least two points on every library axis")
    resolved = shutil.which(str(executable))
    if resolved is None:
        raise ValueError(f"FERRE executable not found or not executable: {executable}")
    executable = Path(resolved).resolve()
    directory = Path(tempfile.mkdtemp(prefix="aspcap-", dir=run_root))
    library = directory / ("grid.hdr" if grid.binary else "grid.dat")
    try:
        if grid.binary:
            library.with_suffix(".hdr").symlink_to(grid.path)
            library.with_suffix(".unf").symlink_to(grid.path.with_suffix(".unf"))
        else:
            library.symlink_to(grid.path)

        def write_row(name, values):
            row = " ".join(f"{x:.12e}" for x in values)
            (directory / name).write_text(row + "\n", encoding="ascii")

        # FERRE uses a half-open normalized interval [0, 1) internally.
        # Nudge exact upper endpoints inwards far enough for its single precision.
        inputs = [
            " ".join(f"{(hi - 2e-6 * (hi - lo) if value == hi else value):.12e}"
                     for value, (lo, hi) in zip(row, grid.ranges))
            for row in rows
        ]
        (directory / "start.ipf").write_text(
            "".join(f"star{i} {row}\n" for i, row in enumerate(inputs)), encoding="ascii"
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
        (directory / "input.nml").write_text(f"""&LISTA
 NDIM={ndim}
 NOV={len(axes)}
 NOBJ={len(rows)}
 SYNTHFILE(1)='{library.name}'
 F_FORMAT={int(grid.binary)}
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
""", encoding="ascii")
        with (directory / "ferre.log").open("w") as log:
            process = subprocess.run([str(executable)], cwd=directory, stdout=log,
                           stderr=subprocess.STDOUT)
        # FERRE 5.1.3 reports a normal completion with STOP 1.
        if process.returncode and not (process.returncode == 1 and
                "quit: Run ended" in (directory / "ferre.log").read_text(errors="replace")):
            process.check_returncode()
        models = tuple(tuple(float(x) for x in line.split())
                       for line in (directory / "model.mdl").read_text().splitlines())
        if len(models) != len(rows) or any(
            len(model) != grid.npix or not all(map(isfinite, model)) for model in models
        ):
            raise ValueError("FERRE returned an invalid model")
        result = models
        if axes:
            fields = (directory / "fit.opf").read_text().split()
            if len(fields) < 2 * ndim + 4 or fields[0] != "star0":
                raise ValueError("Unexpected FERRE parameter output")
            params = tuple(float(x) for x in fields[1:ndim + 1])
            sigma = tuple(float(x) for x in fields[ndim + 1:2 * ndim + 1])
            chi2 = 10 ** float(fields[2 * ndim + 3])
            if not all(map(isfinite, (*params, chi2))) or any(x == -999.999 for x in params):
                raise ValueError("FERRE returned an invalid solution")
            for i, label in enumerate(grid.labels):
                lo, hi = grid.ranges[i]
                tolerance = max(5e-5 * abs(params[i]), 3e-6 * (hi - lo), 1e-8)
                if i not in axes and abs(params[i] - rows[0][i]) > tolerance:
                    raise ValueError(f"FERRE changed frozen label {label!r}")
            result = FerreFit(dict(zip(grid.labels, params)),
                              {label: x if isfinite(x) and x > 0 else None
                               for label, x in zip(grid.labels, sigma)}, chi2, models[0])
    except (OSError, ValueError, subprocess.SubprocessError, OverflowError) as error:
        raise RuntimeError(f"FERRE fit failed. Run files kept in {directory}. {error}") from error
    if not keep_runs:
        shutil.rmtree(directory)
    return result
