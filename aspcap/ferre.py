"""Read a FERRE library header and run a single prepared-spectrum fit."""


from __future__ import annotations
from dataclasses import dataclass
from math import isfinite
from pathlib import Path
import re
import shutil
import subprocess
import tempfile


@dataclass(frozen=True)
class Grid:
    path: Path
    labels: tuple[str, ...]
    lower: tuple[float, ...]
    step: tuple[float, ...]
    counts: tuple[int, ...]
    wavelengths: tuple[float, ...]
    vacuum: bool
    binary: bool

    @classmethod
    def read(cls, path):
        """Read a single-segment library, using its .hdr for binary grids."""
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
        return cls(path, labels, lower, step, counts, wave, bool(vacuum), binary)


@dataclass(frozen=True)
class FerreFit:
    parameters: dict[str, float]
    errors: dict[str, float | None]
    reduced_chi2: float
    model: tuple[float, ...]


def run_ferre(grid, flux, error, start, free_label, weights, *,
              executable="ferre.x", keep_runs=False, run_root=None):
    """Fit one free label on the native grid, with continuum fitting disabled."""
    axis = grid.labels.index(free_label)
    ndim = len(grid.labels)
    if grid.counts[axis] < 2:
        raise ValueError(f"Grid label {free_label!r} has fewer than two points")
    if set(start) != set(grid.labels):
        raise ValueError(f"Starting parameters must match grid labels: {grid.labels}")
    for label, low, step, count in zip(grid.labels, grid.lower, grid.step, grid.counts):
        if not low <= start[label] <= low + (count - 1) * step:
            raise ValueError(f"Starting value for {label!r} is outside the grid")
    if len(flux) != len(error) or len(flux) != len(weights) or len(flux) != len(grid.wavelengths):
        raise ValueError("Flux, error, and weights must match the grid pixel count")
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

        def write_row(name, values, prefix=""):
            row = " ".join(f"{x:.12e}" for x in values)
            (directory / name).write_text(prefix + row + "\n", encoding="ascii")

        write_row("start.ipf", (start[label] for label in grid.labels), "star ")
        write_row("flux.frd", flux)
        write_row("error.err", error)
        write_row("weights.dat", weights)
        (directory / "input.nml").write_text(f"""&LISTA
 NDIM={ndim}
 NOV=1
 INDV={axis + 1}
 NOBJ=1
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
 FFILE='flux.frd'
 ERFILE='error.err'
 FILTERFILE='weights.dat'
 OPFILE='fit.opf'
 OFFILE='model.mdl'
/
""", encoding="ascii")
        with (directory / "ferre.log").open("w") as log:
            process = subprocess.run([str(executable)], cwd=directory, stdout=log,
                           stderr=subprocess.STDOUT)
        # FERRE 5.1.3 reports a normal completion with STOP 1.
        if process.returncode and not (process.returncode == 1 and
                "quit: Run ended" in (directory / "ferre.log").read_text(errors="replace")):
            process.check_returncode()
        fields = (directory / "fit.opf").read_text().split()
        if len(fields) < 2 * ndim + 4 or fields[0] != "star":
            raise ValueError("Unexpected FERRE parameter output")
        params = tuple(float(x) for x in fields[1:ndim + 1])
        sigma = tuple(float(x) for x in fields[ndim + 1:2 * ndim + 1])
        chi2 = 10 ** float(fields[2 * ndim + 3])
        model = tuple(float(x) for x in (directory / "model.mdl").read_text().split())
        if (len(model) != len(flux) or not all(map(isfinite, (*params, *model, chi2)))
                or any(x == -999.999 for x in params)):
            raise ValueError("FERRE returned an invalid solution or model")
        for i, label in enumerate(grid.labels):
            tolerance = max(5e-5 * abs(params[i]), 1e-8)
            if i != axis and abs(params[i] - start[label]) > tolerance:
                raise ValueError(f"FERRE changed frozen label {label!r}")
        fit = FerreFit(dict(zip(grid.labels, params)),
                       {label: x if isfinite(x) and x > 0 else None
                        for label, x in zip(grid.labels, sigma)}, chi2, model)
    except (OSError, ValueError, subprocess.SubprocessError, OverflowError) as error:
        raise RuntimeError(f"FERRE fit failed. Run files kept in {directory}. {error}") from error
    if not keep_runs:
        shutil.rmtree(directory)
    return fit
