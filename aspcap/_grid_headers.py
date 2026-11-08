"""Read FERRE geometry without reading or reconstructing its spectral payload."""

from dataclasses import dataclass
from math import exp, isfinite, prod
from pathlib import Path

import f90nml


@dataclass(frozen=True)
class Segment:
    """Native output segment; stop is exclusive and wavelengths retain file order."""

    start: int
    stop: int
    wavelengths: tuple
    vacuum: bool
    resolution: float


def _header(file):
    """Read one SYNTH namelist, honoring quoted slashes and comments."""
    text, quote = [], None
    while True:
        line = file.readline()
        if not line:
            raise ValueError("Unterminated FERRE SYNTH header")
        line = line.decode("ascii")
        i = 0
        while i < len(line):
            c = line[i]
            if quote:
                text.append(c)
                if c == quote:
                    if i + 1 < len(line) and line[i + 1] == quote:
                        text.append(line[i + 1])
                        i += 1
                    else:
                        quote = None
            elif c in "'\"":
                quote = c
                text.append(c)
            elif c == "!":
                break
            elif c == "/":
                text.append(c)
                return dict(f90nml.reads("".join(text))["synth"])
            elif line[i : i + 4].lower() == "&end":
                text.append("/")
                return dict(f90nml.reads("".join(text))["synth"])
            else:
                text.append(c)
            i += 1
        text.append("\n")


def _array(values, name, default=None):
    result = values.get(name, default)
    if result is None:
        raise ValueError(f"Missing FERRE header field {name.upper()}")
    return result if isinstance(result, (list, tuple)) else [result]


def read_library(path):
    """Return geometry for regular, sparse, PCA, and MULTI libraries.

    Spectral data, PCA bases, and their binary layout remain FERRE's concern.
    NODEFILE coordinates are read only to expose physical axes and bounds.
    """
    path = Path(path).resolve()
    if path.suffix.lower() in (".unf", ".fmt"):
        path = path.with_suffix(".hdr")
    with path.open("rb") as file:
        root = _header(file)
        headers = []
        inherited = dict(root)
        for _ in range(int(root.get("multi", 0))):
            inherited.update(_header(file))
            headers.append(dict(inherited))
    if not headers:
        headers = [root]
    geometry = headers[0]
    ndim = int(geometry.get("n_of_dim", geometry.get("ndim", 0)))
    labels = tuple(str(x).strip() for x in _array(geometry, "label"))
    if ndim < 1 or len(labels) != ndim or len(set(labels)) != ndim:
        raise ValueError("FERRE dimensions and unique labels must agree")
    nodefile = geometry.get("nodefile", root.get("nodefile", ""))
    if nodefile:
        nodes = Path(nodefile)
        if not nodes.is_absolute():
            nodes = path.parent / nodes
        rows = []
        for line in nodes.read_text().splitlines():
            line = line.split("#", 1)[0].split("!", 1)[0].replace(",", " ").strip()
            if not line:
                continue
            row = tuple(
                float(x.replace("D", "E").replace("d", "e")) for x in line.split()
            )
            if len(row) != ndim or not all(map(isfinite, row)):
                raise ValueError("NODEFILE rows must contain NDIM finite coordinates")
            rows.append(row)
        if not rows or len(set(rows)) != len(rows):
            raise ValueError("NODEFILE must contain unique coordinate rows")
        axes = tuple(tuple(sorted({row[d] for row in rows})) for d in range(ndim))
        nspectra = len(rows)
    else:
        counts = tuple(int(x) for x in _array(geometry, "n_p"))
        lower = tuple(float(x) for x in _array(geometry, "llimits"))
        steps = tuple(float(x) for x in _array(geometry, "steps"))
        if not len(counts) == len(lower) == len(steps) == ndim:
            raise ValueError("Axis geometry must contain NDIM values")
        if any(n < 1 for n in counts) or not all(map(isfinite, (*lower, *steps))):
            raise ValueError("Invalid FERRE axis geometry")
        if any(n > 1 and step <= 0 for n, step in zip(counts, steps)):
            raise ValueError("Varying axes must have positive steps")
        axes = tuple(
            tuple(lo + i * dx for i in range(n))
            for lo, dx, n in zip(lower, steps, counts)
        )
        nspectra = prod(counts)
    segments, wavelengths = [], []
    for header in headers:
        # MULTI records may omit repeated geometry; f90nml inheritance above
        # follows Fortran namelist persistence between records.
        if tuple(str(x).strip() for x in _array(header, "label")) != labels:
            raise ValueError("Library segments must use the same label order")
        npix = int(header["npix"])
        resolution = float(header.get("resolution", 0))
        wave = _array(header, "wave", [0, 0])
        logw = int(header.get("logw", 0))
        if npix < 1 or logw not in (0, 1, 2) or not isfinite(resolution):
            raise ValueError("Invalid FERRE segment geometry")
        if resolution == 0:
            samples = (0.0,) * npix  # FERRE photometric module convention.
        else:
            if (
                len(wave) != 2
                or not all(map(isfinite, wave))
                or (npix > 1 and wave[1] <= 0)
            ):
                raise ValueError("Invalid FERRE wavelength sampling")
            samples = tuple(float(wave[0]) + i * float(wave[1]) for i in range(npix))
            if logw == 1:
                samples = tuple(10**x for x in samples)
            elif logw == 2:
                samples = tuple(exp(x) for x in samples)
        start = len(wavelengths)
        wavelengths.extend(samples)
        segments.append(
            Segment(
                start,
                len(wavelengths),
                samples,
                bool(header.get("vacuum", 0)),
                resolution,
            )
        )
    npca = tuple(int(x) for x in _array(root, "npca", [0]))
    compressed = any(npca)
    if compressed and sum(npca) != len(wavelengths):
        # A single PCA header's NPIX is the coefficient count, not output pixels.
        if len(headers) == 1 and not root.get("multi", 0):
            h = dict(root)
            h["npix"] = sum(npca)
            if "wave" not in h:
                raise ValueError("PCA library needs output wavelength metadata")
            # Rebuild the sole spectral segment using the expanded pixel count.
            step = float(h["wave"][1])
            start = float(h["wave"][0])
            samples = tuple(start + i * step for i in range(sum(npca)))
            if int(h.get("logw", 0)) == 1:
                samples = tuple(10**x for x in samples)
            elif int(h.get("logw", 0)) == 2:
                samples = tuple(exp(x) for x in samples)
            wavelengths = list(samples)
            segments = [
                Segment(
                    0,
                    len(samples),
                    samples,
                    bool(h.get("vacuum", 0)),
                    float(h.get("resolution", 0)),
                )
            ]
        else:
            raise ValueError("PCA expanded pixels disagree with spectral segments")
    return dict(
        path=path,
        labels=labels,
        axes=axes,
        nspectra=nspectra,
        segments=tuple(segments),
        wavelengths=tuple(wavelengths),
        compressed=compressed,
        sparse=bool(nodefile),
        header=root,
        stored_npix=int(root.get("npix", len(wavelengths))),
    )
