"""Published DR17 element windows and FERRE pixel weights."""

from importlib.resources import files

import numpy as np


def dr17_species():
    """Return species with bundled DR17 masks, including historical placeholders.

    Availability does not imply a reliable DR17 abundance measurement. In
    particular Ge, Rb, and Yb were not attempted in DR17, and P, Cu, Nd, and
    C13 were judged unreliable.
    """
    return tuple(
        sorted(
            p.name[:-5]
            for p in files("aspcap").joinpath("data/dr17").iterdir()
            if p.name.endswith(".mask") and p.name != "global.mask"
        )
    )


def load_dr17(species, *, kind="mask", wavelength=None):
    """Load original DR17 weights or the edited window table.

    Parameters
    ----------
    species : str
        Species identifier, e.g. ``Mg``, ``CI``, or ``TiII`` (case insensitive).
    kind : {'mask', 'filt', 'wind'}, optional
        Final FERRE weights (default), original weights, or edited intervals.
    wavelength : array-like, optional
        For weights only, select matching native pixels from the DR17 vacuum
        Angstrom scale. Supports the full apStar array (gaps receive zero) or
        a native grid subset. Finite, positive wavelengths must agree within
        0.002 Angstrom; arbitrary resampling is rejected.

    Returns
    -------
    wave, weights : ndarray, ndarray
        With kind='mask' or 'filt', independent arrays of 7514 native wavelengths
        and dimensionless weights, or arrays on the requested wavelength shape.
    table : ndarray
        With kind='wind', columns are lower/upper vacuum Angstrom bounds and
        edit flag: 0 disables a region, 1 retains it, 2 adds unit weights.
        This table alone does not reproduce the final pixel weights.

    Notes
    -----
    No weight normalization, interpolation, or abundance reliability selection
    is applied. Use the final mask for reproducing the archived DR17 weighting.
    """
    if kind not in ("mask", "filt", "wind"):
        raise ValueError("kind must be 'mask', 'filt', or 'wind'")
    names = {name.lower(): name for name in dr17_species()}
    if kind == "mask":
        names["global"] = "global"
    if not isinstance(species, str) or species.lower() not in names:
        raise ValueError("Unknown DR17 species; see dr17_species()")
    root = files("aspcap").joinpath("data/dr17")
    with root.joinpath(names[species.lower()] + "." + kind).open() as stream:
        values = np.loadtxt(stream, ndmin=2 if kind == "wind" else 1)
    if kind == "wind":
        if wavelength is not None:
            raise ValueError("wavelength selection applies only to pixel weights")
        return values
    with root.joinpath("wave.dat").open() as stream:
        reference = np.loadtxt(stream)
    # Archived integer indices have rounding offsets. Recover registration from
    # wavelengths, whose 0.001-Angstrom precision resolves native pixels safely.
    indices = np.rint((np.log10(reference[:, 0]) - 4.179) / 6e-6).astype(int)
    wave = 10.0 ** (4.179 + indices * 6e-6)
    if wavelength is None:
        return wave, values
    target = np.asarray(wavelength, dtype=float)
    if not np.all(np.isfinite(target) & (target > 0)):
        raise ValueError(
            "wavelength must contain finite positive vacuum Angstrom values"
        )
    pixel = np.rint((np.log10(target) - 4.179) / 6e-6).astype(int)
    if np.any((pixel < 0) | (pixel >= 8575)) or not np.allclose(
        target, 10.0 ** (4.179 + pixel * 6e-6), rtol=0, atol=0.002
    ):
        raise ValueError(
            "wavelength does not match DR17 native pixels; resampling is not automatic"
        )
    full = np.zeros(8575)
    full[indices] = values
    return target.copy(), full[pixel]


def _wave(values):
    """Validate a nonempty one-dimensional wavelength scale."""
    result = np.asarray(values, dtype=float)
    if (
        result.ndim != 1
        or not result.size
        or not np.all(np.isfinite(result) & (result > 0))
    ):
        raise ValueError("wavelength must be a nonempty finite positive 1D array")
    return result


def _weights(values, size):
    result = np.asarray(values, dtype=float)
    if result.shape != (size,) or not np.all(np.isfinite(result) & (result >= 0)):
        raise ValueError(
            "weights must have one finite nonnegative value per wavelength"
        )
    return result


def load_weights(filename, *, wavelength=None, scale=None, atol=0.002):
    """Read a wavelength-weight table or legacy one-column pixel weights.

    Parameters
    ----------
    filename : path-like or readable text stream
        Whitespace-separated numeric file; '#' comments are allowed. Two
        columns mean wavelength and weight, one column means legacy weights.
    wavelength : array-like, optional
        For legacy files, the source pixel wavelengths (required unless
        scale='dr17'). For two-column files, an optional target native scale:
        every target pixel must uniquely match a source pixel within atol.
        Reordering/subsetting is allowed; interpolation and gap filling are not.
    scale : {None, 'dr17'}, optional
        Explicitly identify legacy weights as the 7514-pixel DR17 native scale.
        With this option wavelength may instead select matching target pixels.
    atol : float, optional
        Absolute wavelength matching tolerance, default 0.002 in the table's
        units. No wavelength-unit or air/vacuum conversion is performed.

    Returns
    -------
    wave, weights : ndarray, ndarray
        Independent 1D arrays ready for Grid.fit(weights=weights). Zero excludes
        pixels. These are additional weights, not inverse variances.
    """
    if scale not in (None, "dr17"):
        raise ValueError("scale must be None or 'dr17'")
    if not np.isfinite(atol) or atol < 0:
        raise ValueError("atol must be finite and nonnegative")
    table = np.loadtxt(filename, ndmin=2)
    if table.shape[1] == 1:
        if scale == "dr17":
            source, _ = load_dr17("Mg")
        elif wavelength is not None:
            source = _wave(wavelength)
        else:
            raise ValueError("Legacy weights require wavelength= or scale='dr17'")
        values = _weights(table[:, 0], source.size)
        if scale is None:
            return source.copy(), values.copy()
    elif table.shape[1] == 2:
        if scale is not None:
            raise ValueError("scale applies only to legacy one-column weights")
        source = _wave(table[:, 0])
        values = _weights(table[:, 1], source.size)
    else:
        raise ValueError(
            "Weight files must have one or two columns; use load_windows for intervals"
        )
    if np.unique(source).size != source.size:
        raise ValueError("Source wavelengths must be unique")
    if wavelength is None:
        return source.copy(), values.copy()
    target = _wave(wavelength)
    order = np.argsort(source)
    sorted_wave = source[order]
    left = np.searchsorted(sorted_wave, target - atol, side="left")
    right = np.searchsorted(sorted_wave, target + atol, side="right")
    if np.any(right - left != 1):
        raise ValueError(
            "Each target wavelength must match exactly one source pixel; no automatic resampling"
        )
    return target.copy(), values[order[left]].copy()


def load_windows(filename, *, wavelength=None):
    """Read a two-column lower/upper interval file, optionally rasterizing it.

    Intervals are inclusive, finite, positive, and ordered lower <= upper.
    Overlaps form a union. With wavelength=None return an (N, 2) array suitable
    for Grid.fit_abundance(windows=...). Otherwise return (wave, weights), with
    unit weights inside any interval and zero outside. '#' comments are allowed.
    Use matching wavelength units and air/vacuum conventions throughout.
    Three-column legacy editing files require load_legacy_windows instead.
    """
    intervals = np.loadtxt(filename, ndmin=2)
    if (
        intervals.shape[1] != 2
        or not intervals.size
        or not np.all(np.isfinite(intervals) & (intervals > 0))
        or np.any(intervals[:, 0] > intervals[:, 1])
    ):
        raise ValueError("Windows require finite positive (lower, upper) pairs")
    if wavelength is None:
        return intervals.copy()
    wave = _wave(wavelength)
    weights = np.zeros(wave.size)
    for lo, hi in intervals:
        weights[(wave >= lo) & (wave <= hi)] = 1
    return wave.copy(), weights


def load_legacy_windows(filename, *, wavelength, weights):
    """Apply a legacy three-column .wind edit table to supplied base weights.

    Columns are lower wavelength, upper wavelength, and flag: 0 disables,
    1 leaves the base unchanged, 2 sets unit weight. This follows DR17's strict
    interior bounds and applies all disables before all additions, independently
    of row order. Additional exclusions already present in a final .mask are
    not recoverable from the .wind file alone. Inputs are never modified.
    """
    wave = _wave(wavelength)
    result = _weights(weights, wave.size).copy()
    edits = np.loadtxt(filename, ndmin=2)
    if (
        edits.shape[1] != 3
        or not edits.size
        or not np.all(np.isfinite(edits))
        or np.any(edits[:, :2] <= 0)
        or np.any(edits[:, 0] > edits[:, 1])
        or not np.all(np.isin(edits[:, 2], [0, 1, 2]))
    ):
        raise ValueError("Legacy windows require lower, upper, and flag 0, 1, or 2")
    for flag in (0, 2):
        for lo, hi, _ in edits[edits[:, 2] == flag]:
            result[(wave > lo) & (wave < hi)] = 0 if flag == 0 else 1
    return wave.copy(), result


def save_weights(filename, wavelength, weights):
    """Write the recommended two-column wavelength-weight table.

    Wavelength units and air/vacuum convention must be recorded by the caller;
    the file stores numeric values without conversion. Existing files are
    overwritten. Floating-point precision is preserved for round trips.
    """
    wave = _wave(wavelength)
    values = _weights(weights, wave.size)
    if np.unique(wave).size != wave.size:
        raise ValueError("Source wavelengths must be unique")
    np.savetxt(
        filename,
        np.column_stack((wave, values)),
        fmt="%.17g",
        header="wavelength weight",
    )
