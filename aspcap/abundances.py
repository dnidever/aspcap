"""Fit an individual species using the continuum and labels of a global fit."""

from math import isclose, isfinite
import numpy as np


def fit_abundance(
    grid,
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
    ties=None,
    tie_mode="relative",
    add_labels=(),
    executable=None,
    keep_runs=False,
    run_root=None,
    mask=None,
):
    """Return a result row and a model on the native grid.

    This standalone interface is equivalent to grid.fit_abundance(...).
    Both share preparation/reporting and call grid.fit for the numerical fit.
    An omitted executable uses the grid's configured executable.

    Inputs are in the grid's wavelength units and air/vacuum convention.
    The optional mask uses True for unusable pixels.
    Windows are optional inclusive (lower, upper) pairs. Supply windows, weights,
    or both: windows restrict the supplied fractional weights. Weights must be
    finite nonnegative native-pixel values, separate from inverse variance. The optional
    add_labels are coordinates added to report [X/H], using their final values
    if tied. Optional ties and tie_mode follow Grid.fit; none are imposed.
    """
    npix = len(grid.wavelengths)
    badmask = np.zeros(npix, bool) if mask is None else np.asarray(mask, dtype=bool)
    if badmask.shape != (npix,):
        raise ValueError(
            "mask must have one Boolean value per grid pixel (True means bad)"
        )
    if any(len(values) != npix for values in (wavelength, flux, error, continuum)):
        raise ValueError("Spectrum and continuum must have one row per grid pixel")
    if not all(
        isclose(
            a, b, rel_tol=0, abs_tol=1e-4 * (grid.wavelengths[1] - grid.wavelengths[0])
        )
        for a, b in zip(wavelength, grid.wavelengths)
    ):
        raise ValueError("Spectrum wavelengths do not match the native grid")
    if label not in grid.labels:
        raise ValueError(f"Unknown abundance label {label!r}. Available: {grid.labels}")
    if (
        label in add_labels
        or len(set(add_labels)) != len(add_labels)
        or any(name not in parameters or name not in grid.labels for name in add_labels)
    ):
        raise ValueError(
            "Abundance conversion labels must be distinct fixed grid labels"
        )
    if windows is None and weights is None:
        raise ValueError("Provide element windows or weights")
    intervals = None if windows is None else np.asarray(windows, dtype=float)
    if intervals is not None and (
        intervals.ndim != 2
        or intervals.shape[1] != 2
        or not intervals.size
        or not np.all(np.isfinite(intervals))
        or np.any(intervals[:, 0] > intervals[:, 1])
    ):
        raise ValueError("Provide finite wavelength windows (lower, upper)")
    element_weights = (
        np.ones(npix) if weights is None else np.asarray(weights, dtype=float).copy()
    )
    if element_weights.shape != (npix,) or not np.all(
        np.isfinite(element_weights) & (element_weights >= 0)
    ):
        raise ValueError(
            "weights must have one finite nonnegative value per native pixel"
        )
    if intervals is not None:
        selected = np.zeros(npix, dtype=bool)
        for lo, hi in intervals:
            selected |= (np.asarray(wavelength) >= lo) & (np.asarray(wavelength) <= hi)
        element_weights[~selected] = 0

    normalized_flux, normalized_error, weights = [], [], []
    for wave, value, sigma, scale, bad, weight in zip(
        wavelength, flux, error, continuum, badmask, element_weights
    ):
        if not isfinite(scale) or scale <= 0:
            raise ValueError("Global continuum values must be finite and positive")
        value, sigma = value / scale, sigma / scale
        good = not bad and isfinite(value) and isfinite(sigma) and sigma > 0
        normalized_flux.append(value if good else 1.0)
        normalized_error.append(sigma if good else 1e10)
        weights.append(float(weight) if good else 0.0)
    pixels_used = sum(w > 0 for w in weights)
    if pixels_used < 2:
        raise ValueError(
            f"{species}: fewer than two usable pixels in the element windows"
        )

    fit = grid.fit(
        normalized_flux,
        normalized_error,
        parameters,
        free_labels=(label,),
        weights=weights,
        ties=ties,
        tie_mode=tie_mode,
        executable=executable,
        keep_runs=keep_runs,
        run_root=run_root,
    )
    axis = grid.labels.index(label)
    raw = fit.parameters[label]
    coordinates = grid.axis_values[axis]
    low, high = coordinates[0], coordinates[-1]
    low_step = coordinates[1] - low if len(coordinates) > 1 else 0.0
    high_step = high - coordinates[-2] if len(coordinates) > 1 else 0.0
    flags = []
    if raw - low <= 0.125 * low_step or high - raw <= 0.125 * high_step:
        flags.append("GRID_EDGE")
    if fit.errors[label] is None:
        flags.append("NO_UNCERTAINTY")
    conversion_slope = 1.0
    for name in add_labels:
        if name in (ties or {}):
            specification = ties[name]
            coefficients = specification.get("coefficients", specification)
            conversion_slope += float(coefficients.get(label, 0))
    conditional_error = fit.errors[label]
    if conditional_error is not None:
        conditional_error *= abs(conversion_slope)
    result = {
        "star_id": star_id,
        "species": species,
        "grid_label": label,
        "raw_value": raw,
        "abundance_scale": f"[{species}/H]",
        "abundance": raw + sum(fit.parameters[name] for name in add_labels),
        "conditional_error": conditional_error,
        "ferre_reduced_chi2": fit.reduced_chi2,
        "pixels_used": pixels_used,
        "flags": "|".join(flags),
    }
    pixels = zip(wavelength, normalized_flux, normalized_error, fit.model, weights)
    return result, pixels
