"""Fit an individual species using the continuum and labels of a global fit."""

from math import isclose, isfinite
import numpy as np
from .ferre import run_ferre


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
    windows,
    add_labels=(),
    executable="ferre",
    keep_runs=False,
    run_root=None,
    mask=None,
):
    """Return a result row and a model on the native grid.

    Inputs are in the grid's wavelength units and air/vacuum convention.
    The optional mask uses True for unusable pixels.
    Windows are inclusive (lower, upper) pairs with unit weights. The optional
    add_labels are fixed global-fit coordinates added to report [X/H].
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
    if not windows or any(
        not (isfinite(lo) and isfinite(hi) and lo <= hi) for lo, hi in windows
    ):
        raise ValueError("Provide one or more finite wavelength windows (lower, upper)")

    normalized_flux, normalized_error, weights = [], [], []
    for wave, value, sigma, scale, bad in zip(
        wavelength, flux, error, continuum, badmask
    ):
        if not isfinite(scale) or scale <= 0:
            raise ValueError("Global continuum values must be finite and positive")
        value, sigma = value / scale, sigma / scale
        good = not bad and isfinite(value) and isfinite(sigma) and sigma > 0
        normalized_flux.append(value if good else 1.0)
        normalized_error.append(sigma if good else 1e10)
        weights.append(int(good and any(lo <= wave <= hi for lo, hi in windows)))
    pixels_used = sum(weights)
    if pixels_used < 2:
        raise ValueError(
            f"{species}: fewer than two usable pixels in the element windows"
        )

    fit = run_ferre(
        grid,
        normalized_flux,
        normalized_error,
        parameters,
        label,
        weights,
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
    result = {
        "star_id": star_id,
        "species": species,
        "grid_label": label,
        "raw_value": raw,
        "abundance_scale": f"[{species}/H]",
        "abundance": raw + sum(parameters[name] for name in add_labels),
        "conditional_error": fit.errors[label],
        "ferre_reduced_chi2": fit.reduced_chi2,
        "pixels_used": pixels_used,
        "flags": "|".join(flags),
    }
    pixels = zip(wavelength, normalized_flux, normalized_error, fit.model, weights)
    return result, pixels
