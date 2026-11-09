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
