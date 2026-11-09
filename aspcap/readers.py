"""Read prepared CSV spectra and selected apStar FITS spectral rows."""

from pathlib import Path
import csv
import numpy as np
from .spectrum import Spec1D


def read_apstar(filename, *, row=0, grid=None, bad_bits=255):
    """Read one apStar row without velocity shifting or interpolation.

    Parameters
    ----------
    filename : path-like
        apStar FITS with flux, one-sigma error, and integer pixel flags in
        HDUs 1, 2, and 3. Compressed FITS is supported through Astropy.
    row : int, optional
        Zero-based spectral row (default 0, combined). Specify actual FITS
        row indices for alternate combinations or visits; layouts vary.
    grid : Grid, optional
        Select matching native vacuum Angstrom pixels, including subsets.
        Every model pixel must uniquely match an observed pixel within 0.002
        Angstrom. No interpolation or air/vacuum conversion is performed.
    bad_bits : int, optional
        Bit selection to mark unusable pixels, default 255 (bits 0 through 7).
        Original integer flags are preserved, including warning bits.
        Zero applies only finite-flux/positive-error checks.

    Returns
    -------
    Spec1D
        Flat native-pixel spectrum with header/row/file identity and flags.
        The reader does not interpret LSF arrays or perform continuum fitting.
    """
    from astropy.io import fits

    if isinstance(row, bool) or not isinstance(row, (int, np.integer)) or row < 0:
        raise ValueError("row must be a nonnegative integer")
    if (
        isinstance(bad_bits, bool)
        or not isinstance(bad_bits, (int, np.integer))
        or bad_bits < 0
    ):
        raise ValueError("bad_bits must be a nonnegative integer")
    with fits.open(filename, memmap=False) as hdus:
        if len(hdus) < 4:
            raise ValueError("apStar requires flux, error, and flag HDUs 1, 2, 3")
        arrays = [np.asarray(hdus[i].data) for i in (1, 2, 3)]
        if arrays[0].ndim not in (1, 2) or any(
            a.shape != arrays[0].shape for a in arrays
        ):
            raise ValueError(
                "apStar flux, error, and flags must have matching 1D/2D shapes"
            )
        if not np.issubdtype(arrays[2].dtype, np.integer):
            raise ValueError("apStar pixel flags must be integers")
        if arrays[0].ndim == 1:
            if row != 0:
                raise ValueError("Single-row spectrum only has row 0")
            flux, error, flags = [a.copy() for a in arrays]
        else:
            if row >= arrays[0].shape[0]:
                raise ValueError("Requested row is outside the apStar image")
            flux, error, flags = [a[row].copy() for a in arrays]
        primary, header = hdus[0].header.copy(), hdus[1].header.copy()
    if len(flux) != 8575:
        raise ValueError(
            "Expected an 8575-pixel apStar row; other FITS products are not supported"
        )
    crval = float(header["CRVAL1"])
    cdelt = float(header["CDELT1"])
    crpix = float(header.get("CRPIX1", 1))
    if (
        not np.all(np.isfinite([crval, cdelt, crpix]))
        or cdelt <= 0
        or not 3 < crval < 6
    ):
        raise ValueError("Expected an apStar log10 wavelength header")
    wave = 10.0 ** (crval + (np.arange(len(flux)) + 1 - crpix) * cdelt)
    if grid is not None:
        if grid.vacuum is not True:
            raise ValueError("apStar reader requires a vacuum-wavelength grid")
        target = np.asarray(grid.wavelengths, dtype=float)
        left = np.searchsorted(wave, target - 0.002)
        right = np.searchsorted(wave, target + 0.002, side="right")
        if np.any(right - left != 1):
            raise ValueError(
                "Grid wavelengths do not match apStar pixels; no automatic resampling"
            )
        flux, error, flags, wave = flux[left], error[left], flags[left], wave[left]
    mask = (
        ((flags & bad_bits) != 0)
        | ~np.isfinite(flux)
        | ~np.isfinite(error)
        | (error <= 0)
    )
    identity = str(primary.get("APOGEE_ID", primary.get("OBJECT", Path(filename).stem)))
    return Spec1D(
        flux,
        wave=wave,
        err=error,
        mask=mask,
        bitmask=flags,
        wavevac=True,
        filename=str(filename),
        head=primary,
        metadata={
            "star_id": f"{identity}:row{row}",
            "fits_row": int(row),
            "source_file": str(filename),
            "bad_bits": int(bad_bits),
        },
        instrument="APOGEE",
        wave_unit="Angstrom",
    )


def read_spectrum(filename, *, grid=None, row=0, bad_bits=255, wavevac=True):
    """Read apStar FITS or a prepared CSV (wavelength,flux,error[,mask])."""
    suffixes = "".join(Path(filename).suffixes).lower()
    if suffixes.endswith((".fits", ".fit", ".fits.gz", ".fit.gz", ".fits.fz")):
        if not wavevac:
            raise ValueError("apStar wavelengths are vacuum wavelengths")
        return read_apstar(filename, row=row, grid=grid, bad_bits=bad_bits)
    with Path(filename).open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    return Spec1D(
        [float(r["flux"]) for r in rows],
        wave=[float(r["wavelength"]) for r in rows],
        err=[float(r["error"]) for r in rows],
        mask=[r.get("mask", "0").lower() in ("1", "true") for r in rows],
        wavevac=wavevac,
        filename=str(filename),
        metadata={"star_id": Path(filename).stem},
    )
