"""External, order-by-order continuum estimation. No FERRE calls are made.

Methods are deliberately explicit: polynomial/spline fits, a running median,
a binned percentile envelope, or a median correction to a data/model ratio.
These are configurable estimators, not claims of exact DR17 reproduction.
"""

from copy import deepcopy
from dataclasses import dataclass, field
import numpy as np
from .spectrum import Spec1D


def _aligned(a, b):
    return (
        a is None
        and b is None
        or (a is not None and b is not None and np.array_equal(a, b, equal_nan=True))
    )


@dataclass
class ContinuumResult:
    """Continuum estimate and diagnostics, independent of a spectrum object.

    continuum : ndarray
        Estimated divisor in the input flux units; dimensionless for model_ratio.
        Invalid orders and padding contain NaN.
    fit_mask : ndarray of bool
        True where an original sample contributed to the estimate (opposite to
        the spectrum's bad-pixel mask). Interpolated bad pixels are not fit data.
    flags : tuple of tuples
        Per-order quality flags. Empty tuples indicate successful estimation.
    settings, parameters : dict, tuple
        Method settings and per-order coefficients/support diagnostics.

    Pixel errors rescaled by ``apply`` remain conditional on this continuum;
    continuum-estimation uncertainty is not added to the error array.
    """

    continuum: np.ndarray
    fit_mask: np.ndarray
    flags: tuple
    settings: dict
    parameters: tuple
    _wave: object = field(repr=False)
    _numpix: np.ndarray = field(repr=False)

    def apply(self, spectrum, *, inplace=False):
        """Divide flux/errors by the same scale; copy by default.

        Accepts Doppler-compatible objects and preserves their type. Invalid
        scales remain masked with NaN flux/infinite errors. This method does not
        change wavelength, bitmasks, or uncertainty floors. Repeated explicit
        applications compose cached continuum factors and retain the raw data.
        """
        data = Spec1D.from_spectrum(spectrum)
        if (
            data.shape != self.continuum.shape
            or not _aligned(data.wave, self._wave)
            or not np.array_equal(data.numpix, self._numpix)
        ):
            raise ValueError(
                "Continuum result does not match the spectrum wavelengths/order layout"
            )
        out = spectrum if inplace else deepcopy(spectrum)
        scale = self.continuum
        good_scale = np.isfinite(scale) & (scale > 0)
        good = data.valid & good_scale
        old = getattr(spectrum, "_cont", None) if data.normalized else None
        if data.normalized and old is None:
            raise ValueError(
                "Normalized input has no cached continuum for composing corrections"
            )
        flux = np.full(data.shape, np.nan)
        np.divide(data.flux, scale, out=flux, where=good)
        error = None
        if data.err is not None:
            error = np.full(data.shape, np.inf)
            np.divide(data.err, scale, out=error, where=good)
        if not data.normalized:
            out._flux = data.flux.copy()
            out._err = None if data.err is None else data.err.copy()
            out._flux_unit = data.flux_unit
        out.flux, out.err, out.mask = flux, error, ~good
        out.cont = scale.copy() if old is None else np.asarray(old) * scale
        out.normalized = True
        out.flux_unit = (
            data.flux_unit if self.settings["method"] == "model_ratio" else ""
        )
        history = deepcopy(getattr(out, "continuum_history", []))
        history.append(
            dict(settings=deepcopy(self.settings), flags=deepcopy(self.flags))
        )
        out.continuum_history = history
        return out


def from_array(spectrum, continuum):
    """Wrap a supplied/callback continuum without fitting or changing inputs."""
    data = Spec1D.from_spectrum(spectrum)
    array = np.asarray(continuum, dtype=float)
    if array.shape != data.shape:
        raise ValueError("Continuum must match the flux shape")
    valid = np.isfinite(array) & (array > 0) & data.valid
    flags = tuple(
        () if np.any(Spec1D._column(valid, i)[:n]) else ("INVALID_CONTINUUM",)
        for i, n in enumerate(data.numpix)
    )
    return ContinuumResult(
        array.copy(),
        valid,
        flags,
        {"method": "provided"},
        tuple({} for _ in data.numpix),
        None if data.wave is None else data.wave.copy(),
        data.numpix.copy(),
    )


def fit(
    spectrum,
    *,
    method="chebyshev",
    degree=None,
    coordinate="pixel",
    weighted=False,
    clip=None,
    max_iterations=5,
    window=151,
    percentile=90.0,
    bins=20,
    smoothing=None,
    model=None,
    ratio_bounds=(0.1, 10.0),
    mode="mirror",
    on_failure="raise",
):
    """Estimate a continuum independently in each spectral order.

    Parameters
    ----------
    spectrum : Spec1D or Doppler-compatible object
        Supplies flux, wave, optional err/mask, and optional numpix. Not modified.
    method : {'chebyshev', 'polynomial', 'spline', 'median', 'percentile', 'model_ratio'}
        Polynomial is an alias of Chebyshev (same polynomial space). Percentile
        bins positive flux and linearly interpolates bin estimates, with constant
        endpoints. Model_ratio smooths flux/model; both inputs must already use
        comparable normalization, wavelength sampling, frame, and resolution.
    degree : int, optional
        Polynomial degree (default 4) or spline degree (default 3, range 1..5).
    coordinate : {'pixel', 'wavelength'}
        Fitting coordinate, rescaled per order. Pixel is suitable for APOGEE's
        regular log-wavelength grid. Wavelength mode requires distinct samples.
    weighted : bool
        Inverse-variance polynomial/spline least squares; requires errors.
        False by default. Filtering methods do not accept weighted fitting.
    clip : (lower, upper), optional
        Iterative residual rejection in robust MAD sigma units for polynomial or
        spline fits only. Each side can be None. Disabled by default. This is
        independent of FERRE's REJECTCONT; no FERRE setting is changed.
    max_iterations : int
        Maximum rejection/refit rounds after the initial fit.
    window : odd int
        Running filter width in pixels, for median/model_ratio (default 151).
    percentile, bins : float, int
        Percentile (0..100) and number of bins for the percentile estimator.
    smoothing : nonnegative float, optional
        Spline weighted residual sum bound. Default uses supplied errors, or a
        robust adjacent-difference noise estimate; it scales with flux units.
    model : array or spectrum object, optional
        Required for model_ratio; same order layout. A spectrum's masks are used.
    ratio_bounds : (float, float)
        Positive allowed data/model ratio range before filtering.
    mode : str
        scipy.ndimage median-filter boundary mode; default mirror. Filtering
        fills rejected positions by within-order linear interpolation first.
    on_failure : {'raise', 'mask'}
        Raise for insufficient support/nonpositive continuum, or return NaN for
        that order with a quality flag. Never substitute a valid unity continuum.

    Returns
    -------
    ContinuumResult
        Divisor, original-sample fit mask, method parameters, and per-order flags.
        Applying the result is a separate operation with copy-by-default behavior.
    """
    data = Spec1D.from_spectrum(spectrum)
    methods = (
        "chebyshev",
        "polynomial",
        "spline",
        "median",
        "percentile",
        "model_ratio",
    )
    if method not in methods:
        raise ValueError(f"Unknown method {method!r}; choose from {methods}")
    if coordinate not in ("pixel", "wavelength") or on_failure not in ("raise", "mask"):
        raise ValueError("Invalid coordinate or on_failure setting")
    fitted = method in ("chebyshev", "polynomial", "spline")
    if weighted and (not fitted or data.err is None):
        raise ValueError(
            "Weighted fitting requires a polynomial/spline method and uncertainties"
        )
    if degree is None:
        degree = 3 if method == "spline" else 4
    if (
        not isinstance(degree, (int, np.integer))
        or degree < 0
        or (method == "spline" and not 1 <= degree <= 5)
    ):
        raise ValueError("Invalid fitting degree")
    if not isinstance(max_iterations, (int, np.integer)) or max_iterations < 0:
        raise ValueError("max_iterations must be a nonnegative integer")
    if clip is not None:
        if (
            not fitted
            or len(clip) != 2
            or any(v is not None and (not np.isfinite(v) or v <= 0) for v in clip)
        ):
            raise ValueError(
                "clip requires a fit method and positive lower/upper limits or None"
            )
    if method in ("median", "model_ratio"):
        if not isinstance(window, (int, np.integer)) or window < 1 or window % 2 != 1:
            raise ValueError("window must be a positive odd number of pixels")
        if mode not in ("reflect", "constant", "nearest", "mirror", "wrap"):
            raise ValueError("Invalid median-filter boundary mode")
    if method == "percentile" and (
        not np.isfinite(percentile)
        or not 0 <= percentile <= 100
        or not isinstance(bins, (int, np.integer))
        or bins < 2
    ):
        raise ValueError("percentile must be 0..100 and bins an integer >=2")
    if smoothing is not None and (not np.isfinite(smoothing) or smoothing < 0):
        raise ValueError("smoothing must be finite and nonnegative")
    if coordinate == "wavelength" and data.wave is None:
        raise ValueError("Wavelength-coordinate fitting requires wavelengths")
    model_data = None
    model_flux = None
    if method == "model_ratio":
        if model is None:
            raise ValueError("model_ratio requires a reference model")
        if (
            len(ratio_bounds) != 2
            or not np.isfinite(ratio_bounds).all()
            or not 0 < ratio_bounds[0] < ratio_bounds[1]
        ):
            raise ValueError("Invalid ratio bounds")
        if hasattr(model, "flux"):
            model_data = Spec1D.from_spectrum(model)
            if not _aligned(data.wave, model_data.wave) or not np.array_equal(
                data.numpix, model_data.numpix
            ):
                raise ValueError(
                    "Reference model wavelengths/order layout do not match"
                )
            if (
                data.wavevac != model_data.wavevac
                or data.wave_unit != model_data.wave_unit
            ):
                raise ValueError(
                    "Reference model wavelength convention/units do not match"
                )
            if (
                data.flux_unit is not None
                and model_data.flux_unit is not None
                and data.flux_unit != model_data.flux_unit
            ):
                raise ValueError("Reference model flux units do not match")
            model_flux = model_data.flux
        else:
            model_flux = np.asarray(model, dtype=float)
        if model_flux.shape != data.shape:
            raise ValueError("Reference model must match the flux shape")
    continuum = np.full(data.shape, np.nan)
    used = np.zeros(data.shape, bool)
    flags, parameters = [], []
    for order, count in enumerate(data.numpix):
        y = Spec1D._column(data.flux, order)[:count].copy()
        good = Spec1D._column(data.valid, order)[:count].copy()
        errors = None if data.err is None else Spec1D._column(data.err, order)[:count]
        x = (
            np.arange(count, dtype=float)
            if coordinate == "pixel"
            else Spec1D._column(data.wave, order)[:count].copy()
        )
        good &= np.isfinite(x) & np.isfinite(y)
        if method == "model_ratio":
            m = Spec1D._column(model_flux, order)[:count]
            valid_model = np.isfinite(m) & (m > 0)
            if model_data is not None:
                valid_model &= Spec1D._column(model_data.valid, order)[:count]
            ratio = np.full(count, np.nan)
            np.divide(y, m, out=ratio, where=valid_model)
            y = ratio
            good &= valid_model & (y >= ratio_bounds[0]) & (y <= ratio_bounds[1])
        if method == "percentile":
            good &= y > 0
        try:
            minimum = degree + 1 if fitted else 2
            if np.sum(good) < minimum:
                raise ValueError("INSUFFICIENT_PIXELS")
            xmin, xmax = np.min(x[good]), np.max(x[good])
            z = np.zeros(count) if xmax == xmin else 2 * (x - xmin) / (xmax - xmin) - 1
            details = {}
            for iteration in range((max_iterations if clip else 0) + 1):
                if np.sum(good) < minimum:
                    raise ValueError("INSUFFICIENT_PIXELS")
                if method in ("chebyshev", "polynomial"):
                    matrix = np.polynomial.chebyshev.chebvander(z[good], degree)
                    rhs = y[good]
                    if weighted:
                        matrix = matrix / errors[good, None]
                        rhs = rhs / errors[good]
                    coef, _, rank, _ = np.linalg.lstsq(matrix, rhs, rcond=None)
                    if rank < degree + 1:
                        raise ValueError("SINGULAR_FIT")
                    curve = np.polynomial.chebyshev.chebval(z, coef)
                    details = dict(
                        coefficients=coef.tolist(), domain=[float(xmin), float(xmax)]
                    )
                elif method == "spline":
                    from scipy.interpolate import UnivariateSpline

                    select = np.flatnonzero(good)
                    select = select[np.argsort(z[select])]
                    if np.any(np.diff(z[select]) <= 0):
                        raise ValueError("DUPLICATE_COORDINATES")
                    bound = smoothing
                    if bound is None:
                        if weighted:
                            bound = len(select)
                        elif errors is not None:
                            bound = float(np.sum(errors[select] ** 2))
                        else:
                            diff = np.diff(y[select])
                            noise = (
                                1.4826
                                * np.median(np.abs(diff - np.median(diff)))
                                / np.sqrt(2)
                            )
                            bound = len(select) * noise**2
                    spline = UnivariateSpline(
                        z[select],
                        y[select],
                        w=1 / errors[select] if weighted else None,
                        k=degree,
                        s=bound,
                        ext=0,
                    )
                    curve = spline(z)
                    details = dict(
                        knots=spline.get_knots().tolist(),
                        coefficients=spline.get_coeffs().tolist(),
                        domain=[float(xmin), float(xmax)],
                        smoothing=float(bound),
                    )
                elif method == "percentile":
                    edges = np.linspace(xmin, xmax, bins + 1)
                    groups = np.clip(
                        np.searchsorted(edges, x[good], side="right") - 1, 0, bins - 1
                    )
                    centres, levels = [], []
                    for group in np.unique(groups):
                        sel = groups == group
                        centres.append(float(np.mean(x[good][sel])))
                        levels.append(float(np.percentile(y[good][sel], percentile)))
                    if len(centres) < 2:
                        raise ValueError("INSUFFICIENT_BINS")
                    curve = np.interp(x, centres, levels)
                    details = dict(bin_coordinates=centres, bin_levels=levels)
                else:
                    from scipy.ndimage import median_filter

                    # Filtering is in pixel space even when fitting coordinates
                    # are wavelengths; there is never interpolation across orders.
                    pixels = np.arange(count)
                    filled = np.interp(pixels, pixels[good], y[good])
                    curve = median_filter(filled, size=window, mode=mode)
                if clip is None or iteration == max_iterations:
                    break
                residual = y - curve
                centre = np.median(residual[good])
                sigma = max(
                    1.4826 * np.median(np.abs(residual[good] - centre)),
                    32 * np.finfo(float).eps * max(1.0, float(np.max(abs(y[good])))),
                )
                keep = good.copy()
                if clip[0] is not None:
                    keep &= residual - centre >= -clip[0] * sigma
                if clip[1] is not None:
                    keep &= residual - centre <= clip[1] * sigma
                if np.array_equal(keep, good):
                    break
                good = keep
            if not np.all(np.isfinite(curve[good])) or np.any(curve[good] <= 0):
                raise ValueError("INVALID_CONTINUUM")
            # Masked samples cannot become fit data, but a valid polynomial may
            # still be evaluated at their positions. Invalid scales stay NaN.
            curve[~np.isfinite(curve) | (curve <= 0)] = np.nan
            Spec1D._column(continuum, order)[:count] = curve
            Spec1D._column(used, order)[:count] = good
            details["pixels_used"] = int(sum(good))
            flags.append(())
            parameters.append(details)
        except (ValueError, np.linalg.LinAlgError) as error:
            if on_failure == "raise":
                raise ValueError(f"Order {order}: {error}") from error
            flags.append((str(error),))
            parameters.append({"pixels_used": int(sum(good))})
    settings = dict(
        method=method,
        degree=degree,
        coordinate=coordinate,
        weighted=weighted,
        clip=clip,
        max_iterations=max_iterations,
        window=window,
        percentile=percentile,
        bins=bins,
        smoothing=smoothing,
        ratio_bounds=ratio_bounds if method == "model_ratio" else None,
        mode=mode,
        on_failure=on_failure,
    )
    return ContinuumResult(
        continuum,
        used,
        tuple(flags),
        settings,
        tuple(parameters),
        None if data.wave is None else data.wave.copy(),
        data.numpix.copy(),
    )
