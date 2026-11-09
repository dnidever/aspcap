"""Single-spectrum orchestration of existing continuum and FERRE interfaces."""

from collections.abc import Mapping
from copy import deepcopy
from dataclasses import dataclass, field, replace
from math import isfinite
from pathlib import Path

import numpy as np

from . import continuum as continuum_module
from .continuum import ContinuumResult
from .ferre import Grid
from .spectrum import Spec1D


class _ConfiguredDefault:
    def __repr__(self):
        return "<configured>"


_UNSET = _ConfiguredDefault()


@dataclass(frozen=True)
class ElementConfig:
    """An explicit element-fitting prescription.

    species, label : str
        Result name and exact free grid label. Other labels stay fixed unless explicitly tied.
    windows : sequence of (float, float), optional
        Inclusive intervals restricting element weights. Supply windows or weights.
    weights : sequence of float, optional
        Fractional native-pixel weights for this element. Bad pixels and zero
        global weights are excluded; positive global magnitudes are not applied.
    ties : mapping, optional
        Dependent-label affine relations as in Grid.fit.
    tie_mode : {"relative", "absolute"}, optional
        Relative to initial parameters by default.
    add_labels : sequence of str, optional
        Fixed logarithmic coordinates added to report [X/H].
    grid : Grid or path-like, optional
        Supplied element/mini-grid; defaults to the pipeline's global grid.
        It must use matching native wavelengths and compatible flux/LSF scales.
    initial : float, optional
        Starting value for the free label. Required when that label is absent
        from the global solution, for example a dedicated mini-grid abundance.
    """

    species: str
    label: str
    windows: object = None
    add_labels: tuple = ()
    grid: object = field(default=None, repr=False)
    initial: object = None
    weights: object = field(default=None, repr=False)
    ties: object = field(default=None, repr=False)
    tie_mode: str = "relative"

    def __post_init__(self):
        if not isinstance(self.species, str) or not self.species.strip():
            raise ValueError("Element species must be a nonempty name")
        if not isinstance(self.label, str) or not self.label.strip():
            raise ValueError("Element label must be a nonempty name")
        windows = (
            None
            if self.windows is None
            else tuple(tuple(float(v) for v in pair) for pair in self.windows)
        )
        if windows is not None and (
            not windows
            or any(
                len(pair) != 2 or not all(map(isfinite, pair)) or pair[0] > pair[1]
                for pair in windows
            )
        ):
            raise ValueError("Element windows must contain finite (lower, upper) pairs")
        if windows is None and self.weights is None:
            raise ValueError("Supply element windows or weights")
        if self.weights is not None:
            weights = np.asarray(self.weights, dtype=float)
            if (
                weights.ndim != 1
                or not weights.size
                or not np.all(np.isfinite(weights) & (weights >= 0))
            ):
                raise ValueError("Element weights must be a finite nonnegative vector")
            object.__setattr__(self, "weights", tuple(weights))
        if self.tie_mode not in ("relative", "absolute"):
            raise ValueError("Invalid element tie_mode")
        object.__setattr__(self, "ties", deepcopy(self.ties))
        labels = tuple(self.add_labels)
        if len(set(labels)) != len(labels) or self.label in labels:
            raise ValueError("add_labels must be distinct fixed labels")
        if self.initial is not None:
            if not isfinite(float(self.initial)):
                raise ValueError("Element initial value must be finite")
            object.__setattr__(self, "initial", float(self.initial))
        object.__setattr__(self, "windows", windows)
        object.__setattr__(self, "add_labels", labels)


@dataclass
class ElementFitResult:
    """A successful measurement or an independently recorded element failure.

    species : str
        Configured measurement name.
    result : dict or None
        Grid.fit_abundance report, including abundance, uncertainty, and flags.
    pixels : tuple
        Reusable native diagnostic rows: wavelength, normalized flux/error,
        model, and window weight. Empty when fitting failed.
    error : str or None
        Failure message, including retained FERRE log locations when available.
    """

    species: str
    result: object = None
    pixels: tuple = field(default=(), repr=False)
    error: object = None

    @property
    def success(self):
        """Whether this element produced a measurement without an execution error."""
        return self.result is not None and self.error is None

    @property
    def flags(self):
        """Measurement flags, or ELEMENT_FIT_FAILED on a recorded failure."""
        if self.error is not None:
            return ("ELEMENT_FIT_FAILED",)
        return tuple(filter(None, (self.result or {}).get("flags", "").split("|")))


@dataclass
class PipelineResult:
    """Results and intermediate products from one pipeline run.

    original_spectrum, prepared_spectrum, spectrum : Spec1D or None
        Independent input snapshot, validated input, and normalized spectrum.
        Later fields remain None when a preceding stage fails.
    continuum : ContinuumResult or None
        The estimate applied in this run. None means normalization was skipped.
        No additional continuum correction is made during abundance fitting.
    initial_parameters : dict or None
        Complete starting parameters supplied by the caller.
    global_fit : FerreFit or None
        Global parameter fit, uncertainty estimates, statistic, and native model.
    element_fits : dict of str to ElementFitResult
        Requested measurements, including individual failures. Pixel iterators
        from Grid.fit_abundance are materialized into reusable tuples.
    flags, errors, config : tuple, dict, dict
        Processing flags, failed-stage messages, and a configuration snapshot.
        These flags do not constitute a complete scientific quality assessment.
    """

    original_spectrum: Spec1D
    prepared_spectrum: object = None
    spectrum: object = None
    continuum: object = None
    initial_parameters: object = None
    global_fit: object = None
    element_fits: dict = field(default_factory=dict)
    flags: tuple = ()
    errors: dict = field(default_factory=dict)
    config: dict = field(default_factory=dict)

    def __repr__(self):
        """Summarize the run without printing full spectra or pixel tables."""
        return (
            f"PipelineResult(star_id={self.config.get('star_id')!r}, "
            f"success={self.success}, parameters={self.parameters!r}, "
            f"abundances={self.abundances!r}, flags={self.flags!r})"
        )

    @property
    def parameters(self):
        """A copy of the global fitted labels, or None after a global failure."""
        return None if self.global_fit is None else dict(self.global_fit.parameters)

    @property
    def abundances(self):
        """Map each attempted species to [X/H], or None for a failed element."""
        return {
            name: item.result["abundance"] if item.success else None
            for name, item in self.element_fits.items()
        }

    @property
    def residuals(self):
        """Observed-minus-model residuals in native order; excluded pixels are NaN."""
        if self.global_fit is None or self.spectrum is None:
            return None
        native = (
            self.spectrum.to_native()
            if self.spectrum._native_indices is not None
            else self.spectrum.flatten()
        )
        valid = native.valid.copy()
        weights = self.config.get("weights")
        if weights is not None:
            valid &= np.asarray(weights) > 0
        return np.where(valid, native.flux - np.asarray(self.global_fit.model), np.nan)

    @property
    def success(self):
        """Whether the global stage and all attempted element fits completed."""
        return (
            self.global_fit is not None
            and not self.errors
            and all(item.success for item in self.element_fits.values())
        )


class Pipeline:
    """Run normalization, global fitting, and configured element fits.

    Parameters
    ----------
    grid : Grid or path-like
        Supplied global library. The default free labels are all its labels.
    elements : sequence of ElementConfig/dict or mapping, optional
        Element prescriptions. A mapping can use species names as keys and
        dictionaries of label/windows/add_labels as values. Empty means only
        the global stage. Names must be unique.
    continuum : mapping, array-like, ContinuumResult, or None, optional
        Mapping of options for continuum.fit, or supplied divisors matching
        the input order layout. None explicitly means the input is already on
        the model's flux scale and normalization is skipped. Can be overridden
        per spectrum. No implicit continuum estimator is chosen.
    initial : mapping or sequence, optional
        Explicit starting parameters; can be overridden per spectrum. No
        automatic estimator or midpoint fallback is used by the pipeline.
    free_labels : sequence of str, optional
        Global labels to vary. Other labels stay at their initial values.
    weights : array-like, optional
        Global weights in reconstructed native pixel order. Zero weights also
        exclude pixels from element fits; positive magnitudes affect only the
        global fit, since element fits use their own windows and weights.
    ties : mapping, optional
        Global-fit affine relations as in Grid.fit; element ties are separate.
    tie_mode : {"relative", "absolute"}, optional
        Global-fit relation mode, relative by default.
    keep_runs : bool, optional
        Retain successful FERRE run files for both fitting stages.
    run_root : path-like, optional
        Existing parent directory for FERRE run files.

    Notes
    -----
    The input must be a prepared rest-frame Spec1D with uncertainties and native
    wavelengths. Preparation validates and maps existing pixels; it does not
    read FITS, shift velocities, resample arbitrary observations, convolve an
    LSF, select libraries, or enforce microturbulence parameter relations.
    """

    def __init__(
        self,
        grid,
        elements=(),
        *,
        continuum=None,
        initial=None,
        free_labels=None,
        weights=None,
        ties=None,
        tie_mode="relative",
        keep_runs=False,
        run_root=None,
    ):
        self.grid = grid if isinstance(grid, Grid) else Grid(grid)
        self.ties = deepcopy(ties)
        self.tie_mode = tie_mode
        if isinstance(elements, Mapping):
            if "species" in elements and "label" in elements:
                elements = [elements]
            else:
                converted = []
                for name, options in elements.items():
                    options = dict(options)
                    if options.get("species", name) != name:
                        raise ValueError("Element mapping key must match species")
                    options["species"] = name
                    converted.append(options)
                elements = converted
        configurations = tuple(
            item if isinstance(item, ElementConfig) else ElementConfig(**item)
            for item in elements
        )
        self.elements = tuple(
            (
                replace(item, grid=Grid(item.grid, executable=self.grid.executable))
                if item.grid is not None and not isinstance(item.grid, Grid)
                else item
            )
            for item in configurations
        )
        if len({item.species for item in self.elements}) != len(self.elements):
            raise ValueError("Element species names must be unique")
        for item in self.elements:
            target = self.grid if item.grid is None else item.grid
            if item.label not in target.labels or any(
                label not in target.labels for label in item.add_labels
            ):
                raise ValueError(f"Unknown grid labels in element {item.species}")
        self.free_labels = (
            tuple(label for label in self.grid.labels if label not in (ties or {}))
            if free_labels is None
            else (free_labels,) if isinstance(free_labels, str) else tuple(free_labels)
        )
        if (
            not self.free_labels
            or len(set(self.free_labels)) != len(self.free_labels)
            or any(label not in self.grid.labels for label in self.free_labels)
        ):
            raise ValueError("free_labels must be distinct labels in the global grid")
        self.continuum = deepcopy(continuum)
        self.initial = deepcopy(initial)
        self.weights = self._weights(weights)
        self.keep_runs, self.run_root = bool(keep_runs), run_root

    def _weights(self, weights):
        if weights is None:
            return None
        array = np.asarray(weights, dtype=float)
        if array.shape != (self.grid.npix,) or not np.all(
            np.isfinite(array) & (array >= 0)
        ):
            raise ValueError(
                "Pipeline weights must be finite, nonnegative native pixels"
            )
        return array.copy()

    @staticmethod
    def _native(spectrum):
        return (
            spectrum.to_native()
            if spectrum._native_indices is not None
            else spectrum.flatten()
        )

    def prepare(self, spectrum):
        """Return an independent, validated Spec1D with its order layout retained.

        Requires native rest-frame wavelengths, a common matching air/vacuum
        frame, Angstrom units, and errors. Nonfinite pixels/errors remain
        masked through Grid's shared validation. No flux processing is done.
        """
        if not isinstance(spectrum, Spec1D):
            raise TypeError(
                "Pipeline input must be Spec1D; convert other objects explicitly"
            )
        prepared = spectrum.copy()
        self.grid._spectrum_fit_inputs(prepared)
        return prepared

    def normalize(self, spectrum, *, continuum=_UNSET):
        """Return (normalized copy, ContinuumResult), or (copy, None) when skipped.

        Estimation is order-by-order. Supplied divisors must match the original
        order layout. Applying a result divides both flux and errors. The
        caller's spectrum is never modified.
        """
        setting = self.continuum if continuum is _UNSET else continuum
        if setting is None:
            return spectrum.copy(), None
        if isinstance(setting, ContinuumResult):
            estimate = deepcopy(setting)
        elif isinstance(setting, Mapping):
            estimate = continuum_module.fit(spectrum, **dict(setting))
        else:
            if (
                isinstance(setting, (list, tuple))
                and len(setting)
                and np.ndim(setting[0]) > 0
            ):
                if len(setting) != spectrum.norder:
                    raise ValueError("Provide one continuum array per order")
                array = np.full(spectrum.shape, np.nan)
                for i, (values, count) in enumerate(zip(setting, spectrum.numpix)):
                    values = np.asarray(values, dtype=float)
                    if values.shape != (count,):
                        raise ValueError(
                            "Order continuum length must match its pixel extent"
                        )
                    Spec1D._column(array, i)[:count] = values
                setting = array
            estimate = continuum_module.from_array(spectrum, setting)
        return estimate.apply(spectrum), estimate

    def estimate_initial(self, spectrum, *, initial=_UNSET):
        """Validate supplied starts and return a complete label dictionary.

        Spectrum is accepted for a future estimator interface; this version
        does not infer labels from it. Missing estimates are an explicit error.
        """
        values = self.initial if initial is _UNSET else initial
        if values is None:
            raise ValueError(
                "Supply initial parameters; automatic estimates are not implemented"
            )
        point = self.grid._parameters(values)
        if not self.grid.contains(point):
            raise ValueError("Initial parameters are outside grid coverage")
        return dict(zip(self.grid.labels, point))

    def fit_global(self, spectrum, initial, *, weights=_UNSET):
        """Run Grid.fit on a prepared spectrum using the configured free labels."""
        weights = self.weights if weights is _UNSET else self._weights(weights)
        return self.grid.fit(
            spectrum,
            start=initial,
            free_labels=self.free_labels,
            ties=self.ties,
            tie_mode=self.tie_mode,
            weights=weights,
            keep_runs=self.keep_runs,
            run_root=self.run_root,
        )

    def fit_elements(
        self,
        spectrum,
        parameters,
        *,
        star_id="star",
        weights=_UNSET,
        raise_on_error=False,
        species=None,
    ):
        """Fit configured elements independently with the continuum held fixed.

        The spectrum is already normalized, so each Grid.fit_abundance call
        receives a continuum of ones. This avoids double normalization. Extra
        zero global weights are excluded; positive global weights do not become
        element weight curves. One failure is recorded without discarding others.
        Supplied mini-grids require matching native pixels and explicit initial
        values for free labels absent from the global solution. Optional species
        selects one configured name or a sequence of names for a partial rerun.
        """
        selected = self.elements
        if species is not None:
            names = (species,) if isinstance(species, str) else tuple(species)
            if set(names) - {item.species for item in self.elements}:
                raise ValueError("Unknown species requested for element fitting")
            selected = tuple(item for item in self.elements if item.species in names)
        weights = self.weights if weights is _UNSET else self._weights(weights)
        native = self._native(spectrum)
        if native.err is None:
            raise ValueError("Element fitting requires uncertainties")
        mask = ~native.valid
        if weights is not None:
            mask = mask | (weights <= 0)
        fits = {}
        for item in selected:
            try:
                target = self.grid if item.grid is None else item.grid
                target._spectrum_fit_inputs(native)
                start = {}
                for label in target.labels:
                    if label == item.label and item.initial is not None:
                        start[label] = float(item.initial)
                    elif label in parameters:
                        start[label] = parameters[label]
                    elif label == item.label:
                        raise ValueError(
                            f"Supply an initial value for mini-grid label {label!r}"
                        )
                    else:
                        raise ValueError(
                            f"Global solution has no fixed mini-grid label {label!r}"
                        )
                row, pixels = target.fit_abundance(
                    native.wave,
                    native.flux,
                    native.err,
                    np.ones(native.size),
                    start,
                    star_id=str(star_id),
                    species=item.species,
                    label=item.label,
                    windows=item.windows,
                    weights=item.weights,
                    ties=item.ties,
                    tie_mode=item.tie_mode,
                    add_labels=item.add_labels,
                    mask=mask,
                    keep_runs=self.keep_runs,
                    run_root=self.run_root,
                )
                fits[item.species] = ElementFitResult(item.species, row, tuple(pixels))
            except (ValueError, RuntimeError, OSError) as error:
                if raise_on_error:
                    raise
                fits[item.species] = ElementFitResult(item.species, error=str(error))
        return fits

    def fit(
        self,
        spectrum,
        *,
        initial=_UNSET,
        continuum=_UNSET,
        weights=_UNSET,
        star_id=None,
        raise_on_error=False,
    ):
        """Run the stages and return a PipelineResult, retaining partial results.

        Per-spectrum initial/continuum/weights override constructor settings.
        Explicit continuum=None skips normalization; it does not select a default.
        A global failure stops abundance fitting. An element failure preserves
        other measurements. Expected data/FERRE errors are recorded by default;
        raise_on_error=True re-raises them for debugging. Invalid input types
        raise immediately. No automatic initial estimator or FITS reader is used.
        """
        if not isinstance(spectrum, Spec1D):
            raise TypeError("Pipeline input must be Spec1D")
        result = PipelineResult(spectrum.copy())
        stage = "preparation"
        try:
            if star_id is None:
                star_id = spectrum.metadata.get("star_id") or (
                    Path(spectrum.filename).stem if spectrum.filename else "star"
                )
            run_weights = self.weights if weights is _UNSET else self._weights(weights)
            setting = self.continuum if continuum is _UNSET else continuum
            result.config = {
                "star_id": str(star_id),
                "grid_paths": tuple(str(p) for p in self.grid.paths),
                "executable": str(self.grid.executable),
                "free_labels": self.free_labels,
                "continuum": (
                    deepcopy(dict(setting))
                    if isinstance(setting, Mapping)
                    else "skip" if setting is None else "provided"
                ),
                "weights": None if run_weights is None else run_weights.copy(),
                "ties": deepcopy(self.ties),
                "tie_mode": self.tie_mode,
                "elements": tuple(
                    dict(
                        species=e.species,
                        label=e.label,
                        windows=e.windows,
                        weights=e.weights,
                        ties=deepcopy(e.ties),
                        tie_mode=e.tie_mode,
                        add_labels=e.add_labels,
                        initial=e.initial,
                        grid_paths=tuple(
                            str(p)
                            for p in (self.grid if e.grid is None else e.grid).paths
                        ),
                    )
                    for e in self.elements
                ),
                "keep_runs": self.keep_runs,
                "run_root": None if self.run_root is None else str(self.run_root),
            }
            result.prepared_spectrum = self.prepare(spectrum)
            stage = "normalization"
            result.spectrum, result.continuum = self.normalize(
                result.prepared_spectrum, continuum=setting
            )
            if result.continuum is not None and any(result.continuum.flags):
                result.flags += ("CONTINUUM_FLAGS",)
            stage = "initial_estimates"
            result.initial_parameters = self.estimate_initial(
                result.spectrum, initial=initial
            )
            result.config["initial"] = dict(result.initial_parameters)
            stage = "global_fit"
            result.global_fit = self.fit_global(
                result.spectrum, result.initial_parameters, weights=run_weights
            )
            stage = "element_fits"
            result.element_fits = self.fit_elements(
                result.spectrum,
                result.parameters,
                star_id=star_id,
                weights=run_weights,
                raise_on_error=raise_on_error,
            )
            failed = [
                name for name, item in result.element_fits.items() if not item.success
            ]
            if failed:
                result.flags += ("PARTIAL_ABUNDANCES",)
                result.errors.update(
                    {
                        f"element:{name}": result.element_fits[name].error
                        for name in failed
                    }
                )
        except (ValueError, RuntimeError, OSError) as error:
            if raise_on_error:
                raise
            result.flags += (stage.upper() + "_FAILED",)
            result.errors[stage] = str(error)
        return result
