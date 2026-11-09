"""Declarative, versioned ASCII JSON recipes for prepared spectral pipelines."""

from copy import deepcopy
from dataclasses import dataclass, field
import json
from pathlib import Path
from importlib.resources import files

import numpy as np

from .ferre import Grid
from .pipeline import Pipeline
from .batch import batch_spectra
from .spectrum import Spec1D
from . import windows


@dataclass
class RecipeResult:
    """Products and audit trail from sequential recipe execution.

    spectrum is the latest prepared/normalized copy; global_fit is the last
    selected fit. history records every attempted stage and candidate fit.
    Errors preserve earlier products. Element failures retain successful peers.
    """

    spectrum: Spec1D
    global_fit: object = None
    element_fits: dict = field(default_factory=dict)
    history: list = field(default_factory=list)
    errors: dict = field(default_factory=dict)
    config: dict = field(default_factory=dict)
    completed: bool = False

    @property
    def parameters(self):
        return None if self.global_fit is None else dict(self.global_fit.parameters)

    @property
    def abundances(self):
        return {
            name: item.result["abundance"] if item.success else None
            for name, item in self.element_fits.items()
        }

    @property
    def success(self):
        return self.completed and not self.errors


class PipelineRecipe:
    """Execute explicit steps from an ASCII JSON file on one prepared Spec1D.

    Use read(path) to load a recipe. Paths resolve relative to that file.
    Supported steps are prepare, normalize, initial, global_fit, and abundances.
    The first three are required before a fit (normalize may explicitly skip).
    Global fits can repeat, compare compatible grids, and try multiple starts.
    Each fit uses the latest normalized spectrum. A normalization after a fit
    requires a subsequent global fit before abundances.

    The fit method accepts Spec1D objects; the CLI uses read_spectrum for FITS/CSV.
    This runner does not calibrate survey results, select
    LSF groups automatically, or provide historical DR17 quality logic. The
    DR17-style example describes a modern approximation of the fitting stages.
    """

    def __init__(self, config, *, base_dir=".", executable=None):
        self.config = deepcopy(config)
        self.base_dir = Path(base_dir).resolve()
        self.executable = executable
        self._validate()

    @classmethod
    def read(cls, filename, *, executable=None):
        """Read a version-1 JSON recipe; reject unknown stages/options early."""
        path = Path(filename).resolve()
        return cls(
            json.loads(path.read_text(encoding="ascii")),
            base_dir=path.parent,
            executable=executable,
        )

    @classmethod
    def example(cls, name="basic", *, executable=None):
        """Load the installed 'basic' or 'dr17-style' recipe without a checkout.

        The DR17-style template must be copied and edited to supply real grids.
        Files reside in aspcap/data/recipes and may be copied for customization.
        """
        if name not in ("basic", "dr17-style"):
            raise ValueError("Unknown example; choose basic or dr17-style")
        return cls.read(
            files("aspcap").joinpath("data/recipes", "pipeline-" + name + ".json"),
            executable=executable,
        )

    def _validate(self):
        c = self.config
        if not isinstance(c, dict) or c.get("version") != 1:
            raise ValueError("Recipe must be a mapping with version=1")
        if set(c) - {"version", "description", "grids", "steps", "ferre"}:
            raise ValueError("Unknown recipe-level options")
        if not isinstance(c.get("grids"), dict) or not c["grids"]:
            raise ValueError("Provide named grids")
        if any(not isinstance(v, str) for v in c["grids"].values()):
            raise ValueError("Each grid must be a header path")
        if not isinstance(c.get("steps"), list) or not c["steps"]:
            raise ValueError("Provide an ordered steps list")
        allowed = {
            "prepare": {"grid"},
            "normalize": {"method", "options", "value"},
            "initial": {"parameters"},
            "global_fit": {
                "grids",
                "starts",
                "free_labels",
                "weights",
                "ties",
                "tie_mode",
            },
            "abundances": {"elements"},
        }
        prepared = normalized = initialized = fitted = False
        for step in c["steps"]:
            if not isinstance(step, dict) or step.get("type") not in allowed:
                raise ValueError("Unknown or missing step type")
            kind = step["type"]
            if set(step) - allowed[kind] - {"type", "name", "description"}:
                raise ValueError(f"Unknown options for {kind}")
            if kind == "prepare":
                if prepared or step.get("grid") not in c["grids"]:
                    raise ValueError("Prepare once with a named grid")
                prepared = True
            elif not prepared:
                raise ValueError("prepare must be first")
            elif kind == "normalize":
                if step.get("method") not in (
                    "skip",
                    "constant",
                    "chebyshev",
                    "polynomial",
                    "spline",
                    "median",
                    "percentile",
                    "model_ratio",
                ):
                    raise ValueError("Unknown normalization method")
                if step["method"] == "model_ratio" and not fitted:
                    raise ValueError("model_ratio requires a preceding global fit")
                normalized, fitted = True, False
            elif kind == "initial":
                if not isinstance(step.get("parameters"), dict):
                    raise ValueError("initial requires parameters")
                initialized = True
                fitted = False
            elif kind == "global_fit":
                if not (normalized and initialized):
                    raise ValueError(
                        "global_fit requires normalization and initial steps"
                    )
                if any(name not in c["grids"] for name in step.get("grids", [])):
                    raise ValueError("Unknown candidate grid")
                fitted = True
            elif kind == "abundances":
                if (
                    not fitted
                    or not isinstance(step.get("elements"), list)
                    or not step["elements"]
                ):
                    raise ValueError(
                        "abundances requires a current fit and an element list"
                    )
        if not fitted:
            raise ValueError(
                "Recipe must end with a current global fit or abundance stage"
            )
        if not any(s["type"] == "global_fit" for s in c["steps"]):
            raise ValueError("Recipe requires a global_fit step")

    def _path(self, value):
        return self.base_dir / value

    def _grid(self, name):
        executable = self.executable or self.config.get("ferre", "ferre")
        if "/" in str(executable) and not Path(executable).is_absolute():
            executable = self._path(executable)
        return Grid(self._path(self.config["grids"][name]), executable=executable)

    def _weights(self, specification, grid):
        if specification is None:
            return None
        if isinstance(specification, list):
            return specification
        if not isinstance(specification, dict):
            raise ValueError("Weights must be an array or a source mapping")
        if set(specification) == {"dr17"}:
            return windows.load_dr17(
                specification["dr17"], wavelength=grid.wavelengths
            )[1]
        if set(specification) == {"file"}:
            return windows.load_weights(
                self._path(specification["file"]), wavelength=grid.wavelengths
            )[1]
        if set(specification) == {"legacy_file"}:
            return windows.load_weights(
                self._path(specification["legacy_file"]), wavelength=grid.wavelengths
            )[1]
        raise ValueError("Unknown weight source; use file, legacy_file, or dr17")

    @batch_spectra
    def fit(
        self,
        spectrum,
        *,
        initial=None,
        star_id="star",
        star_ids=None,
        raise_on_error=False,
        stop_after=None,
    ):
        """Run the ordered recipe without mutating the input spectrum.

        initial optionally overrides the recipe's supplied starting parameters.
        A failed stage stops execution, recording the error and prior results.
        Failed candidate starts are recorded but another successful start can
        still win. Grid comparisons use unpenalized reduced chi-square and
        require identical wavelength sampling, labels, and wavelength frames.
        """
        if stop_after is not None and (
            isinstance(stop_after, bool)
            or not isinstance(stop_after, int)
            or not 1 <= stop_after <= len(self.config["steps"])
        ):
            raise ValueError("stop_after must be a one-based recipe step number")
        if not isinstance(spectrum, Spec1D):
            raise TypeError("Recipe input must be Spec1D")
        result = RecipeResult(spectrum.copy(), config=deepcopy(self.config))
        grid = parameters = global_weights = None
        for number, step in enumerate(self.config["steps"], 1):
            record = {
                "step": number,
                "type": step["type"],
                "name": step.get("name", step["type"]),
            }
            if stop_after is not None and number > stop_after:
                break
            result.history.append(record)
            try:
                kind = step["type"]
                if kind == "prepare":
                    grid = self._grid(step["grid"])
                    result.spectrum = Pipeline(grid).prepare(result.spectrum)
                elif kind == "normalize":
                    method = step["method"]
                    options = deepcopy(step.get("options", {}))
                    if method == "skip":
                        setting = None
                    elif method == "constant":
                        setting = np.full(
                            result.spectrum.flux.shape, float(step["value"])
                        )
                    else:
                        setting = dict(options, method=method)
                        if method == "model_ratio":
                            if result.spectrum.flux.ndim != 1:
                                raise ValueError(
                                    "model_ratio recipe requires a flat native-pixel spectrum"
                                )
                            model = grid(parameters)
                            # Preserve multi-order pixel layout for continuum fitting.
                            layout = result.spectrum.flatten()
                            if not np.allclose(
                                layout.wave, model.wave, rtol=0, atol=1e-6
                            ):
                                raise ValueError(
                                    "model_ratio recipe requires a flat native-pixel spectrum"
                                )
                            setting["model"] = model.flux.reshape(
                                result.spectrum.flux.shape
                            )
                    result.spectrum, continuum = Pipeline(grid).normalize(
                        result.spectrum, continuum=setting
                    )
                    record["flags"] = None if continuum is None else continuum.flags
                elif kind == "initial":
                    parameters = dict(
                        step["parameters"] if initial is None else initial
                    )
                    record["parameters"] = parameters.copy()
                elif kind == "global_fit":
                    candidates = step.get("grids") or [
                        next(
                            name
                            for name in self.config["grids"]
                            if self._path(self.config["grids"][name]).resolve()
                            == grid.path
                        )
                    ]
                    best = None
                    record["attempts"] = []
                    for name in candidates:
                        candidate = self._grid(name)
                        if (
                            candidate.labels != grid.labels
                            or candidate.vacuum != grid.vacuum
                            or not np.array_equal(
                                candidate.wavelengths, grid.wavelengths
                            )
                        ):
                            raise ValueError(
                                "Candidate grids must share labels, native wavelengths, and frame"
                            )
                        weights = self._weights(step.get("weights"), candidate)
                        starts = step.get("starts") or [parameters]
                        for start in starts:
                            attempt = {"grid": name, "start": deepcopy(start)}
                            record["attempts"].append(attempt)
                            try:
                                fit = candidate.fit(
                                    result.spectrum,
                                    start=start,
                                    free_labels=step.get("free_labels"),
                                    weights=weights,
                                    ties=step.get("ties"),
                                    tie_mode=step.get("tie_mode", "relative"),
                                )
                                attempt["reduced_chi2"] = fit.reduced_chi2
                                if (
                                    best is None
                                    or fit.reduced_chi2 < best[0].reduced_chi2
                                ):
                                    best = fit, candidate, weights, name
                            except (ValueError, RuntimeError, OSError) as error:
                                attempt["error"] = str(error)
                    if best is None:
                        raise ValueError(
                            "All grid/start candidates failed; see history"
                        )
                    result.global_fit, grid, global_weights, record["selected_grid"] = (
                        best
                    )
                    parameters = dict(result.global_fit.parameters)
                    result.element_fits = {}
                    record["parameters"] = parameters.copy()
                elif kind == "abundances":
                    elements = []
                    for original in step["elements"]:
                        item = deepcopy(original)
                        target = grid
                        if "grid" in item:
                            item["grid"] = str(self._path(item["grid"]))
                            target = Grid(item["grid"], executable=grid.executable)
                        if "weights" in item:
                            item["weights"] = self._weights(item["weights"], target)
                        if "windows_file" in item:
                            item["windows"] = windows.load_windows(
                                self._path(item.pop("windows_file"))
                            )
                        elements.append(item)
                    pipeline = Pipeline(grid, elements)
                    fits = pipeline.fit_elements(
                        result.spectrum,
                        parameters,
                        star_id=star_id,
                        weights=global_weights,
                        raise_on_error=raise_on_error,
                    )
                    result.element_fits.update(fits)
                    for name, fit in fits.items():
                        if not fit.success:
                            result.errors[f"{number}:{name}"] = fit.error
                record["status"] = "complete"
            except (ValueError, RuntimeError, OSError, KeyError, TypeError) as error:
                if raise_on_error:
                    raise
                record["status"] = "failed"
                result.errors[str(number)] = str(error)
                break
        result.completed = not any(
            record.get("status") == "failed" for record in result.history
        )
        return result
