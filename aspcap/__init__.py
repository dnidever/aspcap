"""Spectral data, external continuum normalization, and FERRE abundance fits."""

from .abundances import fit_abundance
from .ferre import Grid
from .spectrum import Spec1D
from . import continuum, windows
from .pipeline import Pipeline, PipelineResult, ElementConfig, ElementFitResult

__all__ = [
    "Grid",
    "fit_abundance",
    "Spec1D",
    "continuum",
    "windows",
    "Pipeline",
    "PipelineResult",
    "ElementConfig",
    "ElementFitResult",
]
