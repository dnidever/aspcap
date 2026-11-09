"""Spectral data, external continuum normalization, and FERRE abundance fits."""

from .abundances import fit_abundance
from .ferre import Grid
from .spectrum import Spec1D
from . import continuum, windows, grids
from .pipeline import Pipeline, PipelineResult, ElementConfig, ElementFitResult

__all__ = [
    "Grid",
    "grids",
    "BatchItem",
    "read_apstar",
    "read_spectrum",
    "PipelineRecipe",
    "RecipeResult",
    "fit_abundance",
    "Spec1D",
    "continuum",
    "windows",
    "Pipeline",
    "PipelineResult",
    "ElementConfig",
    "ElementFitResult",
]

from .recipe import PipelineRecipe, RecipeResult

from .batch import BatchItem
from .readers import read_apstar, read_spectrum
