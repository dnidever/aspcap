"""Real FERRE affine ties, model registration, and validation."""

import os
from pathlib import Path
import shutil
import tempfile
import unittest
import numpy as np
from aspcap import Grid, Pipeline, ElementConfig

ROOT = Path(__file__).resolve().parents[1]


class TieTests(unittest.TestCase):
    def setUp(self):
        exe = os.environ.get("ASPCAP_TEST_FERRE") or shutil.which("ferre")
        if not exe:
            self.skipTest("FERRE unavailable")
        self.grid = Grid(ROOT / "examples/toy-grid.dat", executable=exe)

    def test_relative_ties_and_final_model(self):
        truth = {"ALPHA": 0.2, "METALS": -0.2}
        spec = self.grid(truth)
        spec.err = np.full(8, 0.01)
        with tempfile.TemporaryDirectory() as root:
            fit = self.grid.fit(
                spec,
                start={"ALPHA": 0, "METALS": 0},
                free_labels="ALPHA",
                ties={"METALS": {"ALPHA": -1}},
                keep_runs=True,
                run_root=root,
            )
            self.assertAlmostEqual(fit.parameters["ALPHA"], 0.2, delta=0.003)
            self.assertAlmostEqual(fit.parameters["METALS"], -0.2, delta=0.003)
            np.testing.assert_allclose(
                fit.model, self.grid(fit.parameters).flux, atol=1e-6
            )
            self.assertIsNone(fit.errors["METALS"])
            self.assertTrue(
                any(
                    "TYPETIE=1" in p.read_text() for p in Path(root).glob("*/input.nml")
                )
            )

    def test_absolute_ties_and_pipeline(self):
        spec = self.grid({"ALPHA": 0.2, "METALS": -0.25})
        spec.err = np.full(8, 0.01)
        tie = {"METALS": {"coefficients": {"ALPHA": -0.5}, "offset": -0.15}}
        fit = self.grid.fit(
            spec, start={"ALPHA": 0, "METALS": 0}, ties=tie, tie_mode="absolute"
        )
        self.assertAlmostEqual(fit.parameters["ALPHA"], 0.2, delta=0.003)
        self.assertAlmostEqual(fit.parameters["METALS"], -0.25, delta=0.003)
        pipeline = Pipeline(
            self.grid,
            [
                ElementConfig(
                    "Mg",
                    "ALPHA",
                    weights=np.ones(8),
                    ties=tie,
                    tie_mode="absolute",
                    add_labels=("METALS",),
                )
            ],
            initial={"ALPHA": 0, "METALS": 0},
            ties=tie,
            tie_mode="absolute",
        )
        result = pipeline.fit(spec)
        self.assertTrue(result.success, result.errors)
        self.assertAlmostEqual(result.abundances["Mg"], -0.05, delta=0.004)

    def test_invalid_ties(self):
        spec = self.grid({"ALPHA": 0.2, "METALS": -0.2})
        spec.err = np.full(8, 0.01)
        for ties in (
            {"BOGUS": {"ALPHA": 1}},
            {"ALPHA": {"ALPHA": 1}},
            {"METALS": {"ALPHA": float("nan")}},
            {"METALS": {"ALPHA": 1}, "ALPHA": {"METALS": 1}},
        ):
            with self.assertRaises(ValueError):
                self.grid.fit(spec, free_labels="ALPHA", ties=ties)
