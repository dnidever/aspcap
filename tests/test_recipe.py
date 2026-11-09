"""Declarative stages execute real fits and retain candidate/failure histories."""

import csv
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

import numpy as np
from aspcap import PipelineRecipe, Grid, Spec1D

ROOT = Path(__file__).resolve().parents[1]
EXAMPLES = ROOT / "examples"


class RecipeTests(unittest.TestCase):
    def setUp(self):
        self.exe = os.environ.get("ASPCAP_TEST_FERRE") or shutil.which("ferre")
        self.config = json.loads((EXAMPLES / "pipeline-basic.json").read_text())

    def spectrum(self):
        if not self.exe:
            self.skipTest("FERRE unavailable")
        grid = Grid(EXAMPLES / "toy-grid.dat", executable=self.exe)
        spec = grid({"ALPHA": 0.2, "METALS": -0.25})
        spec.err = np.full(8, 0.01)
        return spec

    def test_basic_recipe_relative_paths_and_no_input_mutation(self):
        spec = self.spectrum()
        flux = spec.flux.copy()
        result = PipelineRecipe.read(
            EXAMPLES / "pipeline-basic.json", executable=self.exe
        ).fit(spec)
        self.assertTrue(result.success, result.errors)
        self.assertAlmostEqual(result.abundances["Mg"], -0.05, delta=0.003)
        np.testing.assert_array_equal(spec.flux, flux)
        self.assertEqual(len(result.history), 5)

    def test_repeated_fit_multistart_and_normalization(self):
        self.config["steps"][3]["starts"] = [
            dict(ALPHA=99, METALS=0),
            dict(ALPHA=0, METALS=-0.25),
        ]
        self.config["steps"].insert(
            4, dict(type="normalize", method="model_ratio", options={"window": 3})
        )
        self.config["steps"].insert(5, dict(type="global_fit"))
        result = PipelineRecipe(
            self.config, base_dir=EXAMPLES, executable=self.exe
        ).fit(self.spectrum())
        self.assertTrue(result.success, result.errors)
        self.assertIn("error", result.history[3]["attempts"][0])
        self.assertAlmostEqual(result.parameters["ALPHA"], 0.2, delta=0.004)

    def test_invalid_steps_order_and_unknown_options(self):
        for step in [
            dict(type="shell", command="anything"),
            dict(type="global_fit"),
            dict(type="normalize", method="misspelled"),
        ]:
            config = dict(self.config, steps=[step])
            with self.assertRaises(ValueError):
                PipelineRecipe(config)
        self.config["steps"][0]["typo"] = True
        with self.assertRaises(ValueError):
            PipelineRecipe(self.config)

    def test_failed_stage_retains_products(self):
        self.config["steps"][3]["starts"] = [dict(ALPHA=99, METALS=0)]
        result = PipelineRecipe(
            self.config, base_dir=EXAMPLES, executable=self.exe
        ).fit(self.spectrum())
        self.assertFalse(result.success)
        self.assertIn("4", result.errors)
        self.assertEqual(result.element_fits, {})
        self.assertEqual(result.history[-1]["status"], "failed")

    def test_dr17_template_parses_without_production_grids(self):
        PipelineRecipe.read(EXAMPLES / "pipeline-dr17-style.json")

    def test_bundled_example_and_weight_file(self):
        from aspcap import windows

        spec = self.spectrum()
        bundled = PipelineRecipe.example("basic", executable=self.exe).fit(spec)
        self.assertTrue(bundled.success, bundled.errors)
        with tempfile.TemporaryDirectory() as root:
            self.config["grids"]["toy"] = str(EXAMPLES / "toy-grid.dat")
            weights = np.array([0, 0, 0.2, 0.8, 1, 0.5, 0, 0])
            windows.save_weights(Path(root) / "Mg.txt", spec.wave, weights)
            self.config["steps"][-1]["elements"][0]["weights"] = {"file": "Mg.txt"}
            result = PipelineRecipe(
                self.config, base_dir=root, executable=self.exe
            ).fit(spec)
            self.assertTrue(result.success, result.errors)
            np.testing.assert_array_equal(
                [p[-1] for p in result.element_fits["Mg"].pixels], weights
            )

    def test_initial_override_and_stale_fit_rejection(self):
        recipe = PipelineRecipe(self.config, base_dir=EXAMPLES, executable=self.exe)
        result = recipe.fit(self.spectrum(), initial={"ALPHA": 0.1, "METALS": -0.25})
        self.assertTrue(result.success, result.errors)
        self.assertEqual(result.history[2]["parameters"]["ALPHA"], 0.1)
        self.config["steps"].append(dict(type="normalize", method="skip"))
        with self.assertRaises(ValueError):
            PipelineRecipe(self.config)

    def test_cli_basic(self):
        self.spectrum()
        with tempfile.TemporaryDirectory() as root:
            output = Path(root) / "result.json"
            run = subprocess.run(
                [
                    sys.executable,
                    "-m",
                    "aspcap.cli",
                    "run",
                    str(EXAMPLES / "pipeline-basic.json"),
                    str(EXAMPLES / "pipeline-spectrum.csv"),
                    "--ferre",
                    self.exe,
                    "--output",
                    str(output),
                ],
                capture_output=True,
                text=True,
            )
            self.assertEqual(run.returncode, 0, run.stderr)
            result = json.loads(output.read_text())
            self.assertTrue(result["success"])
            self.assertAlmostEqual(result["abundances"]["Mg"], -0.05, delta=0.003)
