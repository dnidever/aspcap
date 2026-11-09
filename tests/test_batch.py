"""Batch failure isolation, CLI lists, and apStar FITS pixel registration."""

import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
from types import SimpleNamespace
import unittest

import numpy as np
from astropy.io import fits
from aspcap import Grid, Pipeline, PipelineRecipe, Spec1D, read_apstar

ROOT = Path(__file__).resolve().parents[1]
EXAMPLES = ROOT / "examples"


class BatchTests(unittest.TestCase):
    def setUp(self):
        self.exe = os.environ.get("ASPCAP_TEST_FERRE") or shutil.which("ferre")

    def spectrum(self):
        if not self.exe:
            self.skipTest("FERRE unavailable")
        grid = Grid(EXAMPLES / "toy-grid.dat", executable=self.exe)
        spec = grid({"ALPHA": 0.2, "METALS": -0.25})
        spec.err = np.full(8, 0.01)
        return grid, spec

    def test_python_batch_stages_and_full_pipeline(self):
        grid, spec = self.spectrum()
        bad = spec.copy()
        bad.wave += 0.1
        pipeline = Pipeline(grid, initial={"ALPHA": 0, "METALS": -0.25})
        items = pipeline.run_many(
            [spec, bad, spec], step="prepare", star_ids=["a", "b", "c"]
        )
        self.assertEqual([v.success for v in items], [True, False, True])
        self.assertIsNot(items[0].value, spec)
        fits = pipeline.fit_many([spec, bad, spec])
        self.assertEqual([v.success for v in fits], [True, False, True])
        single = pipeline.run_many(
            [spec], step="fit_global", parameters=[{"ALPHA": 0, "METALS": -0.25}]
        )
        self.assertTrue(single[0].success)
        with self.assertRaises(ValueError):
            pipeline.fit_many([spec], star_ids=[])

    def test_recipe_batch_and_prefix(self):
        _, spec = self.spectrum()
        recipe = PipelineRecipe.example("basic", executable=self.exe)
        items = recipe.fit_many([spec, None, spec], stop_after=2)
        self.assertEqual([v.success for v in items], [True, False, True])
        self.assertIsNone(items[0].value.parameters)
        self.assertEqual(len(items[0].value.history), 2)
        fits = recipe.fit_many(
            [spec, spec], initial=[{"ALPHA": 0, "METALS": -0.25}] * 2
        )
        self.assertTrue(all(v.success for v in fits))

    def test_cli_list_relative_paths_and_failure_continuation(self):
        self.spectrum()
        with tempfile.TemporaryDirectory() as root:
            root = Path(root)
            shutil.copy(EXAMPLES / "pipeline-spectrum.csv", root / "good.csv")
            (root / "inputs.txt").write_text(
                "# spectra\ngood.csv\nmissing.fits\ngood.csv\n"
            )
            output = root / "out.json"
            run = subprocess.run(
                [
                    sys.executable,
                    "-m",
                    "aspcap.cli",
                    "run",
                    str(EXAMPLES / "pipeline-basic.json"),
                    "@" + str(root / "inputs.txt"),
                    "--ferre",
                    self.exe,
                    "--output",
                    str(output),
                ],
                capture_output=True,
                text=True,
            )
            self.assertNotEqual(run.returncode, 0)
            result = json.loads(output.read_text())
            self.assertEqual(
                [v["success"] for v in result["results"]], [True, False, True]
            )
            self.assertEqual(result["count"], 3)
            run = subprocess.run(
                [
                    sys.executable,
                    "-m",
                    "aspcap.cli",
                    "run",
                    str(EXAMPLES / "pipeline-basic.json"),
                    str(root / "good.csv") + "," + str(root / "good.csv"),
                    "--ferre",
                    self.exe,
                    "--stop-after",
                    "2",
                    "--output",
                    str(output),
                ],
                capture_output=True,
                text=True,
            )
            self.assertEqual(run.returncode, 0, run.stderr)
            self.assertTrue(json.loads(output.read_text())["success"])


class ApStarReaderTests(unittest.TestCase):
    def write_file(self, path):
        flux = np.ones((3, 8575))
        flux[2] *= 2
        error = np.full_like(flux, 0.1)
        flags = np.zeros(flux.shape, dtype=np.int32)
        flags[:, 5] = 1
        flags[:, 6] = 4096
        error[:, 7] = 0
        primary = fits.PrimaryHDU()
        primary.header["OBJECT"] = "test-star"
        image = fits.ImageHDU(flux)
        image.header["CRVAL1"] = 4.179
        image.header["CDELT1"] = 6e-6
        image.header["CRPIX1"] = 1
        fits.HDUList(
            [primary, image, fits.ImageHDU(error), fits.ImageHDU(flags)]
        ).writeto(path)

    def test_cli_fits_list_and_row_selection(self):
        with tempfile.TemporaryDirectory() as root:
            root = Path(root)
            path = root / "apStar-test.fits"
            self.write_file(path)
            config = {
                "version": 1,
                "grids": {
                    "apogee": str(ROOT / "aspcap/data/apogee-dr17-mgiant-full.hdr")
                },
                "steps": [
                    {"type": "prepare", "grid": "apogee"},
                    {"type": "normalize", "method": "skip"},
                    {
                        "type": "initial",
                        "parameters": {"METALS": 0, "LOGG": 2, "TEFF": 3600},
                    },
                    {"type": "global_fit"},
                ],
            }
            recipe = root / "pipeline.json"
            recipe.write_text(json.dumps(config))
            inputs = root / "spectra.txt"
            inputs.write_text("apStar-test.fits\napStar-test.fits\n")
            output = root / "results.json"
            process = subprocess.run(
                [
                    sys.executable,
                    "-m",
                    "aspcap.cli",
                    "run",
                    str(recipe),
                    "--input-list",
                    str(inputs),
                    "--fits-row",
                    "2",
                    "--stop-after",
                    "1",
                    "--output",
                    str(output),
                ],
                capture_output=True,
                text=True,
            )
            self.assertEqual(process.returncode, 0, process.stderr)
            result = json.loads(output.read_text())
            self.assertTrue(result["success"])
            self.assertEqual(result["results"][0]["star_id"], "test-star:row2")
            self.assertEqual(result["count"], 2)

    def test_rows_masks_and_native_subset(self):
        with tempfile.TemporaryDirectory() as root:
            path = Path(root) / "apStar-test.fits"
            self.write_file(path)
            spec = read_apstar(path, row=2)
            np.testing.assert_array_equal(spec.flux, 2)
            self.assertTrue(spec.mask[5])
            self.assertFalse(spec.mask[6])
            self.assertTrue(spec.mask[7])
            self.assertEqual(spec.bitmask[6], 4096)
            self.assertEqual(spec.metadata["star_id"], "test-star:row2")
            selected = np.array([1, 5, 6, 7, 50])
            grid = SimpleNamespace(vacuum=True, wavelengths=spec.wave[selected])
            native = read_apstar(path, grid=grid)
            np.testing.assert_allclose(native.wave, grid.wavelengths)
            np.testing.assert_array_equal(native.bitmask, spec.bitmask[selected])
            with self.assertRaises(ValueError):
                read_apstar(path, row=99)
            grid.wavelengths = grid.wavelengths + 0.02
            with self.assertRaises(ValueError):
                read_apstar(path, grid=grid)
            allmasked = read_apstar(path, bad_bits=4096)
            self.assertTrue(allmasked.mask[6])
            self.assertFalse(allmasked.mask[5])
