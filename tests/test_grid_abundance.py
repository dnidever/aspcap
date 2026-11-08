"""The method and compatibility function share one Grid.fit execution path."""

import csv
import os
from pathlib import Path
import shutil
import unittest
from unittest.mock import patch

from aspcap import Grid, fit_abundance
from aspcap.ferre import FerreFit

EXAMPLES = Path(__file__).resolve().parents[1] / "examples"


class GridAbundanceTest(unittest.TestCase):
    def setUp(self):
        self.exe = os.environ.get("ASPCAP_TEST_FERRE") or shutil.which("ferre")
        self.grid = Grid(
            EXAMPLES / "toy-grid.dat", executable=self.exe or "custom-ferre"
        )
        with (EXAMPLES / "toy-spectrum.csv").open() as handle:
            rows = list(csv.DictReader(handle))
        self.columns = [
            [float(row[key]) for row in rows]
            for key in ("wavelength", "flux", "error", "continuum")
        ]
        self.parameters = {"ALPHA": 0.0, "METALS": -0.25}
        self.options = dict(
            star_id="test",
            species="Mg",
            label="ALPHA",
            windows=[(15002, 15005)],
            add_labels=("METALS",),
        )

    def test_method_reuses_grid_fit_with_normalized_data_and_one_label(self):
        output = FerreFit(
            {"ALPHA": 0.2, "METALS": -0.25},
            {"ALPHA": 0.01, "METALS": None},
            1.0,
            tuple([1] * self.grid.npix),
        )
        with patch.object(self.grid, "fit", return_value=output) as fitter:
            result, pixels = self.grid.fit_abundance(
                *self.columns,
                self.parameters,
                mask=[False] * 2 + [True] + [False] * 5,
                keep_runs=True,
                **self.options
            )
        fitter.assert_called_once()
        args, kwargs = fitter.call_args
        self.assertEqual(args[2], self.parameters)
        self.assertEqual(kwargs["free_labels"], ("ALPHA",))
        self.assertEqual(kwargs["weights"], [0, 0, 0, 1, 1, 1, 0, 0])
        self.assertTrue(kwargs["keep_runs"])
        self.assertIsNone(kwargs["executable"])  # Grid.fit uses this grid's setting.
        self.assertAlmostEqual(args[0][3], self.columns[1][3] / self.columns[3][3])
        self.assertAlmostEqual(args[1][3], self.columns[2][3] / self.columns[3][3])
        self.assertAlmostEqual(result["abundance"], -0.05)
        self.assertEqual(result["pixels_used"], 3)
        self.assertEqual(len(list(pixels)), self.grid.npix)

    def test_executable_override_and_validation(self):
        output = FerreFit(
            {"ALPHA": 0.2, "METALS": -0.25},
            {"ALPHA": 0.01, "METALS": None},
            1.0,
            tuple([1] * self.grid.npix),
        )
        with patch.object(self.grid, "fit", return_value=output) as fitter:
            self.grid.fit_abundance(
                *self.columns,
                self.parameters,
                executable="/custom/ferre",
                **self.options
            )
        self.assertEqual(fitter.call_args.kwargs["executable"], "/custom/ferre")
        with self.assertRaisesRegex(ValueError, "mask"):
            self.grid.fit_abundance(
                *self.columns, self.parameters, mask=[True], **self.options
            )

    def test_real_method_and_function_agree_and_honor_grid_executable(self):
        if not self.exe:
            self.skipTest("Install FERRE for fitting integration")
        method, method_pixels = self.grid.fit_abundance(
            *self.columns, self.parameters, **self.options
        )
        function, function_pixels = fit_abundance(
            self.grid, *self.columns, self.parameters, **self.options
        )
        self.assertEqual(method, function)
        self.assertEqual(list(method_pixels), list(function_pixels))
        self.assertAlmostEqual(method["abundance"], -0.05, places=3)
