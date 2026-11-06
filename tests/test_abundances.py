import csv
import os
from pathlib import Path
import shutil
import struct
import tempfile
import unittest

from aspcap import Grid, fit_abundance


EXAMPLES = Path(__file__).resolve().parents[1] / "examples"


class AbundanceTest(unittest.TestCase):
    def setUp(self):
        self.grid = Grid.read(EXAMPLES / "toy-grid.dat")
        with (EXAMPLES / "toy-spectrum.csv").open() as file:
            rows = list(csv.DictReader(file))
        self.columns = [[float(row[key]) for row in rows]
                        for key in ("wavelength", "flux", "error", "continuum")]

    def fit(self, grid, windows=((15002, 15005),), **options):
        return fit_abundance(
            grid, *self.columns, {"ALPHA": 0, "METALS": -0.25},
            star_id="synthetic-star", species="Mg", label="ALPHA",
            windows=windows, add_labels=("METALS",), **options,
        )

    def test_rejects_off_grid_wavelengths(self):
        self.columns[0][0] += 0.1
        with self.assertRaisesRegex(ValueError, "native grid"):
            self.fit(self.grid)

    def test_rejects_uncovered_windows(self):
        with self.assertRaisesRegex(ValueError, "fewer than two usable pixels"):
            self.fit(self.grid, windows=((15020, 15030),))

    def test_real_ferre(self):
        executable = os.environ.get("ASPCAP_TEST_FERRE") or shutil.which("ferre.x")
        if not executable:
            self.skipTest("Set ASPCAP_TEST_FERRE to run the FERRE integration test")
        with tempfile.TemporaryDirectory() as root:
            header, spectra = self.grid.path.read_text().split("/\n", 1)
            binary = Path(root) / "toy.hdr"
            binary.write_text(header + "/\n")
            flux = tuple(map(float, spectra.split()))
            binary.with_suffix(".unf").write_bytes(struct.pack(f"={len(flux)}f", *flux))
            for grid in (self.grid, Grid.read(binary)):
                with self.subTest(binary=grid.binary):
                    result, pixels = self.fit(grid, executable=executable, run_root=root)
                    self.assertAlmostEqual(result["raw_value"], 0.2, places=3)
                    self.assertAlmostEqual(result["abundance"], -0.05, places=3)
                    self.assertGreater(result["conditional_error"], 0)
                    rows = list(pixels)
                    self.assertEqual([row[-1] for row in rows], [0, 0, 1, 1, 1, 1, 0, 0])
                    # This pixel depends on frozen METALS and must not pull the fit.
                    self.assertAlmostEqual(rows[1][3], 1.025, places=5)
            self.assertEqual(sorted(p.name for p in Path(root).iterdir()), ["toy.hdr", "toy.unf"])

    def test_grid_model_and_fit(self):
        executable = os.environ.get("ASPCAP_TEST_FERRE") or shutil.which("ferre.x")
        if not executable:
            self.skipTest("Set ASPCAP_TEST_FERRE to run the FERRE integration test")
        grid = Grid(self.grid.path, executable=executable)
        self.assertIn("15 spectra", repr(grid))
        batch = grid([{"ALPHA": 0.2, "METALS": -0.25}, (0.5, 0.0)])
        self.assertEqual(batch["wave"], grid.wavelengths)
        self.assertAlmostEqual(batch["flux"][0][1], 1.025, places=5)
        self.assertAlmostEqual(batch["flux"][1][3], 0.7, places=5)
        self.assertAlmostEqual(grid.model([15002.5], 0.2, -0.25)[0], 0.82, places=5)

        flux = list(grid((0.2, -0.15))["flux"])
        flux[3] -= 0.001
        flux[4] += 0.001
        fit = grid.fit(flux, [0.01] * grid.npix)
        self.assertAlmostEqual(fit.parameters["ALPHA"], 0.2, places=2)
        self.assertAlmostEqual(fit.parameters["METALS"], -0.15, delta=0.015)
