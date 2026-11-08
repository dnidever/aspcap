"""Real FERRE comparisons for storage formats and native output geometry."""

import math
import os
from pathlib import Path
import shutil
import struct
import tempfile
import unittest

from aspcap import Grid, fit_abundance


class GridFormatsTest(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.exe = os.environ.get("ASPCAP_TEST_FERRE") or shutil.which("ferre")

    def require_ferre(self):
        if not self.exe:
            self.skipTest("Install FERRE to test actual library evaluation")

    def library(self, *, pca=False, binary=False, transposed=False, multi=True):
        header = (
            "&SYNTH ID='test / quoted ! header', N_OF_DIM=1, LABEL(1)='X', "
            "N_P=3, LLIMITS=0, STEPS=1, "
            f"NPIX={4 if pca else 8}, MULTI={2 if multi else 0}, "
            f"TRANSPOSED={int(transposed)}, "
        )
        if pca:
            header += "NPCA=2*4, CONSTANT=10, "
        if not multi:
            header += "WAVE=15000,1, LOGW=0, VACUUM=1, RESOLUTION=22500, "
        header += "/\n"
        if multi:
            header += (
                "&SYNTH NPIX=4, WAVE=15000,1, LOGW=0, VACUUM=1, RESOLUTION=22500 /\n"
                "&SYNTH NPIX=4, WAVE=16000,2, LOGW=0, VACUUM=1, RESOLUTION=22500 /\n"
            )
        self.mean = [0.8] * 8
        self.basis = [[0.1, -0.2, 0.3, -0.4] * 2, [0.2, 0.1, -0.1, -0.2] * 2]
        if pca:
            for row in [self.mean, [1] * 8, *self.basis]:
                header += " ".join(map(str, row)) + "\n"
            rows = [[10 + x, 10 + 0.1 * x, 10 + 2 * x, 10 + 0.2 * x] for x in range(3)]
        else:
            rows = [self.expected(x) for x in range(3)]
        if transposed:
            rows = list(zip(*rows))
        path = self.root / ("grid.hdr" if binary else "grid.dat")
        path.write_text(header)
        if binary:
            values = [v for row in rows for v in row]
            path.with_suffix(".unf").write_bytes(
                struct.pack(f"={len(values)}f", *values)
            )
        else:
            with path.open("a") as out:
                out.write("".join(" ".join(map(str, row)) + "\n" for row in rows))
        return Grid(path, executable=self.exe or "ferre")

    def expected(self, x):
        b1 = [0.1, -0.2, 0.3, -0.4] * 2
        b2 = [0.2, 0.1, -0.1, -0.2] * 2
        return [0.8 + x * (1 if i < 4 else 2) * (b1[i] + 0.1 * b2[i]) for i in range(8)]

    def test_pca_headers_and_expanded_pixel_count(self):
        grid = self.library(pca=True, binary=True)
        self.assertTrue(grid.compressed)
        self.assertEqual(grid.npix, 8)
        self.assertEqual(grid.stored_npix, (4,))
        self.assertEqual(
            grid.wavelengths, (15000, 15001, 15002, 15003, 16000, 16002, 16004, 16006)
        )
        self.assertEqual([(s.start, s.stop) for s in grid.segments], [(0, 4), (4, 8)])

    def test_real_pca_ascii_binary_model_fit_and_abundance(self):
        self.require_ferre()
        for binary, direct in ((False, False), (True, False), (True, True)):
            grid = self.library(pca=True, binary=binary)
            if direct:
                grid = Grid(grid.path, executable=self.exe, f_access=1)
            with self.subTest(binary=binary, direct=direct):
                model = grid({"X": 0.4}).flux
                for got, want in zip(model, self.expected(0.4)):
                    self.assertAlmostEqual(got, want, places=5)
                fit = grid.fit(self.expected(0.4), [0.01] * 8, start={"X": 0.8})
                self.assertAlmostEqual(fit.parameters["X"], 0.4, places=3)
                self.assertEqual(len(fit.model), 8)
                result, _ = fit_abundance(
                    grid,
                    grid.wavelengths,
                    self.expected(0.4),
                    [0.01] * 8,
                    [1] * 8,
                    {"X": 0.8},
                    star_id="pca",
                    species="X",
                    label="X",
                    windows=[(15000, 16006)],
                    executable=self.exe,
                )
                self.assertAlmostEqual(result["abundance"], 0.4, places=3)
                # All-zero raw coefficients relative to CONSTANT reconstruct the mean.
                for got in grid({"X": 0}).flux:
                    self.assertAlmostEqual(got, 0.8, places=5)

    def test_multi_transposed_and_no_interpolation_across_gaps(self):
        self.require_ferre()
        for binary in (False, True):
            grid = self.library(binary=binary, transposed=True)
            flux = grid({"X": 0.3}).flux
            for got, want in zip(flux, self.expected(0.3)):
                self.assertAlmostEqual(got, want, places=5)
            self.assertEqual(len(grid({"X": 0.3}, wave=[16001, 15000.5]).flux), 2)
            with self.assertRaisesRegex(ValueError, "gap"):
                grid({"X": 0.3}, wave=[15500])

    def test_separate_synthfiles_and_formatted_direct_access(self):
        self.require_ferre()
        a = self.library(multi=False)
        other = self.root / "other.dat"
        shutil.copy(a.path, other)
        joined = Grid([a.path, other], executable=self.exe)
        self.assertEqual(len(joined({"X": 0.4}).flux), 16)
        h = self.root / "formatted.hdr"
        h.write_text(a.path.read_text().split("/\n", 1)[0] + "/\n")
        rows = [self.expected(x) for x in range(3)]
        h.with_suffix(".fmt").write_text(
            "".join("".join(f"{v:12.5E}" for v in row) for row in rows)
        )
        formatted = Grid(h, executable=self.exe)
        self.assertEqual((formatted.f_format, formatted.f_access), (0, 1))
        for got, want in zip(formatted({"X": 0.4}).flux, self.expected(0.4)):
            self.assertAlmostEqual(got, want, places=5)

    def test_sparse_sidecar_nonuniform_axes_and_invalid_region(self):
        self.require_ferre()
        folder = self.root / "library"
        folder.mkdir()
        path = folder / "sparse.dat"
        path.write_text(
            "&SYNTH N_OF_DIM=2,LABEL(1)='T',LABEL(2)='G',NODEFILE='points.nodes',"
            "NPIX=8,WAVE=15000,1,LOGW=0,VACUUM=1,RESOLUTION=22500 /\n"
            + "".join(
                " ".join(map(str, [1 + 0.0001 * (t - 3000), 0.8 + 0.1 * g, 1, 1] * 2))
                + "\n"
                for t, g in [(3000, 0), (3000, 1), (3100, 0), (3100, 1), (3300, 1)]
            )
        )
        (folder / "points.nodes").write_text("3000 0\n3000 1\n3100 0\n3100 1\n3300 1\n")
        grid = Grid(path, executable=self.exe)
        self.assertTrue(grid.sparse)
        self.assertEqual(grid.nspectra, 5)
        self.assertEqual(grid.axis_values[0], (3000, 3100, 3300))
        out = grid({"T": 3050, "G": 0.5}).flux
        self.assertAlmostEqual(out[0], 1.005, places=5)
        self.assertAlmostEqual(out[1], 0.85, places=5)
        fit = grid.fit(out, [0.01] * 8, start={"T": 3040, "G": 0.4})
        self.assertAlmostEqual(fit.parameters["T"], 3050, delta=1)
        with self.assertRaises(RuntimeError):
            grid({"T": 3200, "G": 0.5})

    def test_natural_log_wavelength_and_photometric_metadata(self):
        grid = self.library(multi=False)
        s = grid.path.read_text().replace(
            "WAVE=15000,1, LOGW=0", f"WAVE={math.log(15000)},{.0001}, LOGW=2"
        )
        grid.path.write_text(s)
        grid = Grid(grid.path)
        self.assertAlmostEqual(grid.wavelengths[0], 15000)
        s = s.replace("RESOLUTION=22500", "RESOLUTION=0")
        grid.path.write_text(s)
        grid = Grid(grid.path)
        self.assertEqual(grid.wavelengths, (0.0,) * 8)

    def test_singleton_sparse_axis_and_photometric_evaluation(self):
        self.require_ferre()
        grid = self.library(multi=False)
        grid.path.write_text(
            grid.path.read_text().replace("RESOLUTION=22500", "RESOLUTION=0")
        )
        photo = Grid(grid.path, executable=self.exe)
        self.assertEqual(len(photo({"X": 0.4}).flux), 8)
        path = self.root / "singleton.dat"
        path.write_text(
            "&SYNTH N_OF_DIM=2,LABEL(1)='X',LABEL(2)='FIXED',"
            "NODEFILE='singleton.nodes',NPIX=4,WAVE=15000,1,LOGW=0,"
            "VACUUM=1,RESOLUTION=22500 /\n1 1 1 1\n1 .8 1 1\n"
        )
        path.with_suffix(".nodes").write_text("0 42\n1 42\n")
        sparse = Grid(path, executable=self.exe)
        self.assertEqual(sparse.counts, (2, 1))
        self.assertAlmostEqual(
            sparse({"X": 0.5, "FIXED": 42}).flux[1], 0.9, places=5
        )
        fit = sparse.fit(
            [1, 0.9, 1, 1], [0.01] * 4, start={"X": 0.4, "FIXED": 42}, free_labels="X"
        )
        self.assertAlmostEqual(fit.parameters["X"], 0.5, places=3)
