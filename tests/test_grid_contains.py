"""Coverage checks separate sparse geometry from optional FERRE validation."""

import os
from pathlib import Path
import shutil
import tempfile
import unittest

from aspcap import Grid


class GridContainsTest(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.exe = os.environ.get("ASPCAP_TEST_FERRE") or shutil.which("ferre")
        self.dense = Grid(
            Path(__file__).resolve().parents[1] / "examples/toy-grid.dat",
            executable="missing-executable-for-metadata-only-check",
        )
        path = self.root / "sparse.dat"
        path.write_text(
            "&SYNTH N_OF_DIM=2, LABEL(1)='X',LABEL(2)='Y',"
            "NODEFILE='points.nodes',NPIX=4,WAVE=15000,1,LOGW=0,"
            "VACUUM=1,RESOLUTION=22500 /\n"
            "1 1 1 1\n1 .9 1 1\n1 1 .9 1\n0\n1 .7 .8 1\n"
        )
        (self.root / "points.nodes").write_text("0 0\n0 1\n1 0\n1 1\n3 1\n")
        self.sparse = Grid(path, executable=self.exe or "ferre")

    def test_regular_bounds_mapping_vector_and_no_executable(self):
        for point in [(0.2, -0.25), (-0.5, -0.5), (0.5, 0)]:
            self.assertTrue(self.dense.contains(point))
        self.assertTrue(self.dense.contains({"ALPHA": 0.2, "METALS": -0.25}))
        for point in [(1, -0.25), (0.2, -1), (float("nan"), 0), (0, float("inf"))]:
            self.assertFalse(self.dense.contains(point))

    def test_malformed_input_is_an_error(self):
        for point in [{"TYPO": 0, "METALS": 0}, [0], [[0, 0]], [None, 0], 42]:
            with self.subTest(point=point), self.assertRaises(ValueError):
                self.dense.contains(point)

    def test_sparse_nonuniform_holes_edges_and_cached_lookup(self):
        for point in [(0.5, 0.5), (0, 0), (1, 0), (3, 1), (2, 1)]:
            self.assertTrue(self.sparse.contains(point), point)
        for point in [(2, 0.5), (3, 0), (4, 1)]:
            self.assertFalse(self.sparse.contains(point), point)
        # Geometry remains available without rereading the sidecar.
        (self.root / "points.nodes").unlink()
        self.assertTrue(self.sparse.contains((0.5, 0.5)))
        self.assertFalse(self.sparse.contains((2, 0.5)))

    def test_model_check_detects_zero_spectrum_exclusions(self):
        if not self.exe:
            self.skipTest("Install FERRE to check invalid spectra")
        self.assertTrue(self.sparse.contains((1, 1)))
        self.assertFalse(self.sparse.contains((1, 1), check_model=True))
        self.assertTrue(self.sparse.contains((0, 0), check_model=True))
        self.assertFalse(self.sparse.contains((2, 0.5), check_model=True))

    def test_model_check_does_not_hide_missing_executable(self):
        with self.assertRaisesRegex(ValueError, "executable"):
            self.dense.contains((0, -0.25), check_model=True)
