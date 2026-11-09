"""Check archived weights and their native/apStar pixel registration."""

import unittest
from importlib.resources import files

import numpy as np
from aspcap import windows


class WindowsTests(unittest.TestCase):
    def test_all_archived_species(self):
        self.assertEqual(len(windows.dr17_species()), 27)
        root = files("aspcap").joinpath("data/dr17")
        with root.joinpath("wave.dat").open() as stream:
            archived_wave = np.loadtxt(stream)[:, 0]
        for species in windows.dr17_species():
            for kind in ("mask", "filt"):
                wave, weight = windows.load_dr17(species, kind=kind)
                self.assertEqual(weight.shape, (7514,))
                self.assertTrue(np.all(np.isfinite(weight) & (weight >= 0)))
                np.testing.assert_allclose(wave, archived_wave, rtol=0, atol=0.00051)
            self.assertEqual(windows.load_dr17(species, kind="wind").shape[1], 3)

    def test_apstar_registration_and_subset(self):
        full_wave = 10.0 ** (4.179 + np.arange(8575) * 6e-6)
        wave, weight = windows.load_dr17("mg")
        _, full = windows.load_dr17("Mg", wavelength=full_wave)
        indices = np.concatenate(
            [np.arange(246, 3274), np.arange(3585, 6080), np.arange(6344, 8335)]
        )
        np.testing.assert_allclose(full_wave[indices], wave, rtol=0, atol=1e-10)
        np.testing.assert_array_equal(full[indices], weight)
        self.assertEqual(np.count_nonzero(full), np.count_nonzero(weight))
        _, subset = windows.load_dr17("Mg", wavelength=wave[100:200])
        np.testing.assert_array_equal(subset, weight[100:200])
        with self.assertRaises(ValueError):
            windows.load_dr17("Mg", wavelength=wave + 0.02)

    def test_final_masks_retain_fractional_weights(self):
        _, mask = windows.load_dr17("Mg")
        _, original = windows.load_dr17("Mg", kind="filt")
        self.assertTrue(np.any((mask > 0) & (mask < 1)))
        self.assertTrue(np.any((original > 0) & (mask == 0)))
        mask[:] = 0
        self.assertTrue(np.any(windows.load_dr17("Mg")[1] > 0))


class CustomWeightsTests(unittest.TestCase):
    def test_table_roundtrip_and_subset(self):
        from io import StringIO

        wave = np.array([15000.0, 15001.0, 15002.0])
        stream = StringIO()
        windows.save_weights(stream, wave, [0, 0.3, 1])
        text = stream.getvalue()
        w, weights = windows.load_weights(StringIO(text))
        np.testing.assert_array_equal(w, wave)
        np.testing.assert_array_equal(weights, [0, 0.3, 1])
        _, subset = windows.load_weights(StringIO(text), wavelength=wave[[2, 0]])
        np.testing.assert_array_equal(subset, [1, 0])
        with self.assertRaises(ValueError):
            windows.load_weights(StringIO(text), wavelength=[15000.1])

    def test_legacy_requires_registration(self):
        from io import StringIO

        with self.assertRaises(ValueError):
            windows.load_weights(StringIO("0\n0.3\n1"))
        _, weight = windows.load_weights(StringIO("0\n0.3\n1"), wavelength=[1, 2, 3])
        np.testing.assert_array_equal(weight, [0, 0.3, 1])
        root = files("aspcap").joinpath("data/dr17")
        with root.joinpath("Mg.mask").open() as stream:
            wave, weight = windows.load_weights(stream, scale="dr17")
        np.testing.assert_array_equal(weight, windows.load_dr17("Mg")[1])
        np.testing.assert_array_equal(wave, windows.load_dr17("Mg")[0])

    def test_intervals_and_legacy_edits(self):
        from io import StringIO

        _, weight = windows.load_windows(
            StringIO("1 3\n2 4"), wavelength=[1, 2, 3, 4, 5]
        )
        np.testing.assert_array_equal(weight, [1, 1, 1, 1, 0])
        _, weight = windows.load_legacy_windows(
            StringIO("1 4 0\n2 5 2"), wavelength=[1, 2, 3, 4, 5], weights=[0.3] * 5
        )
        np.testing.assert_array_equal(weight, [0.3, 0, 1, 1, 0.3])
        with self.assertRaises(ValueError):
            windows.load_windows(StringIO("1 3 0"))

    def test_reject_invalid_and_ambiguous_tables(self):
        from io import StringIO

        for text in ["1 -1", "1 nan", "1 0\n1 1"]:
            with self.assertRaises(ValueError):
                windows.load_weights(StringIO(text))
        with self.assertRaises(ValueError):
            windows.load_weights(StringIO("1 0\n1.001 1"), wavelength=[1], atol=0.002)
        with self.assertRaises(ValueError):
            windows.load_legacy_windows(
                StringIO("1 3 9"), wavelength=[1, 2], weights=[1, 1]
            )
