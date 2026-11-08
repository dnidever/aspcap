import os
from pathlib import Path
import shutil
import unittest
import numpy as np
from numpy.testing import assert_allclose, assert_array_equal
from aspcap import Spec1D, continuum, Grid


class ContinuumTests(unittest.TestCase):
    def test_per_order_polynomial_recovery_and_error_scaling(self):
        curves = [
            2 + np.linspace(-1, 1, n) * 0.2 + 0.1 * np.linspace(-1, 1, n) ** 4
            for n in (100, 77, 64)
        ]
        waves = [
            15000 + np.arange(n) * 0.2 + i * 100 for i, n in enumerate((100, 77, 64))
        ]
        spec = Spec1D(
            curves, wave=waves, err=[np.ones(n) * 0.02 for n in (100, 77, 64)]
        )
        result = continuum.fit(spec, degree=4)
        out = result.apply(spec)
        for i, n in enumerate(spec.numpix):
            assert_allclose(out.flux[:n, i], 1, atol=1e-12)
            assert_allclose(out.err[:n, i], 0.02 / curves[i])
            assert_allclose(result.continuum[:n, i], curves[i])
        self.assertTrue(np.all(np.isnan(result.continuum[77:, 1])))
        self.assertFalse(np.any(result.fit_mask[77:, 1]))
        assert_allclose(spec.flux[:, 0], curves[0])
        self.assertEqual(out.continuum_history[0]["settings"]["degree"], 4)

    def test_weighting_and_clipping(self):
        x = np.linspace(-1, 1, 101)
        truth = 2 + 0.2 * x
        flux = truth.copy()
        flux[50] += 10
        err = np.ones(101) * 0.01
        err[50] = 100
        spec = Spec1D(flux, wave=x + 10, err=err)
        plain = continuum.fit(spec, degree=1)
        weighted = continuum.fit(spec, degree=1, weighted=True)
        self.assertLess(
            np.max(abs(weighted.continuum - truth)),
            np.max(abs(plain.continuum - truth)) * 0.001,
        )
        clipped = continuum.fit(spec, degree=1, clip=(3, 3))
        self.assertFalse(clipped.fit_mask[50])
        assert_allclose(clipped.continuum, truth, atol=1e-12)
        with self.assertRaises(ValueError):
            continuum.fit(Spec1D(flux), weighted=True)

    def test_spline_and_wavelength_coordinates(self):
        wave = np.linspace(15000, 16000, 80)[::-1]
        x = (wave - 15500) / 500
        truth = 2 + 0.2 * x + 0.1 * x * x - 0.03 * x**3
        spec = Spec1D(truth, wave=wave)
        out = continuum.fit(spec, method="spline", coordinate="wavelength", smoothing=0)
        assert_allclose(out.continuum, truth, atol=1e-12)
        with self.assertRaises(ValueError):
            continuum.fit(
                Spec1D([1, 2, 3, 4], wave=[1, 1, 2, 3]),
                method="spline",
                coordinate="wavelength",
            )

    def test_filters_and_model_ratio(self):
        wave = np.arange(301.0) + 15000
        x = np.linspace(0, 2 * np.pi, 301)
        reference = 1 - 0.3 * np.exp(-(((wave - 15150) / 5) ** 2))
        truth = 2 + 0.05 * np.sin(x)
        flux = reference * truth
        flux[70] = np.nan
        spec = Spec1D(flux, wave=wave, err=np.ones(301) * 0.01)
        result = continuum.fit(spec, method="model_ratio", model=reference, window=11)
        self.assertFalse(result.fit_mask[70])
        assert_allclose(result.continuum[10:-10], truth[10:-10], atol=0.0004)
        self.assertTrue(result.apply(spec).mask[70])
        flat = Spec1D(np.ones(100) * 3, wave=np.arange(100) + 1.0)
        for method in ("median", "percentile"):
            assert_allclose(continuum.fit(flat, method=method, window=7).continuum, 3)
        model = Spec1D(reference, wave=wave + 1)
        with self.assertRaises(ValueError):
            continuum.fit(spec, method="model_ratio", model=model)
        with self.assertRaises(ValueError):
            continuum.fit(spec, method="median", window=4)
        with self.assertRaises(ValueError):
            continuum.fit(spec, method="model_ratio")

    def test_failure_flags_and_invalid_scale_application(self):
        spec = Spec1D(
            [np.ones(10), np.ones(6)],
            wave=[np.arange(10) + 1.0, np.arange(6) + 20.0],
            mask=[np.zeros(10, bool), np.ones(6, bool)],
        )
        with self.assertRaisesRegex(ValueError, "Order 1"):
            continuum.fit(spec)
        result = continuum.fit(spec, on_failure="mask")
        self.assertIn("INSUFFICIENT_PIXELS", result.flags[1])
        self.assertTrue(np.isnan(result.continuum[:, 1]).all())
        self.assertTrue(result.apply(spec).mask[:, 1].all())
        scale = np.ones(spec.shape)
        scale[2, 0] = 0
        out = continuum.from_array(spec, scale).apply(spec)
        self.assertTrue(out.mask[2, 0])
        self.assertTrue(np.isnan(out.flux[2, 0]))
        with self.assertRaises(ValueError):
            continuum.fit(spec, method="unknown")

    def test_edited_errors_are_still_excluded(self):
        spec = Spec1D(np.ones(20), wave=np.arange(20) + 1.0, err=np.ones(20) * 0.1)
        spec.err[3] = 0
        spec.flux[4] = np.nan
        result = continuum.fit(spec, weighted=True)
        self.assertFalse(result.fit_mask[3])
        self.assertFalse(result.fit_mask[4])
        self.assertTrue(result.apply(spec).mask[3])

    def test_apply_alignment_composition_and_no_hidden_conversion(self):
        spec = Spec1D(
            np.ones(30) * 6,
            wave=np.arange(30) + 15000.0,
            err=np.ones(30) * 0.6,
            wavevac=False,
            flux_unit="counts",
        )
        first = continuum.from_array(spec, np.ones(30) * 2).apply(spec)
        final = continuum.from_array(first, np.ones(30) * 3).apply(first)
        assert_allclose(final.flux, 1)
        assert_allclose(final.err, 0.1)
        assert_allclose(final.cont, 6)
        assert_allclose(final._flux, 6)
        self.assertEqual(final.flux_unit, "")
        self.assertEqual(final._flux_unit, "counts")
        self.assertFalse(final.wavevac)
        assert_array_equal(final.wave, spec.wave)
        shifted = spec.copy()
        shifted.wave += 1
        with self.assertRaises(ValueError):
            continuum.fit(spec).apply(shifted)

    def test_real_apogee_chip_mapping_and_fit(self):
        root = Path(__file__).resolve().parents[1] / "aspcap/data"
        stem = "apogee-dr17-mgiant-full"
        flux = (
            np.fromfile(root / (stem + ".unf"), dtype="<f4")
            .reshape(125, 8575)[62]
            .astype(float)
        )
        wave = 10 ** (4.179 + np.arange(8575) * 0.000006)
        regions = [
            (wave[a], wave[a + n - 1])
            for a, n in ((246, 3028), (3585, 2495), (6344, 1991))
        ]
        reference = Spec1D.from_native(wave, flux, regions=regions)
        raw = reference.copy()
        for i, n in enumerate(raw.numpix):
            raw.flux[:n, i] *= 2 + i
        raw.err = np.full(raw.shape, 0.003)
        raw.mask[10, 0] = True
        result = continuum.fit(raw, method="model_ratio", model=reference)
        normalized = result.apply(raw)
        native = normalized.to_native()
        good = ~native.mask
        assert_allclose(native.flux[good], flux[good], rtol=1e-12)
        self.assertEqual(native.npix, 8575)
        self.assertEqual(sum(good), 7513)
        executable = os.environ.get("ASPCAP_TEST_FERRE") or shutil.which("ferre")
        if executable:
            grid = Grid(root / (stem + ".hdr"), executable=executable)
            fit = grid.fit(
                native.flux,
                native.err,
                start=(-0.6, 0.8, 3550),
                weights=good.astype(float),
            )
            for label, expected, tol in [
                ("METALS", -0.5, 0.02),
                ("LOGG", 1, 0.03),
                ("TEFF", 3600, 10),
            ]:
                self.assertAlmostEqual(fit.parameters[label], expected, delta=tol)


if __name__ == "__main__":
    unittest.main()
