import unittest
from types import SimpleNamespace
import numpy as np
from numpy.testing import assert_allclose, assert_array_equal
from aspcap import Spec1D, continuum


class SpectrumTests(unittest.TestCase):
    def test_ragged_orders_masks_and_independence(self):
        flux = [np.ones(12), np.ones(9) * 2, np.ones(7) * 3]
        wave = [np.arange(n) + 15000 + i * 100 for i, n in enumerate((12, 9, 7))]
        spec = Spec1D(flux, wave=wave, err=[f * 0.1 for f in flux])
        self.assertEqual(spec.shape, (12, 3))
        assert_array_equal(spec.numpix, [12, 9, 7])
        self.assertTrue(np.all(spec.mask[9:, 1]))
        self.assertTrue(np.all(spec.mask[7:, 2]))
        self.assertEqual(spec.size, 28)
        assert_array_equal(np.asarray(spec), spec.flux)
        child = spec[1]
        child.flux[:] = -1
        assert_array_equal(spec.flux[:9, 1], 2)
        self.assertEqual(spec[2:6, 2].shape, (4,))
        self.assertEqual(spec[1:].norder, 2)
        clone = spec.copy()
        clone.mask[0, 0] = True
        self.assertFalse(spec.mask[0, 0])
        flux[0][:] = 99
        assert_array_equal(spec.flux[:, 0], 1)

    def test_padding_and_invalid_input(self):
        wave = np.array([[1, 10], [2, 11], [3, 0], [4, 0.0]])
        spec = Spec1D(
            np.ones((4, 2)),
            wave=wave,
            err=np.array([[0.1, 0.1], [0.1, 0.1], [0.1, 0], [np.nan, 0.1]]),
        )
        assert_array_equal(spec.numpix, [4, 2])
        self.assertTrue(spec.mask[3, 0])
        self.assertTrue(np.all(spec.mask[2:, 1]))
        with self.assertRaises(ValueError):
            Spec1D([1, 2], wave=[1])
        with self.assertRaises(ValueError):
            Spec1D([1, 2], numpix=[3])
        with self.assertRaises(ValueError):
            Spec1D([1, 2], bitmask=[0.1, 0.2])
        with self.assertRaises(ValueError):
            Spec1D(np.zeros((2, 2, 2)))
        self.assertIsNone(Spec1D([1, 2]).err)
        bits = np.array([2**63, 0], dtype=np.uint64)
        assert_array_equal(Spec1D([1, 2], bitmask=bits).bitmask, bits)

    def test_native_mapping_and_masks(self):
        wave = np.arange(30.0) + 15000
        spec = Spec1D.from_native(
            wave,
            np.arange(30.0) + 1,
            err=np.ones(30),
            bitmask=np.arange(30),
            regions=[(15001, 15008), (15012, 15018), (15022, 15026)],
        )
        assert_array_equal(spec.numpix, [8, 7, 5])
        restored = spec.to_native()
        good = ~restored.mask
        assert_allclose(restored.flux[good], np.arange(30.0)[good] + 1)
        assert_array_equal(restored.bitmask[good], np.arange(30)[good])
        self.assertTrue(np.all(np.isnan(restored.flux[~good])))
        self.assertTrue(np.all(np.isinf(restored.err[~good])))
        assert_array_equal(spec.flatten().to_native().mask, restored.mask)
        partial = spec[0].to_native()
        self.assertEqual(np.sum(~partial.mask), 8)
        with self.assertRaises(ValueError):
            Spec1D.from_native(
                wave, np.ones(30), regions=[(15001, 15005), (15005, 15008)]
            )
        with self.assertRaises(ValueError):
            Spec1D(np.ones(3)).to_native()

    def test_callback_options_cache_and_copy(self):
        calls = []

        def estimator(spec, scale=2):
            calls.append(scale)
            return np.full(spec.shape, scale)

        spec = Spec1D(np.ones(12) * 6, err=np.ones(12) * 0.6, wave=np.arange(12) + 1.0)
        spec.continuum_func = estimator
        out = spec.normalize(scale=3, inplace=False)
        self.assertEqual(calls, [3])
        assert_allclose(out.flux, 2)
        assert_allclose(out.err, 0.2)
        assert_allclose(spec.flux, 6)
        self.assertFalse(spec.normalized)
        assert_allclose(spec.cont, 2)
        self.assertEqual(calls, [3, 2])
        spec.continuum_func = lambda s: np.ones(s.shape) * 4
        assert_allclose(spec.cont, 4)
        out.normalize()  # already normalized: no-op without new options
        with self.assertRaises(ValueError):
            out.normalize(method="median", window=3)
        spec.normalize()
        self.assertTrue(spec.normalized)
        assert_allclose(spec._flux, 6)
        assert_allclose(spec._err, 0.6)

    def test_normalized_state_survives_adaptation_and_order_selection(self):
        source = Spec1D.from_native(
            np.arange(20.0) + 1,
            np.full(20, 6.0),
            err=np.full(20, 0.6),
            regions=[(1, 7), (12, 20)],
            flux_unit="counts",
        )
        normalized = continuum.from_array(source, np.full(source.shape, 2.0)).apply(
            source
        )
        for selected in (normalized[0], normalized.flatten(), normalized.to_native()):
            self.assertEqual(selected._flux_unit, "counts")
            self.assertEqual(selected.flux_unit, "")
            self.assertTrue(selected.normalized)
        foreign = SimpleNamespace(
            flux=normalized.flux,
            err=normalized.err,
            wave=normalized.wave,
            mask=normalized.mask,
            numpix=normalized.numpix,
            normalized=True,
            _cont=normalized.cont,
            _flux=normalized._flux,
            _err=normalized._err,
            flux_unit="",
            _flux_unit="counts",
        )
        adapted = Spec1D.from_spectrum(foreign)
        out = continuum.from_array(adapted, np.full(adapted.shape, 3.0)).apply(adapted)
        assert_allclose(out.flux[out.valid], 1.0)
        assert_allclose(out.cont[out.valid], 6.0)
        self.assertEqual(out._flux_unit, "counts")
        assert_allclose(adapted.cont[adapted.valid], 2.0)

    def test_foreign_spectrum_without_doppler(self):
        spec = SimpleNamespace(
            flux=np.ones((20, 3)) * 4,
            wave=np.tile(np.arange(20.0)[:, None] + 1, (1, 3)),
            err=np.ones((20, 3)) * 0.4,
            mask=np.zeros((20, 3), bool),
            numpix=np.array([20, 20, 20]),
            wavevac=True,
            normalized=False,
        )
        result = continuum.fit(spec, degree=0)
        out = result.apply(spec)
        self.assertIsInstance(out, SimpleNamespace)
        assert_allclose(out.flux, 1)
        assert_allclose(out.err, 0.1)
        assert_allclose(spec.flux, 4)


if __name__ == "__main__":
    unittest.main()
