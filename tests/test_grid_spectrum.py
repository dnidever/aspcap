"""Spectrum model returns and fitting retain native sampling and data masks."""

import os
from pathlib import Path
import shutil
import unittest

import numpy as np
from numpy.testing import assert_allclose, assert_array_equal

from aspcap import Grid, Spec1D

EXAMPLES = Path(__file__).resolve().parents[1] / "examples"


class GridSpectrumTest(unittest.TestCase):
    def setUp(self):
        self.exe = os.environ.get("ASPCAP_TEST_FERRE") or shutil.which("ferre")
        self.grid = Grid(EXAMPLES / "toy-grid.dat", executable=self.exe or "ferre")
        self.truth = {"ALPHA": 0.2, "METALS": -0.25}

    def require_ferre(self):
        if not self.exe:
            self.skipTest("Install FERRE for spectrum integration")

    def test_model_return_types_metadata_and_numeric_callback(self):
        self.require_ferre()
        spec = self.grid(self.truth)
        self.assertIsInstance(spec, Spec1D)
        self.assertIsNone(spec.err)
        self.assertEqual(spec.norder, 1)
        self.assertTrue(spec.wavevac)
        self.assertEqual(spec.metadata["parameters"], self.truth)
        assert_array_equal(spec.wave, self.grid.wavelengths)
        legacy = self.grid(self.truth, as_dict=True)
        self.assertIsInstance(legacy, dict)
        self.assertIsInstance(legacy["flux"], tuple)
        assert_array_equal(spec.flux, legacy["flux"])
        batch = self.grid([self.truth, self.truth])
        self.assertIsInstance(batch, list)
        self.assertTrue(all(isinstance(item, Spec1D) for item in batch))
        batch[0].flux[0] = 99
        self.assertNotEqual(batch[1].flux[0], 99)
        batch[0].metadata["parameters"]["ALPHA"] = 99
        self.assertEqual(batch[1].metadata["parameters"], self.truth)
        self.assertEqual(len(self.grid([self.truth])), 1)
        old_batch = self.grid([self.truth], as_dict=True)
        self.assertEqual(len(old_batch["flux"]), 1)
        self.assertIsInstance(self.grid.model(None, 0.2, -0.25), tuple)
        assert_array_equal(self.grid.model(None, 0.2, -0.25), spec.flux)
        sampled = self.grid(self.truth, wave=[15004.5, 15002.5])
        assert_array_equal(sampled.wave, [15004.5, 15002.5])

    def test_spec_fit_masks_errors_and_array_equivalence(self):
        self.require_ferre()
        spec = self.grid(self.truth)
        spec.err = np.full(self.grid.npix, 0.01)
        spec.mask[3] = True
        spec.flux[2] = np.nan
        spec.err[4] = np.nan
        before = [a.copy() for a in (spec.flux, spec.err, spec.mask, spec.wave)]
        weights = np.ones(self.grid.npix)
        weights[7] = 0
        weight_copy = weights.copy()
        fit = self.grid.fit(
            spec, start=self.truth, free_labels="ALPHA", weights=weights
        )
        native_weights = np.where(spec.valid, weights, 0)
        reference = self.grid.fit(
            spec.flux,
            spec.err,
            start=self.truth,
            free_labels="ALPHA",
            weights=native_weights,
        )
        self.assertEqual(fit.parameters, reference.parameters)
        self.assertAlmostEqual(fit.parameters["ALPHA"], 0.2, places=3)
        for a, b in zip((spec.flux, spec.err, spec.mask, spec.wave), before):
            assert_array_equal(a, b)
        assert_array_equal(weights, weight_copy)

    def test_orders_and_saved_native_mapping(self):
        self.require_ferre()
        model = self.grid(self.truth)
        errors = np.full(self.grid.npix, 0.01)
        mapped = Spec1D.from_native(
            model.wave, model.flux, err=errors, regions=[(15000, 15002), (15005, 15007)]
        )
        fit = self.grid.fit(mapped, start=self.truth, free_labels="ALPHA")
        self.assertAlmostEqual(fit.parameters["ALPHA"], 0.2, places=3)
        orders = Spec1D(
            [model.flux[:3], model.flux[3:]],
            wave=[model.wave[:3], model.wave[3:]],
            err=[errors[:3], errors[3:]],
        )
        self.assertEqual(orders.norder, 2)
        fit = self.grid.fit(orders, start=self.truth, free_labels="ALPHA")
        self.assertAlmostEqual(fit.parameters["ALPHA"], 0.2, places=3)

    def test_explicit_error_override_and_validation(self):
        spec = Spec1D(np.ones(8), wave=self.grid.wavelengths)
        with self.assertRaisesRegex(ValueError, "uncertainties"):
            self.grid.fit(spec)
        with self.assertRaisesRegex(ValueError, "errors"):
            self.grid.fit(np.ones(8))
        for update, message in [
            ("shift", "native grid"),
            ("frame", "frame"),
            ("unit", "Angstrom"),
            ("wave", "wavelengths"),
        ]:
            bad = Spec1D(np.ones(8), wave=self.grid.wavelengths, err=np.ones(8))
            if update == "shift":
                bad.wave += 0.1
            elif update == "frame":
                bad.wavevac = False
            elif update == "unit":
                bad.wave_unit = "nm"
            else:
                bad.wave = None
            with self.subTest(update=update), self.assertRaisesRegex(
                ValueError, message
            ):
                self.grid.fit(bad)
        if self.exe:
            model = self.grid(self.truth)
            fit = self.grid.fit(
                model, np.full(8, 0.01), start=self.truth, free_labels="ALPHA"
            )
            self.assertAlmostEqual(fit.parameters["ALPHA"], 0.2, places=3)

    def test_unknown_frame_metadata_is_preserved(self):
        spec = Spec1D([1, 1], wave=[15000, 15001], wavevac=None)
        self.assertIsNone(spec.wavevac)
        self.assertIsNone(spec.copy().wavevac)
