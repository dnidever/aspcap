"""Pipeline stages agree with explicit component calls and preserve partial results."""

import os
from pathlib import Path
import shutil
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
from numpy.testing import assert_allclose, assert_array_equal

from aspcap import Pipeline, PipelineResult, ElementConfig, Grid, Spec1D
from aspcap.ferre import FerreFit

EXAMPLES = Path(__file__).resolve().parents[1] / "examples"
ELEMENT = dict(
    species="Mg", label="ALPHA", windows=[(15002, 15005)], add_labels=["METALS"]
)
TRUTH = dict(ALPHA=0.2, METALS=-0.25)


class PipelineTest(unittest.TestCase):
    def setUp(self):
        self.exe = os.environ.get("ASPCAP_TEST_FERRE") or shutil.which("ferre")
        self.grid = Grid(EXAMPLES / "toy-grid.dat", executable=self.exe or "ferre")
        self.spec = Spec1D(
            np.ones(8) * 2,
            wave=self.grid.wavelengths,
            err=np.ones(8) * 0.02,
            metadata={"star_id": "one-star"},
        )

    def synthetic(self):
        if not self.exe:
            self.skipTest("Install FERRE for full pipeline integration")
        model = self.grid(TRUTH)
        model.flux *= 2
        model.err = np.full(8, 0.02)
        return model

    def test_real_pipeline_matches_explicit_steps_and_does_not_double_normalize(self):
        spec = self.synthetic()
        original = spec.copy()
        p = Pipeline(
            self.grid,
            [ELEMENT],
            initial={"ALPHA": 0, "METALS": -0.25},
            continuum=np.full(8, 2.0),
        )
        result = p.fit(spec, star_id="synthetic")
        self.assertIsInstance(result, PipelineResult)
        self.assertTrue(result.success, result.errors)
        self.assertAlmostEqual(result.parameters["ALPHA"], 0.2, places=3)
        self.assertAlmostEqual(result.parameters["METALS"], -0.25, delta=0.01)
        self.assertAlmostEqual(result.abundances["Mg"], -0.05, delta=0.01)
        prepared = p.prepare(spec)
        normalized, estimate = p.normalize(prepared)
        initial = p.estimate_initial(normalized)
        glob = p.fit_global(normalized, initial)
        elems = p.fit_elements(normalized, glob.parameters, star_id="synthetic")
        self.assertEqual(result.parameters, glob.parameters)
        self.assertEqual(result.element_fits["Mg"].result, elems["Mg"].result)
        assert_allclose(result.continuum.continuum, 2.0)
        assert_allclose(result.spectrum.err, 0.01)
        assert_array_equal(spec.flux, original.flux)
        assert_array_equal(spec.err, original.err)
        self.assertFalse(spec.normalized)
        self.assertIsNot(result.original_spectrum, spec)
        self.assertEqual(result.residuals.shape, (8,))
        self.assertLess(np.max(np.abs(result.residuals)), 0.001)
        result.spectrum.flux[0] = 99
        assert_array_equal(spec.flux, original.flux)

    def test_fractional_element_weights_and_global_exclusions(self):
        spec = self.synthetic()
        weights = np.array([0, 0, 0.2, 1, 0.7, 0.3, 0, 0])
        config = ElementConfig("Mg", "ALPHA", weights=weights, add_labels=("METALS",))
        weights[2] = 99  # Configuration holds its own immutable copy.
        p = Pipeline(self.grid, [config], initial=TRUTH, continuum=np.full(8, 2.0))
        global_weights = np.ones(8) * 4
        global_weights[3] = 0
        result = p.fit(spec, weights=global_weights)
        self.assertTrue(result.success, result.errors)
        actual = [row[-1] for row in result.element_fits["Mg"].pixels]
        assert_allclose(actual, [0, 0, 0.2, 0, 0.7, 0.3, 0, 0])
        self.assertEqual(result.element_fits["Mg"].result["pixels_used"], 3)
        self.assertEqual(result.config["elements"][0]["weights"], config.weights)

    def test_mapping_elements_skip_normalization_and_global_only(self):
        p = Pipeline(
            self.grid,
            {"Mg": {k: v for k, v in ELEMENT.items() if k != "species"}},
            initial=TRUTH,
        )
        self.assertEqual(p.elements[0].species, "Mg")
        copy, estimate = p.normalize(self.spec)
        self.assertIsNone(estimate)
        assert_array_equal(copy.flux, self.spec.flux)
        fake = FerreFit(TRUTH, {"ALPHA": 0.01, "METALS": 0.01}, 1.0, tuple([1] * 8))
        only = Pipeline(self.grid, initial=TRUTH)
        with patch.object(self.grid, "fit", return_value=fake):
            result = only.fit(self.spec)
        self.assertTrue(result.success)
        self.assertEqual(result.abundances, {})
        self.assertIsNone(result.continuum)

    def test_partial_element_failure_preserves_success(self):
        spec = self.synthetic()
        p = Pipeline(
            self.grid,
            [
                ELEMENT,
                dict(species="Uncovered", label="ALPHA", windows=[(17000, 17001)]),
            ],
            initial=TRUTH,
            continuum=np.full(8, 2.0),
        )
        result = p.fit(spec)
        self.assertFalse(result.success)
        self.assertIsNotNone(result.global_fit)
        self.assertAlmostEqual(result.abundances["Mg"], -0.05, delta=0.01)
        self.assertIsNone(result.abundances["Uncovered"])
        self.assertIn("PARTIAL_ABUNDANCES", result.flags)
        self.assertIn("element:Uncovered", result.errors)
        self.assertEqual(
            result.element_fits["Uncovered"].flags, ("ELEMENT_FIT_FAILED",)
        )
        self.assertGreater(len(result.element_fits["Mg"].pixels), 0)
        repeated = p.fit_elements(result.spectrum, result.parameters, species="Mg")
        self.assertEqual(set(repeated), {"Mg"})
        self.assertTrue(repeated["Mg"].success)

    def test_fatal_global_failure_stops_abundances_and_can_reraise(self):
        p = Pipeline(self.grid, [ELEMENT], initial=TRUTH)
        with patch.object(
            self.grid, "fit", side_effect=RuntimeError("engine failed")
        ), patch.object(self.grid, "fit_abundance") as element:
            result = p.fit(self.spec)
            element.assert_not_called()
            self.assertIn("GLOBAL_FIT_FAILED", result.flags)
            self.assertIsNone(result.parameters)
            self.assertEqual(result.element_fits, {})
            self.assertIsNotNone(result.spectrum)
            with self.assertRaisesRegex(RuntimeError, "engine failed"):
                p.fit(self.spec, raise_on_error=True)

    def test_preparation_initial_and_continuum_errors_are_identified(self):
        p = Pipeline(self.grid, initial=TRUTH)
        bad = self.spec.copy()
        bad.wave += 1
        result = p.fit(bad)
        self.assertIn("PREPARATION_FAILED", result.flags)
        self.assertIsNone(result.global_fit)
        result = Pipeline(self.grid).fit(self.spec)
        self.assertIn("INITIAL_ESTIMATES_FAILED", result.flags)
        with self.assertRaisesRegex(ValueError, "initial"):
            Pipeline(self.grid).fit(self.spec, raise_on_error=True)
        result = p.fit(self.spec, continuum=np.ones(7))
        self.assertIn("NORMALIZATION_FAILED", result.flags)
        with self.assertRaises(TypeError):
            p.fit([1, 2, 3])

    def test_global_zero_weights_are_carried_into_element_masks(self):
        p = Pipeline(
            self.grid, [ELEMENT], initial=TRUTH, weights=[0, 1, 1, 1, 1, 1, 1, 1]
        )
        row = dict(abundance=-0.05, flags="")
        with patch.object(
            self.grid, "fit_abundance", return_value=(row, iter(()))
        ) as fit:
            out = p.fit_elements(self.spec, TRUTH)
        self.assertTrue(out["Mg"].success)
        self.assertTrue(fit.call_args.kwargs["mask"][0])
        assert_array_equal(fit.call_args.args[3], np.ones(8))

    def test_normalization_options_and_per_spectrum_overrides(self):
        p = Pipeline(
            self.grid, initial=TRUTH, continuum={"method": "chebyshev", "degree": 0}
        )
        normalized, estimate = p.normalize(self.spec)
        assert_allclose(estimate.continuum, 2.0)
        assert_allclose(normalized.flux, 1.0)
        assert_allclose(normalized.err, 0.01)
        skipped, estimate = p.normalize(self.spec, continuum=None)
        self.assertIsNone(estimate)
        assert_allclose(skipped.flux, 2.0)
        self.assertEqual(
            p.estimate_initial(self.spec, initial=(0.1, -0.1)),
            dict(ALPHA=0.1, METALS=-0.1),
        )

    def test_orders_and_native_gap_mapping_survive_pipeline(self):
        model = self.synthetic()
        spec = Spec1D.from_native(
            model.wave,
            model.flux,
            err=model.err,
            regions=[(15000, 15002), (15005, 15007)],
        )
        p = Pipeline(
            self.grid,
            [ELEMENT],
            initial=TRUTH,
            free_labels=["ALPHA"],
            continuum=[np.full(3, 2.0), np.full(3, 2.0)],
        )
        result = p.fit(spec)
        self.assertTrue(result.success, result.errors)
        self.assertEqual(result.spectrum.norder, 2)
        native = result.spectrum.to_native()
        self.assertTrue(native.mask[3])
        self.assertTrue(native.mask[4])
        self.assertEqual(result.element_fits["Mg"].result["pixels_used"], 2)
        self.assertAlmostEqual(result.abundances["Mg"], -0.05, places=3)

    def test_static_configuration_errors(self):
        with self.assertRaisesRegex(ValueError, "unique"):
            Pipeline(self.grid, [ELEMENT, ELEMENT])
        with self.assertRaisesRegex(ValueError, "labels"):
            Pipeline(self.grid, [dict(ELEMENT, label="UNKNOWN")])
        with self.assertRaises(ValueError):
            ElementConfig("Mg", "ALPHA", [(2, 1)])
        with self.assertRaises(ValueError):
            Pipeline(self.grid, weights=[1])

    def test_supplied_mini_grid_uses_its_dedicated_label(self):
        spec = self.synthetic()
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "mini.dat"
            path.write_text(
                (EXAMPLES / "toy-grid.dat").read_text().replace("'ALPHA'", "'MG_M'")
            )
            mini = Grid(path, executable=self.exe)
            path_config = ElementConfig(
                "Mg", "MG_M", [(15002, 15005)], ("METALS",), path, 0.0
            )
            path_pipeline = Pipeline(self.grid, [path_config])
            self.assertEqual(
                path_pipeline.elements[0].grid.executable, self.grid.executable
            )
            cfg = ElementConfig("Mg", "MG_M", [(15002, 15005)], ("METALS",), mini, 0.0)
            pipeline = Pipeline(
                self.grid, [cfg], initial=TRUTH, continuum=np.full(8, 2.0)
            )
            result = pipeline.fit(spec)
            self.assertTrue(result.success, result.errors)
            self.assertEqual(result.element_fits["Mg"].result["grid_label"], "MG_M")
            self.assertAlmostEqual(result.abundances["Mg"], -0.05, delta=0.01)
            missing = ElementConfig("Mg", "MG_M", [(15002, 15005)], ("METALS",), mini)
            result = Pipeline(
                self.grid, [missing], initial=TRUTH, continuum=np.full(8, 2.0)
            ).fit(spec)
            self.assertFalse(result.success)
            self.assertIsNotNone(result.global_fit)
            self.assertIn("initial", result.element_fits["Mg"].error)
