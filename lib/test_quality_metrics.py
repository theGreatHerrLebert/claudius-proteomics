#!/usr/bin/env python3
"""Tests for quality_metrics.py - moments, Gaussian fitting, isotope scoring."""

import unittest
import numpy as np

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))

from quality_metrics import (
    calculate_moments,
    gaussian,
    fit_gaussian,
    _simple_isotope_score,
    compute_quality_summary,
)


class TestGaussianFunction(unittest.TestCase):
    """Test the Gaussian function itself."""

    def test_peak_at_center(self):
        """Peak value equals amplitude at mu."""
        result = gaussian(np.array([5.0]), 100.0, 5.0, 1.0)
        self.assertAlmostEqual(result[0], 100.0, places=5)

    def test_symmetric(self):
        """Gaussian is symmetric around mu."""
        x = np.array([4.0, 6.0])
        result = gaussian(x, 100.0, 5.0, 1.0)
        self.assertAlmostEqual(result[0], result[1], places=5)

    def test_fwhm(self):
        """At +/- sigma, value should be ~60.65% of peak."""
        result = gaussian(np.array([4.0]), 100.0, 5.0, 1.0)
        self.assertAlmostEqual(result[0], 100.0 * np.exp(-0.5), places=3)

    def test_zero_amplitude(self):
        result = gaussian(np.array([5.0]), 0.0, 5.0, 1.0)
        self.assertAlmostEqual(result[0], 0.0)

    def test_array_output(self):
        x = np.linspace(0, 10, 100)
        result = gaussian(x, 100.0, 5.0, 1.0)
        self.assertEqual(len(result), 100)


class TestCalculateMoments(unittest.TestCase):
    """Test moment statistics calculation."""

    def test_symmetric_peak(self):
        coords = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
        intensities = np.array([10.0, 50.0, 100.0, 50.0, 10.0])
        m = calculate_moments(coords, intensities)
        self.assertAlmostEqual(m["mean"], 3.0, places=1)
        self.assertAlmostEqual(m["apex"], 3.0)
        self.assertGreater(m["total_intensity"], 0)

    def test_single_point(self):
        m = calculate_moments(np.array([5.0]), np.array([100.0]))
        self.assertAlmostEqual(m["mean"], 5.0)
        self.assertAlmostEqual(m["apex"], 5.0)
        self.assertAlmostEqual(m["total_intensity"], 100.0)

    def test_empty_arrays(self):
        m = calculate_moments(np.array([]), np.array([]))
        self.assertEqual(m["mean"], 0.0)
        self.assertEqual(m["total_intensity"], 0.0)

    def test_zero_intensities(self):
        m = calculate_moments(np.array([1.0, 2.0]), np.array([0.0, 0.0]))
        self.assertEqual(m["mean"], 0.0)

    def test_skewness_symmetric(self):
        """Symmetric peak should have near-zero skewness."""
        coords = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
        intensities = np.array([10.0, 50.0, 100.0, 50.0, 10.0])
        m = calculate_moments(coords, intensities)
        self.assertAlmostEqual(m["skewness"], 0.0, places=1)

    def test_fwhm_reasonable(self):
        coords = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
        intensities = np.array([10.0, 50.0, 100.0, 50.0, 10.0])
        m = calculate_moments(coords, intensities)
        self.assertGreater(m["fwhm"], 0)
        self.assertLess(m["fwhm"], 5.0)  # Can't be wider than full range


class TestFitGaussian(unittest.TestCase):
    """Test Gaussian curve fitting."""

    def test_perfect_gaussian(self):
        """Fit to perfect Gaussian should give R2 ~1.0."""
        x = np.linspace(0, 10, 50)
        y = gaussian(x, 100.0, 5.0, 1.0)
        result = fit_gaussian(x, y)
        self.assertAlmostEqual(result["apex"], 5.0, places=1)
        self.assertAlmostEqual(result["sigma"], 1.0, places=1)
        self.assertGreater(result["r2"], 0.99)

    def test_noisy_gaussian(self):
        """Fit with moderate noise should still give good R2."""
        np.random.seed(42)
        x = np.linspace(0, 10, 50)
        y = gaussian(x, 100.0, 5.0, 1.0) + np.random.normal(0, 3, len(x))
        y = np.maximum(y, 0)
        result = fit_gaussian(x, y)
        self.assertAlmostEqual(result["apex"], 5.0, delta=0.5)
        self.assertGreater(result["r2"], 0.8)

    def test_too_few_points(self):
        result = fit_gaussian(np.array([1.0, 2.0]), np.array([50.0, 100.0]))
        self.assertEqual(result["r2"], 0.0)

    def test_empty_arrays(self):
        result = fit_gaussian(np.array([]), np.array([]))
        self.assertEqual(result["r2"], 0.0)

    def test_zero_intensities(self):
        result = fit_gaussian(np.array([1.0, 2.0, 3.0]), np.array([0.0, 0.0, 0.0]))
        self.assertEqual(result["r2"], 0.0)

    def test_r2_clamped(self):
        """R2 should always be in [0, 1]."""
        x = np.linspace(0, 10, 50)
        y = gaussian(x, 100.0, 5.0, 1.0)
        result = fit_gaussian(x, y)
        self.assertGreaterEqual(result["r2"], 0.0)
        self.assertLessEqual(result["r2"], 1.0)


class TestSimpleIsotopeScore(unittest.TestCase):
    """Test simple isotope pattern scoring."""

    def test_good_pattern(self):
        """Typical peptide isotope pattern should score well."""
        observed = [0.8, 1.0, 0.6, 0.3, 0.1]
        score = _simple_isotope_score(observed)
        self.assertGreater(score, 0.5)

    def test_single_peak(self):
        score = _simple_isotope_score([1.0])
        self.assertEqual(score, 0.0)

    def test_empty(self):
        score = _simple_isotope_score([])
        self.assertEqual(score, 0.0)

    def test_all_zeros(self):
        score = _simple_isotope_score([0.0, 0.0, 0.0])
        self.assertEqual(score, 0.0)

    def test_decreasing_pattern(self):
        """Monotonically decreasing pattern should score well."""
        observed = [1.0, 0.7, 0.4, 0.2, 0.1]
        score = _simple_isotope_score(observed)
        self.assertGreater(score, 0.5)


class TestComputeQualitySummary(unittest.TestCase):
    """Test quality summary computation."""

    def test_high_quality(self):
        result = compute_quality_summary(rt_r2=0.95, im_r2=0.90, isotope_cosim=0.95)
        self.assertTrue(result["is_high_quality"])
        self.assertEqual(result["rt_quality"], "good")
        self.assertEqual(result["im_quality"], "good")
        self.assertEqual(result["isotope_quality"], "good")

    def test_low_quality(self):
        result = compute_quality_summary(rt_r2=0.5, im_r2=0.3, isotope_cosim=0.4)
        self.assertFalse(result["is_high_quality"])
        self.assertEqual(result["rt_quality"], "poor")

    def test_mixed_quality(self):
        result = compute_quality_summary(rt_r2=0.95, im_r2=0.5, isotope_cosim=0.95)
        self.assertFalse(result["is_high_quality"])
        self.assertEqual(result["rt_quality"], "good")
        self.assertEqual(result["im_quality"], "poor")

    def test_custom_thresholds(self):
        thresholds = {"rt_r2_min": 0.5, "im_r2_min": 0.5, "isotope_cosim_min": 0.5}
        result = compute_quality_summary(
            rt_r2=0.6, im_r2=0.6, isotope_cosim=0.6, thresholds=thresholds
        )
        self.assertTrue(result["is_high_quality"])

    def test_combined_score_range(self):
        result = compute_quality_summary(rt_r2=1.0, im_r2=1.0, isotope_cosim=1.0)
        self.assertAlmostEqual(result["combined_score"], 1.0)

        result = compute_quality_summary(rt_r2=0.0, im_r2=0.0, isotope_cosim=0.0)
        self.assertAlmostEqual(result["combined_score"], 0.0)

    def test_boundary_values(self):
        """Exactly at threshold is 'good'."""
        result = compute_quality_summary(rt_r2=0.8, im_r2=0.8, isotope_cosim=0.9)
        self.assertTrue(result["is_high_quality"])


if __name__ == "__main__":
    unittest.main(verbosity=2)
