#!/usr/bin/env python3
"""
Tests for dashboard backend utility functions.

Tests pure functions extracted from main.py (no server startup needed):
- sage_to_observed_mz: Sage deconvolution reversal
- standardize_modified_sequence: Multi-engine mod normalization
- normalize_peptide_sequence: Strip mods to plain AA
- compute_spectrum_cosine: Spectrum similarity
- one_over_k0_to_ccs: Ion mobility -> CCS conversion
- _clean_str / _to_float: Robust type coercion
"""

import unittest
import numpy as np
import pandas as pd

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))

from main import (
    sage_to_observed_mz,
    standardize_modified_sequence,
    normalize_peptide_sequence,
    compute_spectrum_cosine,
    one_over_k0_to_ccs,
    _clean_str,
    _to_float,
    PROTON_MASS,
)


# =============================================================================
# sage_to_observed_mz
# =============================================================================

class TestSageToObservedMz(unittest.TestCase):
    """Test Sage deconvolved m/z -> observed m/z conversion."""

    def test_below_max_unchanged(self):
        """m/z below max passes through unchanged."""
        mz, z = sage_to_observed_mz(500.0, 1)
        self.assertAlmostEqual(mz, 500.0)
        self.assertEqual(z, 1)

    def test_at_boundary(self):
        """m/z exactly at max passes through."""
        mz, z = sage_to_observed_mz(1750.0, 1)
        self.assertAlmostEqual(mz, 1750.0)
        self.assertEqual(z, 1)

    def test_above_max_reconverted(self):
        """High m/z (deconvolved) gets reconverted to lower charge state."""
        # Sage reports neutral_mass + proton for a z=2 fragment
        # neutral_mass = 2000, sage_mz = 2000 + 1.007276 = 2001.007276
        sage_mz = 2001.007276
        mz, z = sage_to_observed_mz(sage_mz, 1)
        # Should reconvert: (2000 + 2*1.007276) / 2 = 1001.007276
        self.assertAlmostEqual(mz, 1001.007276, places=3)
        self.assertEqual(z, 2)

    def test_triply_charged(self):
        """Large neutral mass reconverted to z=3."""
        neutral = 4000.0
        sage_mz = neutral + PROTON_MASS  # Sage reports at z=1
        mz, z = sage_to_observed_mz(sage_mz, 1)
        expected = (neutral + 3 * PROTON_MASS) / 3
        if expected <= 1750.0:
            self.assertAlmostEqual(mz, expected, places=3)
            self.assertEqual(z, 3)

    def test_custom_max(self):
        """Custom mz_max works."""
        mz, z = sage_to_observed_mz(600.0, 1, mz_max=500.0)
        # 600 > 500 -> should reconvert
        self.assertLessEqual(mz, 500.0)

    def test_no_valid_charge_returns_original(self):
        """If no z in [2,4] gives obs <= max, return original."""
        # Very large mass that can't be brought down
        mz, z = sage_to_observed_mz(100000.0, 1, mz_max=1750.0)
        self.assertAlmostEqual(mz, 100000.0)


# =============================================================================
# standardize_modified_sequence
# =============================================================================

class TestStandardizeModifiedSequence(unittest.TestCase):
    """Test multi-format modification standardization."""

    def test_sage_format(self):
        result = standardize_modified_sequence("C[+57.021465]PEPTIDEK")
        self.assertIn("[UNIMOD:4]", result)

    def test_diann_format(self):
        result = standardize_modified_sequence("PEPTIDE(UniMod:35)K")
        self.assertIn("[UNIMOD:35]", result)

    def test_fragpipe_absolute_mass(self):
        result = standardize_modified_sequence("M[147.0354]PEPTIDEK")
        self.assertIn("[UNIMOD:35]", result)

    def test_unmodified(self):
        result = standardize_modified_sequence("PEPTIDEK")
        self.assertEqual(result, "PEPTIDEK")

    def test_none_input(self):
        result = standardize_modified_sequence(None)
        self.assertIsNone(result)

    def test_nan_input(self):
        result = standardize_modified_sequence(float("nan"))
        # Should handle gracefully
        self.assertTrue(result is None or pd.isna(result))

    def test_nterm_sage(self):
        result = standardize_modified_sequence("[+42.010567]-MPEPTIDEK")
        self.assertIn("[UNIMOD:1]", result)

    def test_unknown_mass_preserved(self):
        result = standardize_modified_sequence("PEPTIDE[+999.999]K")
        self.assertIn("[+999.999]", result)


# =============================================================================
# normalize_peptide_sequence
# =============================================================================

class TestNormalizePeptideSequence(unittest.TestCase):
    """Test stripping mods to bare amino acid sequence."""

    def test_plain_sequence(self):
        self.assertEqual(normalize_peptide_sequence("PEPTIDEK"), "PEPTIDEK")

    def test_unimod_stripped(self):
        self.assertEqual(normalize_peptide_sequence("C[UNIMOD:4]PEPTIDEK"), "CPEPTIDEK")

    def test_sage_format_stripped(self):
        self.assertEqual(normalize_peptide_sequence("C[+57.021]PEPTIDEK"), "CPEPTIDEK")

    def test_diann_format_stripped(self):
        self.assertEqual(normalize_peptide_sequence("C(UniMod:4)PEPTIDEK"), "CPEPTIDEK")

    def test_terminal_format_stripped(self):
        self.assertEqual(normalize_peptide_sequence("[UNIMOD:1]-MPEPTIDEK-[]"), "MPEPTIDEK")

    def test_empty_input(self):
        self.assertEqual(normalize_peptide_sequence(""), "")

    def test_nterm_stripped(self):
        self.assertEqual(normalize_peptide_sequence("[+42.0106]MPEPTIDEK"), "MPEPTIDEK")


# =============================================================================
# compute_spectrum_cosine
# =============================================================================

class TestComputeSpectrumCosine(unittest.TestCase):
    """Test spectrum cosine similarity."""

    def test_perfect_match(self):
        obs_mz = [100.0, 200.0, 300.0]
        obs_int = [50.0, 100.0, 75.0]
        frags = [
            {"mz_observed": 100.0, "intensity": 50.0},
            {"mz_observed": 200.0, "intensity": 100.0},
            {"mz_observed": 300.0, "intensity": 75.0},
        ]
        cos = compute_spectrum_cosine(obs_mz, obs_int, frags)
        self.assertAlmostEqual(cos, 1.0, places=3)

    def test_orthogonal(self):
        """No matching peaks -> None."""
        obs_mz = [100.0, 200.0]
        obs_int = [50.0, 100.0]
        frags = [{"mz_observed": 999.0, "intensity": 100.0}]
        cos = compute_spectrum_cosine(obs_mz, obs_int, frags)
        self.assertIsNone(cos)  # <2 matches

    def test_empty_inputs(self):
        self.assertIsNone(compute_spectrum_cosine([], [], []))
        self.assertIsNone(compute_spectrum_cosine([100.0], [50.0], []))

    def test_within_ppm_tolerance(self):
        """Peaks within 20ppm should match."""
        mz = 500.0
        shifted = mz + mz * 10e-6  # 10 ppm shift
        obs_mz = [mz, 600.0]
        obs_int = [100.0, 50.0]
        frags = [
            {"mz_observed": shifted, "intensity": 100.0},
            {"mz_observed": 600.0, "intensity": 50.0},
        ]
        cos = compute_spectrum_cosine(obs_mz, obs_int, frags)
        self.assertIsNotNone(cos)
        self.assertGreater(cos, 0.9)

    def test_cosine_range(self):
        """Result should be in [0, 1]."""
        obs_mz = [100.0, 200.0, 300.0]
        obs_int = [50.0, 100.0, 75.0]
        frags = [
            {"mz_observed": 100.0, "intensity": 10.0},
            {"mz_observed": 200.0, "intensity": 200.0},
            {"mz_observed": 300.0, "intensity": 5.0},
        ]
        cos = compute_spectrum_cosine(obs_mz, obs_int, frags)
        self.assertGreaterEqual(cos, 0.0)
        self.assertLessEqual(cos, 1.0)


# =============================================================================
# one_over_k0_to_ccs
# =============================================================================

class TestOneOverK0ToCCS(unittest.TestCase):
    """Test 1/K0 -> CCS conversion (Mason-Schamp)."""

    def test_reasonable_ccs_range(self):
        """Typical peptide CCS should be 200-600 A^2."""
        ccs = one_over_k0_to_ccs(one_over_k0=1.0, mz=500.0, charge=2)
        self.assertGreater(ccs, 100)
        self.assertLess(ccs, 1000)

    def test_higher_mobility_higher_ccs(self):
        """Higher 1/K0 (lower mobility) -> higher CCS."""
        ccs_low = one_over_k0_to_ccs(0.8, 500.0, 2)
        ccs_high = one_over_k0_to_ccs(1.2, 500.0, 2)
        self.assertGreater(ccs_high, ccs_low)

    def test_higher_charge_higher_ccs(self):
        """Higher charge -> higher CCS (same 1/K0)."""
        ccs_z2 = one_over_k0_to_ccs(1.0, 500.0, 2)
        ccs_z3 = one_over_k0_to_ccs(1.0, 500.0, 3)
        self.assertGreater(ccs_z3, ccs_z2)

    def test_zero_mobility_returns_zero(self):
        self.assertEqual(one_over_k0_to_ccs(0.0, 500.0, 2), 0.0)

    def test_zero_charge_returns_zero(self):
        self.assertEqual(one_over_k0_to_ccs(1.0, 500.0, 0), 0.0)

    def test_negative_mobility_returns_zero(self):
        self.assertEqual(one_over_k0_to_ccs(-1.0, 500.0, 2), 0.0)


# =============================================================================
# _clean_str / _to_float
# =============================================================================

class TestCleanStr(unittest.TestCase):
    """Test robust string cleaning."""

    def test_normal_string(self):
        self.assertEqual(_clean_str("hello"), "hello")

    def test_whitespace_stripped(self):
        self.assertEqual(_clean_str("  hello  "), "hello")

    def test_none_returns_none(self):
        self.assertIsNone(_clean_str(None))

    def test_nan_returns_none(self):
        self.assertIsNone(_clean_str(float("nan")))

    def test_empty_string_returns_none(self):
        self.assertIsNone(_clean_str(""))

    def test_nan_string_returns_none(self):
        self.assertIsNone(_clean_str("nan"))
        self.assertIsNone(_clean_str("NaN"))

    def test_numeric_converted(self):
        self.assertEqual(_clean_str(42), "42")


class TestToFloat(unittest.TestCase):
    """Test robust float conversion."""

    def test_normal_float(self):
        self.assertAlmostEqual(_to_float(3.14), 3.14)

    def test_int_to_float(self):
        self.assertAlmostEqual(_to_float(42), 42.0)

    def test_string_to_float(self):
        self.assertAlmostEqual(_to_float("3.14"), 3.14)

    def test_none_returns_none(self):
        self.assertIsNone(_to_float(None))

    def test_nan_returns_none(self):
        self.assertIsNone(_to_float(float("nan")))

    def test_invalid_string_returns_none(self):
        self.assertIsNone(_to_float("not_a_number"))

    def test_numpy_nan(self):
        self.assertIsNone(_to_float(np.nan))


if __name__ == "__main__":
    unittest.main(verbosity=2)
