#!/usr/bin/env python3
"""Tests for sequence matcher, coordinate matcher, and consensus calculation."""

import unittest
import pandas as pd
import numpy as np

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

from precursor_merging.config import MatchConfig, MatchTier, MatchResult
from precursor_merging.sequence_matcher import SequenceMatcher, normalize_sequence_for_matching
from precursor_merging.coordinate_matcher import CoordinateMatcher
from precursor_merging.consensus import count_engines, get_consensus_peptide, calculate_consensus


# =============================================================================
# SequenceMatcher tests
# =============================================================================


class TestNormalizeSequenceForMatching(unittest.TestCase):
    """Test sequence normalization for matching."""

    def test_plain_sequence(self):
        self.assertEqual(normalize_sequence_for_matching("PEPTIDEK"), "PEPTLDEK")

    def test_with_unimod(self):
        self.assertEqual(normalize_sequence_for_matching("[]-PEPTIDEK-[]"), "PEPTLDEK")

    def test_il_normalization(self):
        self.assertEqual(
            normalize_sequence_for_matching("PEPTIDEK", normalize_il=True),
            "PEPTLDEK"
        )
        self.assertEqual(
            normalize_sequence_for_matching("PEPTIDEK", normalize_il=False),
            "PEPTIDEK"
        )

    def test_empty_input(self):
        self.assertEqual(normalize_sequence_for_matching(""), "")
        self.assertEqual(normalize_sequence_for_matching(None), "")

    def test_nan_input(self):
        self.assertEqual(normalize_sequence_for_matching(float("nan")), "")


class TestSequenceMatcher(unittest.TestCase):
    """Test sequence-based matching."""

    def setUp(self):
        self.config = MatchConfig()
        self.matcher = SequenceMatcher(self.config)

    def test_exact_match(self):
        source = pd.DataFrame({
            "seq": ["PEPTIDEK"], "charge": [2], "raw_file": ["run1"]
        })
        target = pd.DataFrame({
            "seq": ["PEPTIDEK"], "charge": [2], "raw_file": ["run1"]
        })
        matches = self.matcher.match(
            source, target, "seq", "charge", "seq", "charge"
        )
        self.assertEqual(len(matches), 1)
        self.assertEqual(matches[0][2], MatchTier.SEQUENCE_IL_NORM)

    def test_il_normalized_match(self):
        """I and L should match."""
        source = pd.DataFrame({
            "seq": ["PEPTIDEK"], "charge": [2], "raw_file": ["run1"]
        })
        target = pd.DataFrame({
            "seq": ["PEPTLDEK"], "charge": [2], "raw_file": ["run1"]
        })
        matches = self.matcher.match(
            source, target, "seq", "charge", "seq", "charge"
        )
        self.assertEqual(len(matches), 1)

    def test_charge_mismatch(self):
        source = pd.DataFrame({
            "seq": ["PEPTIDEK"], "charge": [2], "raw_file": ["run1"]
        })
        target = pd.DataFrame({
            "seq": ["PEPTIDEK"], "charge": [3], "raw_file": ["run1"]
        })
        matches = self.matcher.match(
            source, target, "seq", "charge", "seq", "charge"
        )
        self.assertEqual(len(matches), 0)

    def test_raw_file_mismatch(self):
        source = pd.DataFrame({
            "seq": ["PEPTIDEK"], "charge": [2], "raw_file": ["run1"]
        })
        target = pd.DataFrame({
            "seq": ["PEPTIDEK"], "charge": [2], "raw_file": ["run2"]
        })
        matches = self.matcher.match(
            source, target, "seq", "charge", "seq", "charge"
        )
        self.assertEqual(len(matches), 0)

    def test_multiple_matches(self):
        source = pd.DataFrame({
            "seq": ["PEPTIDEK", "ACDEFGH"], "charge": [2, 3], "raw_file": ["run1", "run1"]
        })
        target = pd.DataFrame({
            "seq": ["PEPTIDEK", "ACDEFGH"], "charge": [2, 3], "raw_file": ["run1", "run1"]
        })
        matches = self.matcher.match(
            source, target, "seq", "charge", "seq", "charge"
        )
        self.assertEqual(len(matches), 2)

    def test_no_double_matching(self):
        """Each target should match at most once."""
        source = pd.DataFrame({
            "seq": ["PEPTIDEK", "PEPTIDEK"], "charge": [2, 2], "raw_file": ["run1", "run1"]
        })
        target = pd.DataFrame({
            "seq": ["PEPTIDEK"], "charge": [2], "raw_file": ["run1"]
        })
        matches = self.matcher.match(
            source, target, "seq", "charge", "seq", "charge"
        )
        self.assertEqual(len(matches), 1)

    def test_nan_handling(self):
        source = pd.DataFrame({
            "seq": [None, "PEPTIDEK"], "charge": [2, float("nan")], "raw_file": ["run1", "run1"]
        })
        target = pd.DataFrame({
            "seq": ["PEPTIDEK"], "charge": [2], "raw_file": ["run1"]
        })
        matches = self.matcher.match(
            source, target, "seq", "charge", "seq", "charge"
        )
        self.assertEqual(len(matches), 0)

    def test_create_index(self):
        df = pd.DataFrame({
            "seq": ["PEPTIDEK", "ACDEFGH"], "charge": [2, 3], "raw_file": ["run1", "run1"]
        })
        index = self.matcher.create_index(df, "seq", "charge")
        self.assertEqual(len(index), 2)

    def test_match_to_index(self):
        target = pd.DataFrame({
            "seq": ["PEPTIDEK"], "charge": [2], "raw_file": ["run1"]
        })
        index = self.matcher.create_index(target, "seq", "charge")

        source_row = pd.Series({"seq": "PEPTIDEK", "charge": 2, "raw_file": "run1"})
        result = self.matcher.match_to_index(source_row, index, "seq", "charge")
        self.assertIsNotNone(result)

        miss_row = pd.Series({"seq": "XXXXXXX", "charge": 2, "raw_file": "run1"})
        result = self.matcher.match_to_index(miss_row, index, "seq", "charge")
        self.assertIsNone(result)


# =============================================================================
# CoordinateMatcher tests
# =============================================================================


class TestCoordinateMatcher(unittest.TestCase):
    """Test coordinate-based matching."""

    def setUp(self):
        self.config = MatchConfig(mz_tol_ppm=20.0, rt_tol_sec=5.0, im_tol=0.05)
        self.matcher = CoordinateMatcher(self.config)

    def test_exact_coordinate_match(self):
        source = pd.DataFrame({
            "mz": [500.25], "charge": [2], "rt": [120.0], "im": [0.95], "raw_file": ["run1"]
        })
        target = pd.DataFrame({
            "mz": [500.25], "charge": [2], "rt": [120.0], "im": [0.95], "raw_file": ["run1"]
        })
        matches = self.matcher.match(
            source, target, "mz", "charge", "rt", "im", "mz", "charge", "rt", "im"
        )
        self.assertEqual(len(matches), 1)
        _, _, result = matches[0]
        self.assertTrue(result.matched)
        self.assertEqual(result.tier, MatchTier.COORDINATE_FULL)

    def test_within_ppm_tolerance(self):
        """m/z within 20 ppm should match."""
        mz = 500.0
        mz_shifted = mz + mz * 15e-6  # 15 ppm shift
        source = pd.DataFrame({
            "mz": [mz], "charge": [2], "rt": [120.0], "im": [0.95], "raw_file": ["run1"]
        })
        target = pd.DataFrame({
            "mz": [mz_shifted], "charge": [2], "rt": [120.0], "im": [0.95], "raw_file": ["run1"]
        })
        matches = self.matcher.match(
            source, target, "mz", "charge", "rt", "im", "mz", "charge", "rt", "im"
        )
        self.assertEqual(len(matches), 1)

    def test_outside_ppm_tolerance(self):
        """m/z outside 20 ppm should not match."""
        mz = 500.0
        mz_shifted = mz + mz * 25e-6  # 25 ppm shift
        source = pd.DataFrame({
            "mz": [mz], "charge": [2], "rt": [120.0], "im": [0.95], "raw_file": ["run1"]
        })
        target = pd.DataFrame({
            "mz": [mz_shifted], "charge": [2], "rt": [120.0], "im": [0.95], "raw_file": ["run1"]
        })
        matches = self.matcher.match(
            source, target, "mz", "charge", "rt", "im", "mz", "charge", "rt", "im"
        )
        self.assertEqual(len(matches), 0)

    def test_rt_tolerance(self):
        """RT within tolerance matches, outside doesn't."""
        source = pd.DataFrame({
            "mz": [500.0], "charge": [2], "rt": [120.0], "im": [0.95], "raw_file": ["run1"]
        })
        # Within 5s RT tolerance
        target_ok = pd.DataFrame({
            "mz": [500.0], "charge": [2], "rt": [124.0], "im": [0.95], "raw_file": ["run1"]
        })
        matches = self.matcher.match(
            source, target_ok, "mz", "charge", "rt", "im", "mz", "charge", "rt", "im"
        )
        self.assertEqual(len(matches), 1)

        # Outside 5s RT tolerance
        target_bad = pd.DataFrame({
            "mz": [500.0], "charge": [2], "rt": [126.0], "im": [0.95], "raw_file": ["run1"]
        })
        matches = self.matcher.match(
            source, target_bad, "mz", "charge", "rt", "im", "mz", "charge", "rt", "im"
        )
        self.assertEqual(len(matches), 0)

    def test_im_tolerance(self):
        """IM within tolerance matches, outside doesn't."""
        source = pd.DataFrame({
            "mz": [500.0], "charge": [2], "rt": [120.0], "im": [0.95], "raw_file": ["run1"]
        })
        # Within 0.05 IM tolerance
        target_ok = pd.DataFrame({
            "mz": [500.0], "charge": [2], "rt": [120.0], "im": [0.99], "raw_file": ["run1"]
        })
        matches = self.matcher.match(
            source, target_ok, "mz", "charge", "rt", "im", "mz", "charge", "rt", "im"
        )
        self.assertEqual(len(matches), 1)

        # Outside 0.05 IM tolerance
        target_bad = pd.DataFrame({
            "mz": [500.0], "charge": [2], "rt": [120.0], "im": [1.01], "raw_file": ["run1"]
        })
        matches = self.matcher.match(
            source, target_bad, "mz", "charge", "rt", "im", "mz", "charge", "rt", "im"
        )
        self.assertEqual(len(matches), 0)

    def test_charge_must_match(self):
        source = pd.DataFrame({
            "mz": [500.0], "charge": [2], "rt": [120.0], "im": [0.95], "raw_file": ["run1"]
        })
        target = pd.DataFrame({
            "mz": [500.0], "charge": [3], "rt": [120.0], "im": [0.95], "raw_file": ["run1"]
        })
        matches = self.matcher.match(
            source, target, "mz", "charge", "rt", "im", "mz", "charge", "rt", "im"
        )
        self.assertEqual(len(matches), 0)

    def test_best_match_selected(self):
        """When multiple candidates, closest should be selected."""
        source = pd.DataFrame({
            "mz": [500.0], "charge": [2], "rt": [120.0], "im": [0.95], "raw_file": ["run1"]
        })
        target = pd.DataFrame({
            "mz": [500.005, 500.001], "charge": [2, 2],
            "rt": [120.0, 120.0], "im": [0.95, 0.95], "raw_file": ["run1", "run1"]
        })
        matches = self.matcher.match(
            source, target, "mz", "charge", "rt", "im", "mz", "charge", "rt", "im"
        )
        self.assertEqual(len(matches), 1)
        # Should match index 1 (500.001 is closer to 500.0)
        self.assertEqual(matches[0][1], 1)

    def test_nan_coordinates_skipped(self):
        source = pd.DataFrame({
            "mz": [float("nan")], "charge": [2], "rt": [120.0], "im": [0.95], "raw_file": ["run1"]
        })
        target = pd.DataFrame({
            "mz": [500.0], "charge": [2], "rt": [120.0], "im": [0.95], "raw_file": ["run1"]
        })
        matches = self.matcher.match(
            source, target, "mz", "charge", "rt", "im", "mz", "charge", "rt", "im"
        )
        self.assertEqual(len(matches), 0)

    def test_match_result_diagnostics(self):
        """MatchResult should contain diagnostic values."""
        source = pd.DataFrame({
            "mz": [500.0], "charge": [2], "rt": [120.0], "im": [0.95], "raw_file": ["run1"]
        })
        target = pd.DataFrame({
            "mz": [500.005], "charge": [2], "rt": [121.0], "im": [0.96], "raw_file": ["run1"]
        })
        matches = self.matcher.match(
            source, target, "mz", "charge", "rt", "im", "mz", "charge", "rt", "im"
        )
        self.assertEqual(len(matches), 1)
        _, _, result = matches[0]
        self.assertIsNotNone(result.mz_diff_ppm)
        self.assertIsNotNone(result.rt_diff_sec)
        self.assertIsNotNone(result.im_diff)
        self.assertAlmostEqual(result.rt_diff_sec, 1.0, places=1)
        self.assertAlmostEqual(result.im_diff, 0.01, places=2)


# =============================================================================
# Consensus tests
# =============================================================================


class TestCountEngines(unittest.TestCase):
    """Test engine counting."""

    def test_all_three(self):
        row = pd.Series({
            "fragpipe_peptide": "PEP", "diann_peptide": "PEP", "sage_peptide": "PEP"
        })
        self.assertEqual(count_engines(row), 3)

    def test_two_engines(self):
        row = pd.Series({
            "fragpipe_peptide": "PEP", "diann_peptide": None, "sage_peptide": "PEP"
        })
        self.assertEqual(count_engines(row), 2)

    def test_one_engine(self):
        row = pd.Series({
            "fragpipe_peptide": "PEP", "diann_peptide": None, "sage_peptide": None
        })
        self.assertEqual(count_engines(row), 1)

    def test_zero_engines(self):
        row = pd.Series({
            "fragpipe_peptide": None, "diann_peptide": None, "sage_peptide": None
        })
        self.assertEqual(count_engines(row), 0)

    def test_empty_string_not_counted(self):
        row = pd.Series({
            "fragpipe_peptide": "", "diann_peptide": None, "sage_peptide": "PEP"
        })
        self.assertEqual(count_engines(row), 1)


class TestGetConsensusPeptide(unittest.TestCase):
    """Test consensus peptide selection."""

    def test_all_agree(self):
        row = pd.Series({
            "fragpipe_peptide": "PEPTIDEK", "fragpipe_probability": 0.99,
            "diann_peptide": "PEPTIDEK", "diann_qvalue": 0.01,
            "sage_peptide": "PEPTIDEK", "sage_qvalue": 0.01,
        })
        peptide, weight = get_consensus_peptide(row)
        self.assertEqual(weight, 1.0)
        self.assertNotEqual(peptide, "")

    def test_disagreement_picks_best(self):
        row = pd.Series({
            "fragpipe_peptide": "PEPTIDEK", "fragpipe_probability": 0.99,
            "diann_peptide": "DIFFERENTK", "diann_qvalue": 0.5,
            "sage_peptide": None, "sage_qvalue": None,
        })
        peptide, weight = get_consensus_peptide(row)
        self.assertLess(weight, 1.0)
        self.assertGreater(weight, 0.0)

    def test_no_engines(self):
        row = pd.Series({
            "fragpipe_peptide": None, "diann_peptide": None, "sage_peptide": None,
        })
        peptide, weight = get_consensus_peptide(row)
        self.assertEqual(peptide, "")
        self.assertEqual(weight, 0.0)

    def test_il_normalization_in_consensus(self):
        """I and L should be treated as same for consensus."""
        row = pd.Series({
            "fragpipe_peptide": "PEPTIDEK", "fragpipe_probability": 0.99,
            "diann_peptide": "PEPTLDEK", "diann_qvalue": 0.01,
            "sage_peptide": None, "sage_qvalue": None,
        })
        peptide, weight = get_consensus_peptide(row)
        # I/L normalized -> same sequence -> agreement
        self.assertEqual(weight, 1.0)


class TestCalculateConsensus(unittest.TestCase):
    """Test full consensus calculation on DataFrame."""

    def test_adds_columns(self):
        df = pd.DataFrame({
            "fragpipe_peptide": ["PEP", None],
            "diann_peptide": ["PEP", "PEP"],
            "sage_peptide": ["PEP", None],
            "fragpipe_probability": [0.99, None],
            "diann_qvalue": [0.01, 0.01],
            "sage_qvalue": [0.01, None],
        })
        result = calculate_consensus(df)
        self.assertIn("n_engines", result.columns)
        self.assertIn("consensus_peptide", result.columns)
        self.assertIn("confidence_weight", result.columns)
        self.assertEqual(result.iloc[0]["n_engines"], 3)
        self.assertEqual(result.iloc[1]["n_engines"], 1)


# =============================================================================
# MatchConfig / MatchResult tests
# =============================================================================


class TestMatchConfig(unittest.TestCase):
    """Test MatchConfig defaults."""

    def test_defaults(self):
        config = MatchConfig()
        self.assertEqual(config.mz_tol_ppm, 20.0)
        self.assertEqual(config.rt_tol_sec, 0.01)
        self.assertEqual(config.im_tol, 0.05)

    def test_custom_values(self):
        config = MatchConfig(mz_tol_ppm=10.0, rt_tol_sec=1.0)
        self.assertEqual(config.mz_tol_ppm, 10.0)
        self.assertEqual(config.rt_tol_sec, 1.0)


class TestMatchResult(unittest.TestCase):
    """Test MatchResult factory methods."""

    def test_no_match(self):
        r = MatchResult.no_match()
        self.assertFalse(r.matched)
        self.assertEqual(r.tier, MatchTier.NO_MATCH)

    def test_sequence_match(self):
        r = MatchResult.sequence_match(normalize_il=True)
        self.assertTrue(r.matched)
        self.assertEqual(r.tier, MatchTier.SEQUENCE_IL_NORM)
        self.assertTrue(r.sequence_match)

    def test_coordinate_match_full(self):
        r = MatchResult.coordinate_match(score=0.5, full=True, mz_diff_ppm=5.0)
        self.assertTrue(r.matched)
        self.assertEqual(r.tier, MatchTier.COORDINATE_FULL)

    def test_coordinate_match_partial(self):
        r = MatchResult.coordinate_match(score=0.5, full=False, mz_diff_ppm=5.0)
        self.assertEqual(r.tier, MatchTier.COORDINATE_PARTIAL)


if __name__ == "__main__":
    unittest.main(verbosity=2)
