#!/usr/bin/env python3
"""
Tests for step5_merge.py - merge logic, consensus, quality columns.

Tests internal functions:
- _merge_datasets: Join precursor_index with raw_features
- _add_consensus_columns: n_engines, consensus_peptide, confidence_weight
- _add_quality_columns: is_high_quality flag
- _compute_quality_summary: Engine distribution + quality stats
"""

import unittest
import pandas as pd
import numpy as np

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from runner.steps.step5_merge import (
    _merge_datasets,
    _add_consensus_columns,
    _add_quality_columns,
    _compute_quality_summary,
)


def _make_index(n=3):
    """Create synthetic precursor_index."""
    return pd.DataFrame({
        "precursor_id": range(1, n + 1),
        "raw_file": ["run1"] * n,
        "mz": [500.0 + i * 10 for i in range(n)],
        "charge": [2] * n,
        "fragpipe_peptide": ["AAAK", "BBBK", None],
        "diann_peptide": ["AAAK", None, "CCCK"],
        "sage_peptide": ["AAAK", "BBBK", "CCCK"],
        "fragpipe_modified": ["[]-AAAK-[]", "[]-BBBK-[]", None],
        "diann_modified": ["[]-AAAK-[]", None, "[]-CCCK-[]"],
        "sage_modified": ["[]-AAAK-[]", "[]-BBBK-[]", "[]-CCCK-[]"],
        "n_engines": [3, 2, 2],
    })


def _make_features(n=3):
    """Create synthetic raw_features."""
    return pd.DataFrame({
        "precursor_id": range(1, n + 1),
        "raw_file": ["run1"] * n,
        "charge": [2] * n,
        "mz": [500.0 + i * 10 for i in range(n)],
        "rt_seconds": [100.0 + i * 50 for i in range(n)],
        "mobility": [0.9 + i * 0.02 for i in range(n)],
        "ms1_rt_r2": [0.95, 0.50, 0.85],
        "ms1_im_r2": [0.90, 0.40, 0.82],
        "isotope_cosim": [0.95, 0.60, 0.92],
    })


class TestMergeDatasets(unittest.TestCase):
    """Test merging precursor_index with raw_features."""

    def test_join_on_precursor_id(self):
        index = _make_index()
        features = _make_features()
        merged = _merge_datasets(index, features)
        self.assertEqual(len(merged), 3)
        self.assertIn("ms1_rt_r2", merged.columns)
        self.assertIn("fragpipe_peptide", merged.columns)

    def test_empty_features(self):
        """Empty raw_features -> use index only."""
        index = _make_index()
        merged = _merge_datasets(index, pd.DataFrame())
        self.assertEqual(len(merged), 3)

    def test_empty_index(self):
        """Empty index -> use features only."""
        features = _make_features()
        merged = _merge_datasets(pd.DataFrame(), features)
        self.assertEqual(len(merged), 3)
        self.assertIn("n_engines", merged.columns)

    def test_raw_file_normalization(self):
        """'.d' stripped for matching."""
        index = pd.DataFrame({
            "precursor_id": [1], "raw_file": ["run1.d"],
            "fragpipe_peptide": ["PEP"], "n_engines": [1],
        })
        features = pd.DataFrame({
            "precursor_id": [1], "raw_file": ["run1"],
            "ms1_rt_r2": [0.9],
        })
        merged = _merge_datasets(index, features)
        self.assertEqual(len(merged), 1)
        self.assertIn("ms1_rt_r2", merged.columns)

    def test_n_engines_filled(self):
        """Unmatched raw precursors get n_engines=0."""
        index = pd.DataFrame({
            "precursor_id": [1], "raw_file": ["run1"],
            "fragpipe_peptide": ["PEP"], "n_engines": [1],
        })
        features = pd.DataFrame({
            "precursor_id": [1, 2], "raw_file": ["run1", "run1"],
        })
        merged = _merge_datasets(index, features)
        # pid=2 has no engine match -> n_engines=0
        pid2 = merged[merged["precursor_id"] == 2]
        if not pid2.empty:
            self.assertEqual(pid2.iloc[0]["n_engines"], 0)


class TestAddConsensusColumns(unittest.TestCase):
    """Test consensus column generation."""

    def test_all_three_engines(self):
        df = pd.DataFrame({
            "fragpipe_peptide": ["PEP"], "fragpipe_modified": ["[]-PEP-[]"],
            "diann_peptide": ["PEP"], "diann_modified": ["[]-PEP-[]"],
            "sage_peptide": ["PEP"], "sage_modified": ["[]-PEP-[]"],
        })
        result = _add_consensus_columns(df)
        self.assertEqual(result.iloc[0]["n_engines"], 3)
        self.assertEqual(result.iloc[0]["confidence_weight"], 1.0)
        self.assertNotEqual(result.iloc[0]["consensus_peptide"], "")

    def test_two_engines(self):
        df = pd.DataFrame({
            "fragpipe_peptide": ["PEP"], "fragpipe_modified": ["[]-PEP-[]"],
            "diann_peptide": [None], "diann_modified": [None],
            "sage_peptide": ["PEP"], "sage_modified": ["[]-PEP-[]"],
        })
        result = _add_consensus_columns(df)
        self.assertEqual(result.iloc[0]["n_engines"], 2)
        self.assertEqual(result.iloc[0]["confidence_weight"], 0.9)

    def test_one_engine(self):
        df = pd.DataFrame({
            "fragpipe_peptide": ["PEP"], "fragpipe_modified": ["[]-PEP-[]"],
            "diann_peptide": [None], "diann_modified": [None],
            "sage_peptide": [None], "sage_modified": [None],
        })
        result = _add_consensus_columns(df)
        self.assertEqual(result.iloc[0]["n_engines"], 1)
        self.assertEqual(result.iloc[0]["confidence_weight"], 0.7)

    def test_zero_engines(self):
        df = pd.DataFrame({
            "fragpipe_peptide": [None], "fragpipe_modified": [None],
            "diann_peptide": [None], "diann_modified": [None],
            "sage_peptide": [None], "sage_modified": [None],
        })
        result = _add_consensus_columns(df)
        self.assertEqual(result.iloc[0]["n_engines"], 0)
        self.assertEqual(result.iloc[0]["confidence_weight"], 0.0)
        self.assertEqual(result.iloc[0]["consensus_peptide"], "")

    def test_consensus_prefers_fragpipe(self):
        """When all available, FragPipe modified seq is used first."""
        df = pd.DataFrame({
            "fragpipe_peptide": ["PEP"], "fragpipe_modified": ["FP_MOD"],
            "diann_peptide": ["PEP"], "diann_modified": ["DN_MOD"],
            "sage_peptide": ["PEP"], "sage_modified": ["SG_MOD"],
        })
        result = _add_consensus_columns(df)
        self.assertEqual(result.iloc[0]["consensus_peptide"], "FP_MOD")


class TestAddQualityColumns(unittest.TestCase):
    """Test quality flag generation."""

    def test_high_quality(self):
        df = pd.DataFrame({
            "ms1_rt_r2": [0.95], "ms1_im_r2": [0.90], "isotope_cosim": [0.95],
        })
        result = _add_quality_columns(df)
        self.assertTrue(result.iloc[0]["is_high_quality"])

    def test_low_quality(self):
        df = pd.DataFrame({
            "ms1_rt_r2": [0.5], "ms1_im_r2": [0.3], "isotope_cosim": [0.4],
        })
        result = _add_quality_columns(df)
        self.assertFalse(result.iloc[0]["is_high_quality"])

    def test_mixed_quality(self):
        """One metric below threshold -> not high quality."""
        df = pd.DataFrame({
            "ms1_rt_r2": [0.95], "ms1_im_r2": [0.5], "isotope_cosim": [0.95],
        })
        result = _add_quality_columns(df)
        self.assertFalse(result.iloc[0]["is_high_quality"])

    def test_boundary_values(self):
        """Exactly at threshold is high quality."""
        df = pd.DataFrame({
            "ms1_rt_r2": [0.8], "ms1_im_r2": [0.8], "isotope_cosim": [0.9],
        })
        result = _add_quality_columns(df)
        self.assertTrue(result.iloc[0]["is_high_quality"])

    def test_missing_quality_columns(self):
        """No quality columns -> all False."""
        df = pd.DataFrame({"precursor_id": [1, 2]})
        result = _add_quality_columns(df)
        self.assertTrue(all(result["is_high_quality"] == False))

    def test_nan_treated_as_low(self):
        df = pd.DataFrame({
            "ms1_rt_r2": [float("nan")], "ms1_im_r2": [0.9], "isotope_cosim": [0.95],
        })
        result = _add_quality_columns(df)
        self.assertFalse(result.iloc[0]["is_high_quality"])


class TestComputeQualitySummary(unittest.TestCase):
    """Test quality summary statistics."""

    def test_engine_distribution(self):
        df = pd.DataFrame({
            "n_engines": [0, 1, 1, 2, 2, 2, 3, 3, 3, 3],
            "is_high_quality": [False] * 6 + [True] * 4,
        })
        summary = _compute_quality_summary(df)
        self.assertEqual(summary["n_0_engines"], 1)
        self.assertEqual(summary["n_1_engines"], 2)
        self.assertEqual(summary["n_2_engines"], 3)
        self.assertEqual(summary["n_3_engines"], 4)
        self.assertEqual(summary["n_high_quality"], 4)
        self.assertAlmostEqual(summary["pct_high_quality"], 40.0)

    def test_quality_metrics_summary(self):
        df = pd.DataFrame({
            "n_engines": [1, 1],
            "ms1_rt_r2": [0.8, 0.9],
            "ms1_im_r2": [0.7, 0.8],
            "isotope_cosim": [0.9, 0.95],
            "is_high_quality": [False, True],
        })
        summary = _compute_quality_summary(df)
        self.assertAlmostEqual(summary["ms1_rt_r2_mean"], 0.85, places=2)
        self.assertIn("ms1_im_r2_median", summary)

    def test_empty_df(self):
        df = pd.DataFrame({"n_engines": [], "is_high_quality": []})
        summary = _compute_quality_summary(df)
        self.assertEqual(summary["n_0_engines"], 0)


if __name__ == "__main__":
    unittest.main(verbosity=2)
