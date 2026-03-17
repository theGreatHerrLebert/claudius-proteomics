#!/usr/bin/env python3
"""
Tests for step3_stratify.py - precursor index building and engine matching.

Tests the core stratification logic:
- _build_precursor_index_anchored: anchor-based 4-step merge
- _build_precursor_index_engines_only: fallback when no raw features
- _match_engine_to_index: tiered sequence + coordinate matching
- _compute_overlap_stats_from_parsers: 3-way Venn overlap calculation
"""

import unittest
import pandas as pd
import numpy as np

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from runner.steps.step3_stratify import (
    _build_precursor_index_anchored,
    _build_precursor_index_engines_only,
    _compute_overlap_stats_from_parsers,
    _match_engine_to_index,
)
from scripts.precursor_merging.config import MatchConfig


def _make_raw(n=5, raw_file="run1"):
    """Create synthetic raw precursors."""
    return pd.DataFrame({
        "precursor_id": range(1, n + 1),
        "raw_file": [raw_file] * n,
        "mz": [500.0 + i * 10 for i in range(n)],
        "charge": [2] * n,
        "rt_seconds": [100.0 + i * 50 for i in range(n)],
        "mobility": [0.90 + i * 0.02 for i in range(n)],
    })


def _make_fragpipe(ids, raw_file="run1", peptides=None):
    """Create synthetic FragPipe results."""
    n = len(ids)
    if peptides is None:
        peptides = [f"PEPTIDE{i}K" for i in ids]
    return pd.DataFrame({
        "raw_file": [raw_file] * n,
        "precursor_id": list(ids),
        "fragpipe_peptide": peptides,
        "fragpipe_modified": [f"[]-{p}-[]" for p in peptides],
        "fragpipe_protein": [f"P{i}" for i in ids],
        "fragpipe_probability": [0.99] * n,
        "fragpipe_pep": [0.01] * n,
        "fragpipe_charge": [2] * n,
        "fragpipe_mz": [500.0 + (i - 1) * 10 for i in ids],
        "fragpipe_rt": [100.0 + (i - 1) * 50 for i in ids],
        "fragpipe_mobility": [0.90 + (i - 1) * 0.02 for i in ids],
    })


def _make_diann(peptides, raw_file="run1", charges=None, mzs=None, rts=None, ims=None):
    """Create synthetic DIA-NN results."""
    n = len(peptides)
    return pd.DataFrame({
        "raw_file": [raw_file] * n,
        "diann_peptide": peptides,
        "diann_modified": [f"[]-{p}-[]" for p in peptides],
        "diann_protein": [f"D{i}" for i in range(n)],
        "diann_charge": charges or [2] * n,
        "diann_mz": mzs or [500.0 + i * 10 for i in range(n)],
        "diann_rt": rts or [100.0 + i * 50 for i in range(n)],
        "diann_mobility": ims or [0.90 + i * 0.02 for i in range(n)],
        "diann_qvalue": [0.01] * n,
        "diann_pep": [0.001] * n,
        "diann_ccs": [350.0] * n,
    })


def _make_sage(ids, raw_file="run1", peptides=None):
    """Create synthetic Sage results."""
    n = len(ids)
    if peptides is None:
        peptides = [f"PEPTIDE{i}K" for i in ids]
    return pd.DataFrame({
        "raw_file": [raw_file] * n,
        "precursor_id": list(ids),
        "sage_peptide": peptides,
        "sage_modified": [f"[]-{p}-[]" for p in peptides],
        "sage_protein": [f"S{i}" for i in ids],
        "sage_charge": [2] * n,
        "sage_mz": [500.0 + (i - 1) * 10 for i in ids],
        "sage_qvalue": [0.01] * n,
        "sage_pep": [0.001] * n,
        "sage_hyperscore": [25.0] * n,
        "sage_rt": [100.0 + (i - 1) * 50 for i in ids],
        "sage_mobility": [0.90 + (i - 1) * 0.02 for i in ids],
    })


class TestBuildPrecursorIndexAnchored(unittest.TestCase):
    """Test the main anchored index builder."""

    def test_basic_merge(self):
        """All three engines join to raw precursors."""
        raw = _make_raw(3)
        fp = _make_fragpipe([1, 2, 3])
        dn = _make_diann(["PEPTIDE1K", "PEPTIDE2K", "PEPTIDE3K"])
        sg = _make_sage([1, 2, 3])

        result = _build_precursor_index_anchored(raw, fp, dn, sg, {})
        self.assertEqual(len(result), 3)
        self.assertIn("n_engines", result.columns)
        self.assertTrue(all(result["n_engines"] == 3))

    def test_fragpipe_direct_join(self):
        """FragPipe joins by (raw_file, precursor_id)."""
        raw = _make_raw(5)
        fp = _make_fragpipe([1, 3, 5])  # Only 3 of 5

        result = _build_precursor_index_anchored(raw, fp, None, None, {})
        self.assertEqual(len(result), 5)
        n_fp = result["fragpipe_peptide"].notna().sum()
        self.assertEqual(n_fp, 3)

    def test_sage_direct_join(self):
        """Sage joins by (raw_file, precursor_id)."""
        raw = _make_raw(5)
        sg = _make_sage([2, 4])

        result = _build_precursor_index_anchored(raw, None, None, sg, {})
        self.assertEqual(len(result), 5)
        n_sg = result["sage_peptide"].notna().sum()
        self.assertEqual(n_sg, 2)

    def test_diann_sequence_match(self):
        """DIA-NN matches by normalized sequence + charge."""
        raw = _make_raw(3)
        fp = _make_fragpipe([1, 2, 3])
        # DIA-NN has same peptides -> should match by sequence
        dn = _make_diann(["PEPTIDE1K", "PEPTIDE2K"])

        result = _build_precursor_index_anchored(raw, fp, dn, None, {})
        n_dn = result["diann_peptide"].notna().sum()
        self.assertEqual(n_dn, 2)

    def test_engine_count(self):
        """n_engines should reflect actual matches."""
        raw = _make_raw(4)
        fp = _make_fragpipe([1, 2, 3, 4])
        dn = _make_diann(["PEPTIDE1K", "PEPTIDE2K"])
        sg = _make_sage([1])

        result = _build_precursor_index_anchored(raw, fp, dn, sg, {})
        engine_counts = result["n_engines"].value_counts().to_dict()
        self.assertEqual(engine_counts.get(3, 0), 1)  # pid=1: all three
        self.assertEqual(engine_counts.get(2, 0), 1)  # pid=2: fp+dn
        self.assertEqual(engine_counts.get(1, 0), 2)  # pid=3,4: fp only

    def test_empty_raw_precursors(self):
        """Empty raw -> falls back to engines-only."""
        fp = _make_fragpipe([1, 2])
        result = _build_precursor_index_anchored(pd.DataFrame(), fp, None, None, {})
        self.assertGreater(len(result), 0)

    def test_no_engines(self):
        """Only raw precursors, no engines."""
        raw = _make_raw(3)
        result = _build_precursor_index_anchored(raw, None, None, None, {})
        self.assertEqual(len(result), 3)
        self.assertTrue(all(result["n_engines"] == 0))

    def test_raw_file_normalization(self):
        """'.d' extension should be stripped from raw_file."""
        raw = _make_raw(2, raw_file="run1.d")
        fp = _make_fragpipe([1, 2], raw_file="run1")

        result = _build_precursor_index_anchored(raw, fp, None, None, {})
        n_fp = result["fragpipe_peptide"].notna().sum()
        self.assertEqual(n_fp, 2)

    def test_preserves_raw_columns(self):
        """Anchor columns from raw precursors are preserved."""
        raw = _make_raw(2)
        result = _build_precursor_index_anchored(raw, None, None, None, {})
        for col in ["precursor_id", "raw_file", "mz", "charge", "rt_seconds", "mobility"]:
            self.assertIn(col, result.columns)

    def test_missing_precursor_id_raises(self):
        """Should raise if raw_precursors lacks precursor_id."""
        raw = pd.DataFrame({"raw_file": ["run1"], "mz": [500.0]})
        with self.assertRaises(ValueError):
            _build_precursor_index_anchored(raw, None, None, None, {})


class TestBuildPrecursorIndexEnginesOnly(unittest.TestCase):
    """Test engine-only fallback index builder."""

    def test_single_engine(self):
        fp = _make_fragpipe([1, 2])
        result = _build_precursor_index_engines_only(fp, None, None)
        self.assertGreater(len(result), 0)

    def test_all_engines(self):
        fp = _make_fragpipe([1, 2], peptides=["AAAK", "BBBK"])
        dn = _make_diann(["AAAK", "CCCK"])
        sg = _make_sage([1, 2], peptides=["AAAK", "DDDK"])
        result = _build_precursor_index_engines_only(fp, dn, sg)
        self.assertGreater(len(result), 0)
        self.assertIn("n_engines", result.columns)

    def test_no_engines(self):
        result = _build_precursor_index_engines_only(None, None, None)
        self.assertTrue(result.empty or len(result) == 0)


class TestComputeOverlapStats(unittest.TestCase):
    """Test 3-way overlap statistics."""

    def test_full_overlap(self):
        """All three engines identify same peptide."""
        fp = pd.DataFrame({"fragpipe_modified": ["PEPTIDEK"], "fragpipe_charge": [2]})
        dn = pd.DataFrame({"diann_modified": ["PEPTIDEK"], "diann_charge": [2]})
        sg = pd.DataFrame({"sage_modified": ["PEPTIDEK"], "sage_charge": [2]})
        stats = _compute_overlap_stats_from_parsers(fp, dn, sg)
        self.assertEqual(stats["n_all_three"], 1)
        self.assertEqual(stats["n_union"], 1)
        self.assertAlmostEqual(stats["three_way_rate"], 1.0)

    def test_no_overlap(self):
        """Each engine identifies different peptide."""
        fp = pd.DataFrame({"fragpipe_modified": ["AAAK"], "fragpipe_charge": [2]})
        dn = pd.DataFrame({"diann_modified": ["BBBK"], "diann_charge": [2]})
        sg = pd.DataFrame({"sage_modified": ["CCCK"], "sage_charge": [2]})
        stats = _compute_overlap_stats_from_parsers(fp, dn, sg)
        self.assertEqual(stats["n_all_three"], 0)
        self.assertEqual(stats["n_union"], 3)
        self.assertEqual(stats["n_fragpipe_only"], 1)
        self.assertEqual(stats["n_diann_only"], 1)
        self.assertEqual(stats["n_sage_only"], 1)

    def test_partial_overlap(self):
        """FragPipe+Sage overlap, DIA-NN unique."""
        fp = pd.DataFrame({"fragpipe_modified": ["AAAK", "BBBK"], "fragpipe_charge": [2, 2]})
        dn = pd.DataFrame({"diann_modified": ["CCCK"], "diann_charge": [2]})
        sg = pd.DataFrame({"sage_modified": ["AAAK"], "sage_charge": [2]})
        stats = _compute_overlap_stats_from_parsers(fp, dn, sg)
        self.assertEqual(stats["n_all_three"], 0)
        self.assertEqual(stats["n_fp_sg_only"], 1)  # AAAK in FP+SG
        self.assertEqual(stats["n_at_least_two"], 1)

    def test_empty_engines(self):
        stats = _compute_overlap_stats_from_parsers(None, None, None)
        self.assertEqual(stats["n_union"], 0)
        self.assertEqual(stats["n_all_three"], 0)

    def test_charge_matters(self):
        """Same sequence with different charge = different precursors."""
        fp = pd.DataFrame({"fragpipe_modified": ["AAAK"], "fragpipe_charge": [2]})
        dn = pd.DataFrame({"diann_modified": ["AAAK"], "diann_charge": [3]})
        stats = _compute_overlap_stats_from_parsers(fp, dn, None)
        self.assertEqual(stats["n_all_three"], 0)
        self.assertEqual(stats["n_union"], 2)  # Different charge = different


class TestMatchEngineToIndex(unittest.TestCase):
    """Test tiered engine matching to existing index."""

    def test_sequence_match(self):
        """DIA-NN matched by sequence when FragPipe ID available."""
        raw = _make_raw(2)
        fp = _make_fragpipe([1, 2])
        index = _build_precursor_index_anchored(raw, fp, None, None, {})

        dn = _make_diann(["PEPTIDE1K"])
        config = MatchConfig()
        result, stats = _match_engine_to_index(
            index, dn, "diann", "diann_modified", "diann_charge",
            "diann_mz", "diann_rt", "diann_mobility", config,
        )
        self.assertEqual(stats["sequence"], 1)
        self.assertEqual(result["diann_peptide"].notna().sum(), 1)

    def test_coordinate_match(self):
        """Coordinate fallback when no sequence available."""
        # Raw precursors without FragPipe IDs
        raw = _make_raw(2)
        index = raw.copy()
        index["sequence_normalized"] = ""
        index["fragpipe_peptide"] = None

        # DIA-NN with matching coordinates
        dn = _make_diann(
            ["NEWPEPTK"],
            mzs=[500.0],
            rts=[100.0],
            ims=[0.90],
        )
        config = MatchConfig(rt_tol_sec=5.0, im_tol=0.05)
        result, stats = _match_engine_to_index(
            index, dn, "diann", "diann_modified", "diann_charge",
            "diann_mz", "diann_rt", "diann_mobility", config,
        )
        self.assertEqual(stats["coordinate"], 1)

    def test_no_match_outside_tolerance(self):
        """No match when coordinates outside tolerance."""
        raw = _make_raw(1)
        index = raw.copy()
        index["sequence_normalized"] = ""
        index["fragpipe_peptide"] = None

        dn = _make_diann(
            ["FARAWAYK"],
            mzs=[999.0],  # Way off
            rts=[100.0],
            ims=[0.90],
        )
        config = MatchConfig(rt_tol_sec=5.0)
        result, stats = _match_engine_to_index(
            index, dn, "diann", "diann_modified", "diann_charge",
            "diann_mz", "diann_rt", "diann_mobility", config,
        )
        self.assertEqual(stats["sequence"], 0)
        self.assertEqual(stats["coordinate"], 0)


if __name__ == "__main__":
    unittest.main(verbosity=2)
