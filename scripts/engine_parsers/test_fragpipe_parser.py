#!/usr/bin/env python3
"""Tests for FragPipe parser - parse_spectrum_id and FragPipeParser.parse."""

import os
import tempfile
import unittest
from pathlib import Path

import pandas as pd

import sys
sys.path.insert(0, str(Path(__file__).parent.parent))

from engine_parsers.fragpipe_parser import parse_spectrum_id, FragPipeParser


class TestParseSpectrumId(unittest.TestCase):
    """Test FragPipe Spectrum column parsing."""

    def test_standard_format(self):
        """Standard format: rawfile.scannum.scannum.charge"""
        raw, pid, charge = parse_spectrum_id("sample01.12345.12345.2")
        self.assertEqual(raw, "sample01")
        self.assertEqual(pid, 12345)
        self.assertEqual(charge, 2)

    def test_charge_3(self):
        raw, pid, charge = parse_spectrum_id("run1.999.999.3")
        self.assertEqual(charge, 3)

    def test_dots_in_filename(self):
        """Raw file name with dots (e.g., file.name.d)."""
        raw, pid, charge = parse_spectrum_id("20190504_TIMS1_FlMe_SA_HeLa_frac01_A10_1_93.5678.5678.2")
        self.assertEqual(raw, "20190504_TIMS1_FlMe_SA_HeLa_frac01_A10_1_93")
        self.assertEqual(pid, 5678)
        self.assertEqual(charge, 2)

    def test_large_scan_numbers(self):
        raw, pid, charge = parse_spectrum_id("sample.999999.999999.4")
        self.assertEqual(pid, 999999)
        self.assertEqual(charge, 4)

    def test_none_input(self):
        self.assertEqual(parse_spectrum_id(None), (None, None, None))

    def test_empty_string(self):
        self.assertEqual(parse_spectrum_id(""), (None, None, None))

    def test_nan_input(self):
        self.assertEqual(parse_spectrum_id(float('nan')), (None, None, None))

    def test_too_few_parts(self):
        """Only 2 dots -> not enough parts."""
        self.assertEqual(parse_spectrum_id("file.123.2"), (None, None, None))

    def test_non_numeric_scan(self):
        """Non-numeric scan number -> raw_file returned, others None."""
        raw, pid, charge = parse_spectrum_id("file.abc.abc.2")
        self.assertEqual(raw, "file")
        self.assertIsNone(pid)
        self.assertIsNone(charge)

    def test_single_part(self):
        self.assertEqual(parse_spectrum_id("nodelimiters"), (None, None, None))


class TestFragPipeParser(unittest.TestCase):
    """Test FragPipeParser with synthetic psm.tsv data."""

    def setUp(self):
        """Create temp directory with minimal psm.tsv."""
        self.tmpdir = tempfile.mkdtemp()
        self.accession = "test_acc"
        acc_dir = Path(self.tmpdir) / self.accession
        acc_dir.mkdir()

        psm_data = pd.DataFrame({
            "Spectrum": [
                "run1.100.100.2",
                "run1.200.200.3",
                "run1.100.100.2",  # Duplicate precursor_id (lower prob should be dropped)
            ],
            "Peptide": ["PEPTIDEK", "ACMEK", "PEPTIDEK"],
            "Assigned Modifications": ["", "2C(57.0215)", ""],
            "Protein": ["P12345", "P67890", "P12345"],
            "Probability": [0.99, 0.85, 0.50],
            "Charge": [2, 3, 2],
            "Observed M/Z": [500.25, 300.15, 500.25],
            "Retention": [120.5, 240.3, 120.5],
            "Ion Mobility": [0.95, 1.05, 0.95],
        })
        psm_data.to_csv(acc_dir / "psm.tsv", sep="\t", index=False)

    def tearDown(self):
        import shutil
        shutil.rmtree(self.tmpdir)

    def test_parse_returns_dataframe(self):
        parser = FragPipeParser()
        result = parser.parse(Path(self.tmpdir), self.accession)
        self.assertIsInstance(result, pd.DataFrame)
        self.assertFalse(result.empty)

    def test_deduplication(self):
        """Duplicate precursor_id keeps highest probability."""
        parser = FragPipeParser()
        result = parser.parse(Path(self.tmpdir), self.accession)
        # run1.100.100.2 appears twice, should keep prob=0.99
        run1_100 = result[result["precursor_id"] == 100]
        self.assertEqual(len(run1_100), 1)
        self.assertAlmostEqual(run1_100.iloc[0]["fragpipe_probability"], 0.99)

    def test_correct_columns(self):
        parser = FragPipeParser()
        result = parser.parse(Path(self.tmpdir), self.accession)
        expected_cols = [
            "raw_file", "precursor_id", "fragpipe_peptide", "fragpipe_modified",
            "fragpipe_protein", "fragpipe_probability", "fragpipe_pep",
            "fragpipe_mz", "fragpipe_rt", "fragpipe_mobility", "fragpipe_charge",
        ]
        for col in expected_cols:
            self.assertIn(col, result.columns, f"Missing column: {col}")

    def test_pep_is_1_minus_prob(self):
        parser = FragPipeParser()
        result = parser.parse(Path(self.tmpdir), self.accession)
        for _, row in result.iterrows():
            self.assertAlmostEqual(
                row["fragpipe_pep"] + row["fragpipe_probability"], 1.0, places=5
            )

    def test_raw_file_extraction(self):
        parser = FragPipeParser()
        result = parser.parse(Path(self.tmpdir), self.accession)
        self.assertTrue(all(result["raw_file"] == "run1"))

    def test_no_results_returns_empty(self):
        parser = FragPipeParser()
        result = parser.parse(Path(self.tmpdir), "nonexistent")
        self.assertTrue(result.empty)

    def test_raw_file_filter(self):
        parser = FragPipeParser()
        result = parser.parse(Path(self.tmpdir), self.accession, raw_file_filter="run1")
        self.assertFalse(result.empty)
        result2 = parser.parse(Path(self.tmpdir), self.accession, raw_file_filter="nonexistent")
        self.assertTrue(result2.empty)

    def test_engine_name(self):
        parser = FragPipeParser()
        self.assertEqual(parser.engine_name, "fragpipe")

    def test_has_precursor_id(self):
        parser = FragPipeParser()
        self.assertTrue(parser.has_precursor_id)


if __name__ == "__main__":
    unittest.main(verbosity=2)
