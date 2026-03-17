#!/usr/bin/env python3
"""Tests for DIA-NN parser - RT conversion, report parsing."""

import tempfile
import unittest
from pathlib import Path

import pandas as pd
import numpy as np

import sys
sys.path.insert(0, str(Path(__file__).parent.parent))

from engine_parsers.diann_parser import DiannParser


class TestDiannParser(unittest.TestCase):
    """Test DIA-NN report.parquet parsing."""

    def setUp(self):
        """Create temp directory with synthetic DIA-NN report.parquet."""
        self.tmpdir = tempfile.mkdtemp()
        self.accession = "test_acc"
        diann_dir = Path(self.tmpdir) / self.accession / "diann"
        diann_dir.mkdir(parents=True)

        report = pd.DataFrame({
            "Run": [
                "data/raw/run1.d/run1",
                "data/raw/run1.d/run1",
                "data/raw/run2.d/run2",
            ],
            "Modified.Sequence": [
                "PEPTIDE(UniMod:35)K",
                "AC(UniMod:4)DEFGHK",
                "(UniMod:1)MPEPTIDEK",
            ],
            "Stripped.Sequence": ["PEPTIDEK", "ACDEFGHK", "MPEPTIDEK"],
            "Protein.Ids": ["P12345", "P67890", "P11111"],
            "Precursor.Charge": [2, 3, 2],
            "Precursor.Mz": [500.25, 300.15, 450.50],
            "Q.Value": [0.01, 0.05, 0.001],
            "PEP": [0.001, 0.01, 0.0005],
            "Global.Q.Value": [0.01, 0.05, 0.001],
            "PG.Q.Value": [0.02, 0.06, 0.002],
            "RT": [2.0, 4.0, 6.5],  # In MINUTES
            "IM": [0.95, 1.05, 0.88],
            "CCS": [350.5, 420.3, 310.8],
        })
        report.to_parquet(diann_dir / "report.parquet", index=False)

    def tearDown(self):
        import shutil
        shutil.rmtree(self.tmpdir)

    def test_parse_returns_dataframe(self):
        parser = DiannParser()
        result = parser.parse(Path(self.tmpdir), self.accession)
        self.assertIsInstance(result, pd.DataFrame)
        self.assertFalse(result.empty)

    def test_rt_converted_to_seconds(self):
        """CRITICAL: RT must be converted from minutes to seconds."""
        parser = DiannParser()
        result = parser.parse(Path(self.tmpdir), self.accession)
        # Input RT was 2.0 minutes -> expected 120.0 seconds
        rt_values = result["diann_rt"].tolist()
        self.assertAlmostEqual(rt_values[0], 120.0)
        self.assertAlmostEqual(rt_values[1], 240.0)
        self.assertAlmostEqual(rt_values[2], 390.0)

    def test_correct_columns(self):
        parser = DiannParser()
        result = parser.parse(Path(self.tmpdir), self.accession)
        expected_cols = [
            "raw_file", "diann_peptide", "diann_modified", "diann_protein",
            "diann_charge", "diann_mz", "diann_qvalue", "diann_pep",
            "diann_global_qvalue", "diann_pg_qvalue", "diann_rt",
            "diann_mobility", "diann_ccs",
        ]
        for col in expected_cols:
            self.assertIn(col, result.columns, f"Missing column: {col}")

    def test_raw_file_extraction(self):
        """Raw file should be extracted from Run column (stem, no .d)."""
        parser = DiannParser()
        result = parser.parse(Path(self.tmpdir), self.accession)
        raw_files = result["raw_file"].unique()
        self.assertIn("run1", raw_files)
        self.assertIn("run2", raw_files)

    def test_sequence_standardization(self):
        """Modified sequences should be converted to UNIMOD terminal format."""
        parser = DiannParser()
        result = parser.parse(Path(self.tmpdir), self.accession)
        mods = result["diann_modified"].tolist()
        # (UniMod:35) -> [UNIMOD:35] with terminal brackets
        self.assertIn("[UNIMOD:35]", mods[0])
        self.assertTrue(mods[0].startswith("[]-") or mods[0].startswith("[UNIMOD:"))
        self.assertTrue(mods[0].endswith("-[]"))

    def test_nterm_acetyl(self):
        """N-terminal acetyl should be in terminal bracket."""
        parser = DiannParser()
        result = parser.parse(Path(self.tmpdir), self.accession)
        nterm_row = result[result["diann_peptide"] == "MPEPTIDEK"].iloc[0]
        self.assertTrue(nterm_row["diann_modified"].startswith("[UNIMOD:1]-"))

    def test_ccs_preserved(self):
        parser = DiannParser()
        result = parser.parse(Path(self.tmpdir), self.accession)
        self.assertTrue(all(result["diann_ccs"].notna()))
        self.assertAlmostEqual(result.iloc[0]["diann_ccs"], 350.5)

    def test_no_results_returns_empty(self):
        parser = DiannParser()
        result = parser.parse(Path(self.tmpdir), "nonexistent")
        self.assertTrue(result.empty)

    def test_raw_file_filter(self):
        parser = DiannParser()
        result = parser.parse(Path(self.tmpdir), self.accession, raw_file_filter="run1")
        self.assertTrue(all(result["raw_file"] == "run1"))
        self.assertEqual(len(result), 2)

    def test_engine_name(self):
        self.assertEqual(DiannParser().engine_name, "diann")

    def test_has_no_precursor_id(self):
        self.assertFalse(DiannParser().has_precursor_id)


if __name__ == "__main__":
    unittest.main(verbosity=2)
