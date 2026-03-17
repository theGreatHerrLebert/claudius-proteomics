#!/usr/bin/env python3
"""Tests for Sage parser - scannr→precursor_id, mass→m/z, RT conversion."""

import tempfile
import unittest
from pathlib import Path

import pandas as pd
import numpy as np

import sys
sys.path.insert(0, str(Path(__file__).parent.parent))

from engine_parsers.sage_parser import SageParser

PROTON_MASS = 1.007276


class TestSageParser(unittest.TestCase):
    """Test Sage results.sage.parquet parsing."""

    def setUp(self):
        """Create temp directory with synthetic Sage results."""
        self.tmpdir = tempfile.mkdtemp()
        self.accession = "test_acc"
        sage_dir = Path(self.tmpdir) / self.accession / "sage"
        sage_dir.mkdir(parents=True)

        # Sage reports:
        # - expmass (neutral mass, not m/z)
        # - rt in MINUTES
        # - scannr (0-indexed, precursor_id = scannr+1)
        results = pd.DataFrame({
            "peptide": [
                "PEPTIDEK",
                "C[+57.021465]PEPTIDEK",
                "[+42.010567]-M[+15.994915]DLAAAAEPGAGSQHLEVR",
            ],
            "proteins": ["P12345", "P67890", "P11111"],
            "filename": [
                "data/raw/run1.d/run1.mzML",
                "data/raw/run1.d/run1.mzML",
                "data/raw/run2.d/run2.mzML",
            ],
            "scannr": [99, 199, 499],  # 0-indexed
            "charge": [2, 3, 2],
            "expmass": [
                900.0,   # neutral mass
                1200.0,
                2100.0,
            ],
            "rt": [2.0, 4.0, 6.5],  # MINUTES
            "ion_mobility": [0.95, 1.05, 0.88],
            "hyperscore": [25.5, 18.3, 32.1],
            "spectrum_q": [0.01, 0.05, 0.001],
            "peptide_q": [0.01, 0.05, 0.001],
            "protein_q": [0.02, 0.06, 0.002],
            "posterior_error": [np.log(0.001), np.log(0.01), np.log(0.0005)],
            "is_decoy": [False, False, False],
        })
        results.to_parquet(sage_dir / "results.sage.parquet", index=False)

        # Also create one with decoys
        results_with_decoys = results.copy()
        decoy = pd.DataFrame({
            "peptide": ["DECOYSEQK"],
            "proteins": ["DECOY_P99"],
            "filename": ["data/raw/run1.d/run1.mzML"],
            "scannr": [999], "charge": [2], "expmass": [800.0],
            "rt": [3.0], "ion_mobility": [0.90], "hyperscore": [5.0],
            "spectrum_q": [0.9], "peptide_q": [0.9], "protein_q": [0.9],
            "posterior_error": [np.log(0.9)], "is_decoy": [True],
        })
        results_with_decoys = pd.concat([results_with_decoys, decoy], ignore_index=True)
        sage_dir2 = Path(self.tmpdir) / "test_decoy" / "sage"
        sage_dir2.mkdir(parents=True)
        results_with_decoys.to_parquet(sage_dir2 / "results.sage.parquet", index=False)

    def tearDown(self):
        import shutil
        shutil.rmtree(self.tmpdir)

    def test_parse_returns_dataframe(self):
        parser = SageParser()
        result = parser.parse(Path(self.tmpdir), self.accession)
        self.assertIsInstance(result, pd.DataFrame)
        self.assertEqual(len(result), 3)

    def test_precursor_id_from_scannr(self):
        """CRITICAL: precursor_id = scannr + 1 (1-indexed)."""
        parser = SageParser()
        result = parser.parse(Path(self.tmpdir), self.accession)
        pids = result["precursor_id"].tolist()
        self.assertEqual(pids[0], 100)   # scannr=99 -> pid=100
        self.assertEqual(pids[1], 200)   # scannr=199 -> pid=200
        self.assertEqual(pids[2], 500)   # scannr=499 -> pid=500

    def test_mz_calculated_from_mass(self):
        """CRITICAL: m/z = (expmass + charge * proton) / charge."""
        parser = SageParser()
        result = parser.parse(Path(self.tmpdir), self.accession)

        for _, row in result.iterrows():
            # Reverse-calculate expected m/z
            charge = row["sage_charge"]
            # Find corresponding expmass
            idx = result.index.get_loc(row.name)
            expmass = [900.0, 1200.0, 2100.0][idx]
            expected_mz = (expmass + charge * PROTON_MASS) / charge
            self.assertAlmostEqual(row["sage_mz"], expected_mz, places=3)

    def test_rt_converted_to_seconds(self):
        """CRITICAL: RT must be converted from minutes to seconds."""
        parser = SageParser()
        result = parser.parse(Path(self.tmpdir), self.accession)
        rt_values = result["sage_rt"].tolist()
        self.assertAlmostEqual(rt_values[0], 120.0)  # 2.0 min -> 120s
        self.assertAlmostEqual(rt_values[1], 240.0)
        self.assertAlmostEqual(rt_values[2], 390.0)

    def test_correct_columns(self):
        parser = SageParser()
        result = parser.parse(Path(self.tmpdir), self.accession)
        expected_cols = [
            "raw_file", "precursor_id", "sage_psm_id",
            "sage_peptide", "sage_modified", "sage_protein", "sage_charge",
            "sage_mz", "sage_qvalue", "sage_peptide_qvalue", "sage_protein_qvalue",
            "sage_pep", "sage_hyperscore", "sage_rt", "sage_mobility",
        ]
        for col in expected_cols:
            self.assertIn(col, result.columns, f"Missing column: {col}")

    def test_raw_file_extraction(self):
        parser = SageParser()
        result = parser.parse(Path(self.tmpdir), self.accession)
        raw_files = result["raw_file"].unique()
        self.assertIn("run1", raw_files)
        self.assertIn("run2", raw_files)

    def test_decoys_filtered(self):
        """Decoy hits should be removed."""
        parser = SageParser()
        result = parser.parse(Path(self.tmpdir), "test_decoy")
        # Should have 3 targets, no decoy
        self.assertEqual(len(result), 3)
        self.assertNotIn("DECOYSEQK", result["sage_peptide"].values)

    def test_sequence_standardization(self):
        """Mass deltas should be converted to UNIMOD terminal format."""
        parser = SageParser()
        result = parser.parse(Path(self.tmpdir), self.accession)
        mods = result["sage_modified"].tolist()
        # C[+57.021465] -> C[UNIMOD:4] in terminal format
        self.assertIn("[UNIMOD:4]", mods[1])
        self.assertTrue(mods[1].endswith("-[]"))

    def test_nterm_mod_standardized(self):
        """Sage N-term [+42.010567]-M... should become [UNIMOD:1]-M...-[]."""
        parser = SageParser()
        result = parser.parse(Path(self.tmpdir), self.accession)
        nterm_row = result.iloc[2]
        self.assertTrue(nterm_row["sage_modified"].startswith("[UNIMOD:1]-"))
        self.assertIn("[UNIMOD:35]", nterm_row["sage_modified"])

    def test_pep_from_posterior_error(self):
        """PEP should be exp(posterior_error)."""
        parser = SageParser()
        result = parser.parse(Path(self.tmpdir), self.accession)
        self.assertAlmostEqual(result.iloc[0]["sage_pep"], 0.001, places=5)
        self.assertAlmostEqual(result.iloc[1]["sage_pep"], 0.01, places=5)

    def test_stripped_peptide(self):
        """sage_peptide should be plain sequence without mods."""
        parser = SageParser()
        result = parser.parse(Path(self.tmpdir), self.accession)
        self.assertEqual(result.iloc[0]["sage_peptide"], "PEPTIDEK")
        self.assertEqual(result.iloc[1]["sage_peptide"], "CPEPTIDEK")

    def test_no_results_returns_empty(self):
        parser = SageParser()
        result = parser.parse(Path(self.tmpdir), "nonexistent")
        self.assertTrue(result.empty)

    def test_raw_file_filter(self):
        parser = SageParser()
        result = parser.parse(Path(self.tmpdir), self.accession, raw_file_filter="run1")
        self.assertTrue(all(result["raw_file"] == "run1"))
        self.assertEqual(len(result), 2)

    def test_engine_name(self):
        self.assertEqual(SageParser().engine_name, "sage")

    def test_has_precursor_id(self):
        self.assertTrue(SageParser().has_precursor_id)


if __name__ == "__main__":
    unittest.main(verbosity=2)
