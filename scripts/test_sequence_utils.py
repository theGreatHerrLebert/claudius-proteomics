#!/usr/bin/env python3
"""
Tests for sequence_utils.py

Run with: python -m pytest test_sequence_utils.py -v
Or simply: python test_sequence_utils.py
"""

import unittest
from sequence_utils import (
    mass_to_unimod,
    standardize_diann_sequence,
    standardize_sage_sequence,
    standardize_fragpipe_sequence,
    standardize_fragpipe_modified_peptide,
    remove_modifications,
    normalize_sequence_il,
    sequences_match,
    standardize_modified_sequence,
)


class TestMassToUnimod(unittest.TestCase):
    """Test mass to UNIMOD conversion."""

    def test_carbamidomethyl(self):
        """57 Da -> UNIMOD:4 (Carbamidomethyl)"""
        self.assertEqual(mass_to_unimod(57.021), "[UNIMOD:4]")
        self.assertEqual(mass_to_unimod(57.0), "[UNIMOD:4]")
        self.assertEqual(mass_to_unimod(57), "[UNIMOD:4]")

    def test_oxidation(self):
        """16 Da -> UNIMOD:35 (Oxidation)"""
        self.assertEqual(mass_to_unimod(15.9949), "[UNIMOD:35]")
        self.assertEqual(mass_to_unimod(16.0), "[UNIMOD:35]")

    def test_acetyl(self):
        """42 Da -> UNIMOD:1 (Acetyl)"""
        self.assertEqual(mass_to_unimod(42.0106), "[UNIMOD:1]")
        self.assertEqual(mass_to_unimod(42), "[UNIMOD:1]")

    def test_phospho(self):
        """80 Da -> UNIMOD:21 (Phospho)"""
        self.assertEqual(mass_to_unimod(79.9663), "[UNIMOD:21]")
        self.assertEqual(mass_to_unimod(80), "[UNIMOD:21]")

    def test_unknown_mass(self):
        """Unknown mass returns None."""
        self.assertIsNone(mass_to_unimod(999))
        self.assertIsNone(mass_to_unimod(123.456))


class TestDiannStandardization(unittest.TestCase):
    """Test DIA-NN sequence standardization."""

    def test_simple_oxidation(self):
        """DIA-NN oxidation format."""
        result = standardize_diann_sequence("PEPTIDE(UniMod:35)K")
        self.assertEqual(result, "PEPTIDE[UNIMOD:35]K")

    def test_carbamidomethyl(self):
        """DIA-NN carbamidomethyl format."""
        result = standardize_diann_sequence("C(UniMod:4)PEPTIDEK")
        self.assertEqual(result, "C[UNIMOD:4]PEPTIDEK")

    def test_multiple_mods(self):
        """Multiple modifications."""
        result = standardize_diann_sequence("C(UniMod:4)PEPTIM(UniMod:35)DEK")
        self.assertEqual(result, "C[UNIMOD:4]PEPTIM[UNIMOD:35]DEK")

    def test_nterm_acetyl(self):
        """N-terminal acetylation."""
        result = standardize_diann_sequence("(UniMod:1)PEPTIDEK")
        self.assertEqual(result, "[UNIMOD:1]PEPTIDEK")

    def test_unmodified(self):
        """Unmodified sequence passes through."""
        result = standardize_diann_sequence("PEPTIDEK")
        self.assertEqual(result, "PEPTIDEK")

    def test_empty_input(self):
        """Empty/None input returns empty string."""
        self.assertEqual(standardize_diann_sequence(""), "")
        self.assertEqual(standardize_diann_sequence(None), "")

    def test_case_insensitive_unimod(self):
        """Handle different UniMod casing."""
        result = standardize_diann_sequence("PEPTIDE(unimod:35)K")
        self.assertEqual(result, "PEPTIDE[UNIMOD:35]K")


class TestSageStandardization(unittest.TestCase):
    """Test Sage sequence standardization."""

    def test_oxidation(self):
        """Sage oxidation format [+15.9949]."""
        result = standardize_sage_sequence("PEPTIDE[+15.9949]K")
        self.assertEqual(result, "PEPTIDE[UNIMOD:35]K")

    def test_carbamidomethyl(self):
        """Sage carbamidomethyl format [+57.021465]."""
        result = standardize_sage_sequence("C[+57.021465]PEPTIDEK")
        self.assertEqual(result, "C[UNIMOD:4]PEPTIDEK")

    def test_multiple_mods(self):
        """Multiple modifications."""
        result = standardize_sage_sequence("C[+57.021]PEPTIM[+15.995]DEK")
        self.assertEqual(result, "C[UNIMOD:4]PEPTIM[UNIMOD:35]DEK")

    def test_nterm_acetyl(self):
        """N-terminal acetylation."""
        result = standardize_sage_sequence("[+42.0106]PEPTIDEK")
        self.assertEqual(result, "[UNIMOD:1]PEPTIDEK")

    def test_unmodified(self):
        """Unmodified sequence passes through."""
        result = standardize_sage_sequence("PEPTIDEK")
        self.assertEqual(result, "PEPTIDEK")

    def test_unknown_mass_preserved(self):
        """Unknown mass modifications are preserved."""
        result = standardize_sage_sequence("PEPTIDE[+999.999]K")
        self.assertEqual(result, "PEPTIDE[+999.999]K")

    def test_empty_input(self):
        """Empty/None input returns empty string."""
        self.assertEqual(standardize_sage_sequence(""), "")
        self.assertEqual(standardize_sage_sequence(None), "")


class TestFragpipeStandardization(unittest.TestCase):
    """Test FragPipe sequence standardization."""

    def test_position_based_carbamidomethyl(self):
        """FragPipe position-based format: 1C(57.0215)."""
        result = standardize_fragpipe_sequence("CPEPTIDEK", "1C(57.0215)")
        self.assertEqual(result, "C[UNIMOD:4]PEPTIDEK")

    def test_position_based_oxidation(self):
        """FragPipe position-based oxidation."""
        result = standardize_fragpipe_sequence("PEPTIMEK", "6M(15.9949)")
        self.assertEqual(result, "PEPTIM[UNIMOD:35]EK")

    def test_multiple_mods(self):
        """Multiple position-based modifications."""
        result = standardize_fragpipe_sequence("CPEPTIMEK", "1C(57.0215), 7M(15.9949)")
        self.assertEqual(result, "C[UNIMOD:4]PEPTIM[UNIMOD:35]EK")

    def test_nterm_acetyl(self):
        """N-terminal acetylation."""
        result = standardize_fragpipe_sequence("PEPTIDEK", "N-term(42.0106)")
        self.assertEqual(result, "[UNIMOD:1]PEPTIDEK")

    def test_nterm_plus_internal(self):
        """N-terminal plus internal modification."""
        result = standardize_fragpipe_sequence("CPEPTIDEK", "N-term(42.0106), 1C(57.0215)")
        self.assertEqual(result, "[UNIMOD:1]C[UNIMOD:4]PEPTIDEK")

    def test_multiple_cys_no_space(self):
        """Multiple C modifications without space after comma."""
        result = standardize_fragpipe_sequence("CCDNDMER", "1C(57.0215),2C(57.0215),6M(15.9949)")
        self.assertEqual(result, "C[UNIMOD:4]C[UNIMOD:4]DNDM[UNIMOD:35]ER")

    def test_multiple_mods_mixed_spacing(self):
        """Multiple modifications with mixed spacing."""
        result = standardize_fragpipe_sequence("ACMEK", "1A(42.0106),2C(57.0215), 3M(15.9949)")
        self.assertEqual(result, "A[UNIMOD:1]C[UNIMOD:4]M[UNIMOD:35]EK")

    def test_unmodified(self):
        """Unmodified sequence."""
        result = standardize_fragpipe_sequence("PEPTIDEK", "")
        self.assertEqual(result, "PEPTIDEK")
        result = standardize_fragpipe_sequence("PEPTIDEK", None)
        self.assertEqual(result, "PEPTIDEK")

    def test_empty_input(self):
        """Empty input returns empty string."""
        self.assertEqual(standardize_fragpipe_sequence("", ""), "")
        self.assertEqual(standardize_fragpipe_sequence(None, None), "")


class TestFragpipeModifiedPeptide(unittest.TestCase):
    """Test FragPipe Modified Peptide column standardization."""

    def test_oxidized_met(self):
        """Oxidized methionine M[147.0354]."""
        result = standardize_fragpipe_modified_peptide("M[147.0354]PEPTIDEK")
        self.assertEqual(result, "M[UNIMOD:35]PEPTIDEK")

    def test_nterm_lowercase(self):
        """N-terminal with lowercase 'n' marker."""
        result = standardize_fragpipe_modified_peptide("n[42.0106]PEPTIDEK")
        self.assertEqual(result, "[UNIMOD:1]PEPTIDEK")

    def test_unmodified(self):
        """Unmodified sequence passes through."""
        result = standardize_fragpipe_modified_peptide("PEPTIDEK")
        self.assertEqual(result, "PEPTIDEK")

    def test_empty_input(self):
        """Empty/None input returns empty string."""
        self.assertEqual(standardize_fragpipe_modified_peptide(""), "")
        self.assertEqual(standardize_fragpipe_modified_peptide(None), "")


class TestRemoveModifications(unittest.TestCase):
    """Test modification removal."""

    def test_unimod_format(self):
        """Remove [UNIMOD:X] format."""
        result = remove_modifications("C[UNIMOD:4]PEPTIM[UNIMOD:35]DEK")
        self.assertEqual(result, "CPEPTIMDEK")

    def test_mass_format(self):
        """Remove [+X.XXX] format."""
        result = remove_modifications("C[+57.021]PEPTIM[+15.995]DEK")
        self.assertEqual(result, "CPEPTIMDEK")

    def test_diann_format(self):
        """Remove (UniMod:X) format."""
        result = remove_modifications("C(UniMod:4)PEPTIM(UniMod:35)DEK")
        self.assertEqual(result, "CPEPTIMDEK")

    def test_mixed_formats(self):
        """Handle mixed formats."""
        result = remove_modifications("[UNIMOD:1]C[+57.021]PEPTIDEK")
        self.assertEqual(result, "CPEPTIDEK")

    def test_unmodified(self):
        """Unmodified sequence passes through."""
        result = remove_modifications("PEPTIDEK")
        self.assertEqual(result, "PEPTIDEK")

    def test_empty_input(self):
        """Empty/None input returns empty string."""
        self.assertEqual(remove_modifications(""), "")
        self.assertEqual(remove_modifications(None), "")


class TestNormalizeSequenceIL(unittest.TestCase):
    """Test I/L normalization."""

    def test_i_to_l(self):
        """I is replaced with L."""
        result = normalize_sequence_il("PEPTIDEK")
        self.assertEqual(result, "PEPTLDEK")

    def test_multiple_i(self):
        """Multiple I's are replaced."""
        result = normalize_sequence_il("IIIIII")
        self.assertEqual(result, "LLLLLL")

    def test_with_modifications(self):
        """Modifications are removed, I→L."""
        result = normalize_sequence_il("PEPTIDE[UNIMOD:35]K")
        self.assertEqual(result, "PEPTLDEK")

    def test_lowercase(self):
        """Lowercase is converted to uppercase."""
        result = normalize_sequence_il("peptidek")
        self.assertEqual(result, "PEPTLDEK")

    def test_empty_input(self):
        """Empty input returns empty string."""
        self.assertEqual(normalize_sequence_il(""), "")


class TestSequencesMatch(unittest.TestCase):
    """Test sequence matching."""

    def test_identical(self):
        """Identical sequences match."""
        self.assertTrue(sequences_match("PEPTIDEK", "PEPTIDEK"))

    def test_il_difference(self):
        """I and L are considered equivalent."""
        self.assertTrue(sequences_match("PEPTIDEK", "PEPTLDEK"))
        self.assertTrue(sequences_match("IIIIII", "LLLLLL"))

    def test_with_modifications(self):
        """Sequences with different mod formats match if base sequence same."""
        self.assertTrue(sequences_match(
            "PEPTIDE[UNIMOD:35]K",
            "PEPTIDE(UniMod:35)K"
        ))

    def test_different_sequences(self):
        """Different sequences don't match."""
        self.assertFalse(sequences_match("PEPTIDEK", "PEPTIDEA"))

    def test_different_lengths(self):
        """Different length sequences don't match."""
        self.assertFalse(sequences_match("PEPTIDEK", "PEPTIDE"))


class TestStandardizeModifiedSequence(unittest.TestCase):
    """Test the unified standardization function."""

    def test_diann_source(self):
        """DIA-NN source."""
        result = standardize_modified_sequence("PEPTIDE(UniMod:35)K", "diann")
        self.assertEqual(result, "PEPTIDE[UNIMOD:35]K")

    def test_sage_source(self):
        """Sage source."""
        result = standardize_modified_sequence("PEPTIDE[+15.9949]K", "sage")
        self.assertEqual(result, "PEPTIDE[UNIMOD:35]K")

    def test_fragpipe_source(self):
        """FragPipe source with modifications string."""
        result = standardize_modified_sequence("CPEPTIDEK", "fragpipe", "1C(57.0215)")
        self.assertEqual(result, "C[UNIMOD:4]PEPTIDEK")

    def test_fragpipe_modified_source(self):
        """FragPipe Modified Peptide column."""
        result = standardize_modified_sequence("M[147.0354]PEPTIDEK", "fragpipe_modified")
        self.assertEqual(result, "M[UNIMOD:35]PEPTIDEK")

    def test_case_insensitive_source(self):
        """Source name is case insensitive."""
        result = standardize_modified_sequence("PEPTIDE(UniMod:35)K", "DIANN")
        self.assertEqual(result, "PEPTIDE[UNIMOD:35]K")


class TestCrossEngineConsistency(unittest.TestCase):
    """Test that the same peptide from different engines standardizes to the same format."""

    def test_oxidized_met_all_engines(self):
        """Oxidized methionine should be consistent across engines."""
        diann = standardize_diann_sequence("PEPTIM(UniMod:35)EK")
        sage = standardize_sage_sequence("PEPTIM[+15.9949]EK")
        fragpipe = standardize_fragpipe_sequence("PEPTIMEK", "6M(15.9949)")

        # All should produce the same result
        self.assertEqual(diann, "PEPTIM[UNIMOD:35]EK")
        self.assertEqual(sage, "PEPTIM[UNIMOD:35]EK")
        self.assertEqual(fragpipe, "PEPTIM[UNIMOD:35]EK")

    def test_carbamidomethyl_cys_all_engines(self):
        """Carbamidomethyl cysteine should be consistent across engines."""
        diann = standardize_diann_sequence("C(UniMod:4)PEPTIDEK")
        sage = standardize_sage_sequence("C[+57.021465]PEPTIDEK")
        fragpipe = standardize_fragpipe_sequence("CPEPTIDEK", "1C(57.0215)")

        # All should produce the same result
        self.assertEqual(diann, "C[UNIMOD:4]PEPTIDEK")
        self.assertEqual(sage, "C[UNIMOD:4]PEPTIDEK")
        self.assertEqual(fragpipe, "C[UNIMOD:4]PEPTIDEK")

    def test_nterm_acetyl_all_engines(self):
        """N-terminal acetylation should be consistent across engines."""
        diann = standardize_diann_sequence("(UniMod:1)PEPTIDEK")
        sage = standardize_sage_sequence("[+42.0106]PEPTIDEK")
        fragpipe = standardize_fragpipe_sequence("PEPTIDEK", "N-term(42.0106)")

        # All should produce the same result
        self.assertEqual(diann, "[UNIMOD:1]PEPTIDEK")
        self.assertEqual(sage, "[UNIMOD:1]PEPTIDEK")
        self.assertEqual(fragpipe, "[UNIMOD:1]PEPTIDEK")

    def test_complex_peptide_all_engines(self):
        """Complex peptide with multiple mods should be consistent."""
        # N-term acetyl + Carbamidomethyl C + Oxidized M
        diann = standardize_diann_sequence("(UniMod:1)C(UniMod:4)PEPTIM(UniMod:35)EK")
        sage = standardize_sage_sequence("[+42.0106]C[+57.021]PEPTIM[+15.995]EK")
        fragpipe = standardize_fragpipe_sequence("CPEPTIMEK", "N-term(42.0106), 1C(57.0215), 7M(15.9949)")

        expected = "[UNIMOD:1]C[UNIMOD:4]PEPTIM[UNIMOD:35]EK"
        self.assertEqual(diann, expected)
        self.assertEqual(sage, expected)
        self.assertEqual(fragpipe, expected)


def run_tests():
    """Run all tests and print summary."""
    # Create test suite
    loader = unittest.TestLoader()
    suite = unittest.TestSuite()

    # Add all test classes
    suite.addTests(loader.loadTestsFromTestCase(TestMassToUnimod))
    suite.addTests(loader.loadTestsFromTestCase(TestDiannStandardization))
    suite.addTests(loader.loadTestsFromTestCase(TestSageStandardization))
    suite.addTests(loader.loadTestsFromTestCase(TestFragpipeStandardization))
    suite.addTests(loader.loadTestsFromTestCase(TestFragpipeModifiedPeptide))
    suite.addTests(loader.loadTestsFromTestCase(TestRemoveModifications))
    suite.addTests(loader.loadTestsFromTestCase(TestNormalizeSequenceIL))
    suite.addTests(loader.loadTestsFromTestCase(TestSequencesMatch))
    suite.addTests(loader.loadTestsFromTestCase(TestStandardizeModifiedSequence))
    suite.addTests(loader.loadTestsFromTestCase(TestCrossEngineConsistency))

    # Run tests with verbosity
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(suite)

    # Print summary
    print("\n" + "=" * 70)
    print(f"Tests run: {result.testsRun}")
    print(f"Failures: {len(result.failures)}")
    print(f"Errors: {len(result.errors)}")
    print(f"Success: {result.wasSuccessful()}")

    return result.wasSuccessful()


if __name__ == "__main__":
    import sys
    success = run_tests()
    sys.exit(0 if success else 1)
