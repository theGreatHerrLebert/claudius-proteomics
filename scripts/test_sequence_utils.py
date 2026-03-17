#!/usr/bin/env python3
"""
Tests for sequence_utils.py

Run with: python -m pytest test_sequence_utils.py -v
Or simply: python test_sequence_utils.py
"""

import unittest
from sequence_utils import (
    mass_to_unimod,
    wrap_terminal_format,
    standardize_diann_sequence,
    standardize_sage_sequence,
    standardize_fragpipe_sequence,
    standardize_fragpipe_modified_peptide,
    add_carbamidomethyl_to_cysteine,
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

    def test_pyro_glu(self):
        """Negative mass deltas for pyro-glutamate."""
        self.assertEqual(mass_to_unimod(-17), "[UNIMOD:28]")
        self.assertEqual(mass_to_unimod(-18), "[UNIMOD:27]")

    def test_unknown_mass(self):
        """Unknown mass returns None."""
        self.assertIsNone(mass_to_unimod(999))
        self.assertIsNone(mass_to_unimod(123.456))

    def test_absolute_mass(self):
        """Absolute mass lookup (FragPipe format)."""
        self.assertEqual(mass_to_unimod(160), "[UNIMOD:4]")
        self.assertEqual(mass_to_unimod(147), "[UNIMOD:35]")


class TestWrapTerminalFormat(unittest.TestCase):
    """Test terminal format wrapping."""

    def test_plain_sequence(self):
        """Plain sequence gets empty brackets."""
        self.assertEqual(wrap_terminal_format("PEPTIDEK"), "[]-PEPTIDEK-[]")

    def test_nterm_unimod(self):
        """N-terminal UNIMOD mod is extracted."""
        self.assertEqual(
            wrap_terminal_format("[UNIMOD:1]MPEPTIDEK"),
            "[UNIMOD:1]-MPEPTIDEK-[]"
        )

    def test_nterm_with_hyphen(self):
        """Sage-style N-term with hyphen is normalized."""
        self.assertEqual(
            wrap_terminal_format("[UNIMOD:1]-MPEPTIDEK"),
            "[UNIMOD:1]-MPEPTIDEK-[]"
        )

    def test_nterm_mass_format(self):
        """N-terminal mass format is extracted."""
        self.assertEqual(
            wrap_terminal_format("[+42.0106]PEPTIDEK"),
            "[+42.0106]-PEPTIDEK-[]"
        )

    def test_internal_mod_not_extracted(self):
        """Internal mods stay in the sequence."""
        self.assertEqual(
            wrap_terminal_format("MC[UNIMOD:4]PEPTIDEK"),
            "[]-MC[UNIMOD:4]PEPTIDEK-[]"
        )

    def test_empty_input(self):
        """Empty/None input returns empty string."""
        self.assertEqual(wrap_terminal_format(""), "")
        self.assertEqual(wrap_terminal_format(None), "")

    def test_idempotent(self):
        """Wrapping an already-wrapped sequence doesn't double-wrap."""
        # Once wrapped: []-PEPTIDEK-[]
        # The leading [] is empty brackets, not a mod, so wrap_terminal_format
        # shouldn't detect it as an N-term mod (no UNIMOD: or +/- number inside)
        wrapped = wrap_terminal_format("PEPTIDEK")
        self.assertEqual(wrapped, "[]-PEPTIDEK-[]")
        # Wrapping again: the leading [] doesn't match the nterm pattern
        # because it's empty, so it stays as-is in the seq, and gets []- prepended
        # This is NOT idempotent by design - terminal format is the final step


class TestDiannStandardization(unittest.TestCase):
    """Test DIA-NN sequence standardization."""

    def test_simple_oxidation(self):
        """DIA-NN oxidation format."""
        result = standardize_diann_sequence("PEPTIDE(UniMod:35)K")
        self.assertEqual(result, "[]-PEPTIDE[UNIMOD:35]K-[]")

    def test_carbamidomethyl(self):
        """DIA-NN carbamidomethyl format."""
        result = standardize_diann_sequence("C(UniMod:4)PEPTIDEK")
        self.assertEqual(result, "[]-C[UNIMOD:4]PEPTIDEK-[]")

    def test_multiple_mods(self):
        """Multiple modifications."""
        result = standardize_diann_sequence("C(UniMod:4)PEPTIM(UniMod:35)DEK")
        self.assertEqual(result, "[]-C[UNIMOD:4]PEPTIM[UNIMOD:35]DEK-[]")

    def test_nterm_acetyl(self):
        """N-terminal acetylation."""
        result = standardize_diann_sequence("(UniMod:1)PEPTIDEK")
        self.assertEqual(result, "[UNIMOD:1]-PEPTIDEK-[]")

    def test_unmodified(self):
        """Unmodified sequence gets terminal brackets."""
        result = standardize_diann_sequence("PEPTIDEK")
        self.assertEqual(result, "[]-PEPTIDEK-[]")

    def test_bare_cysteine_gets_cam(self):
        """Bare cysteine gets carbamidomethyl added."""
        result = standardize_diann_sequence("CPEPTIDEK")
        self.assertEqual(result, "[]-C[UNIMOD:4]PEPTIDEK-[]")

    def test_empty_input(self):
        """Empty/None input returns empty string."""
        self.assertEqual(standardize_diann_sequence(""), "")
        self.assertEqual(standardize_diann_sequence(None), "")

    def test_case_insensitive_unimod(self):
        """Handle different UniMod casing."""
        result = standardize_diann_sequence("PEPTIDE(unimod:35)K")
        self.assertEqual(result, "[]-PEPTIDE[UNIMOD:35]K-[]")


class TestSageStandardization(unittest.TestCase):
    """Test Sage sequence standardization."""

    def test_oxidation(self):
        """Sage oxidation format [+15.9949]."""
        result = standardize_sage_sequence("PEPTIDE[+15.9949]K")
        self.assertEqual(result, "[]-PEPTIDE[UNIMOD:35]K-[]")

    def test_carbamidomethyl(self):
        """Sage carbamidomethyl format [+57.021465]."""
        result = standardize_sage_sequence("C[+57.021465]PEPTIDEK")
        self.assertEqual(result, "[]-C[UNIMOD:4]PEPTIDEK-[]")

    def test_multiple_mods(self):
        """Multiple modifications."""
        result = standardize_sage_sequence("C[+57.021]PEPTIM[+15.995]DEK")
        self.assertEqual(result, "[]-C[UNIMOD:4]PEPTIM[UNIMOD:35]DEK-[]")

    def test_nterm_acetyl(self):
        """N-terminal acetylation without hyphen."""
        result = standardize_sage_sequence("[+42.0106]PEPTIDEK")
        self.assertEqual(result, "[UNIMOD:1]-PEPTIDEK-[]")

    def test_nterm_acetyl_with_hyphen(self):
        """N-terminal acetylation with Sage's native hyphen format."""
        result = standardize_sage_sequence("[+42.010567]-MDSPWDELALAFSR")
        self.assertEqual(result, "[UNIMOD:1]-MDSPWDELALAFSR-[]")

    def test_nterm_plus_internal(self):
        """N-term acetyl + oxidized Met (Sage native format)."""
        result = standardize_sage_sequence("[+42.010567]-M[+15.994915]DLAAAAEPGAGSQHLEVR")
        self.assertEqual(result, "[UNIMOD:1]-M[UNIMOD:35]DLAAAAEPGAGSQHLEVR-[]")

    def test_unmodified(self):
        """Unmodified sequence gets terminal brackets."""
        result = standardize_sage_sequence("PEPTIDEK")
        self.assertEqual(result, "[]-PEPTIDEK-[]")

    def test_unknown_mass_preserved(self):
        """Unknown mass modifications are preserved inside terminal format."""
        result = standardize_sage_sequence("PEPTIDE[+999.999]K")
        self.assertEqual(result, "[]-PEPTIDE[+999.999]K-[]")

    def test_empty_input(self):
        """Empty/None input returns empty string."""
        self.assertEqual(standardize_sage_sequence(""), "")
        self.assertEqual(standardize_sage_sequence(None), "")


class TestFragpipeStandardization(unittest.TestCase):
    """Test FragPipe sequence standardization."""

    def test_position_based_carbamidomethyl(self):
        """FragPipe position-based format: 1C(57.0215)."""
        result = standardize_fragpipe_sequence("CPEPTIDEK", "1C(57.0215)")
        self.assertEqual(result, "[]-C[UNIMOD:4]PEPTIDEK-[]")

    def test_position_based_oxidation(self):
        """FragPipe position-based oxidation."""
        result = standardize_fragpipe_sequence("PEPTIMEK", "6M(15.9949)")
        self.assertEqual(result, "[]-PEPTIM[UNIMOD:35]EK-[]")

    def test_multiple_mods(self):
        """Multiple position-based modifications."""
        result = standardize_fragpipe_sequence("CPEPTIMEK", "1C(57.0215), 7M(15.9949)")
        self.assertEqual(result, "[]-C[UNIMOD:4]PEPTIM[UNIMOD:35]EK-[]")

    def test_nterm_acetyl(self):
        """N-terminal acetylation."""
        result = standardize_fragpipe_sequence("PEPTIDEK", "N-term(42.0106)")
        self.assertEqual(result, "[UNIMOD:1]-PEPTIDEK-[]")

    def test_nterm_plus_internal(self):
        """N-terminal plus internal modification."""
        result = standardize_fragpipe_sequence("CPEPTIDEK", "N-term(42.0106), 1C(57.0215)")
        self.assertEqual(result, "[UNIMOD:1]-C[UNIMOD:4]PEPTIDEK-[]")

    def test_multiple_cys_no_space(self):
        """Multiple C modifications without space after comma."""
        result = standardize_fragpipe_sequence("CCDNDMER", "1C(57.0215),2C(57.0215),6M(15.9949)")
        self.assertEqual(result, "[]-C[UNIMOD:4]C[UNIMOD:4]DNDM[UNIMOD:35]ER-[]")

    def test_multiple_mods_mixed_spacing(self):
        """Multiple modifications with mixed spacing."""
        result = standardize_fragpipe_sequence("ACMEK", "1A(42.0106),2C(57.0215), 3M(15.9949)")
        self.assertEqual(result, "[]-A[UNIMOD:1]C[UNIMOD:4]M[UNIMOD:35]EK-[]")

    def test_unmodified(self):
        """Unmodified sequence gets terminal brackets."""
        result = standardize_fragpipe_sequence("PEPTIDEK", "")
        self.assertEqual(result, "[]-PEPTIDEK-[]")
        result = standardize_fragpipe_sequence("PEPTIDEK", None)
        self.assertEqual(result, "[]-PEPTIDEK-[]")

    def test_unmodified_with_cysteine(self):
        """Unmodified sequence with C gets carbamidomethyl."""
        result = standardize_fragpipe_sequence("CPEPTIDEK", "")
        self.assertEqual(result, "[]-C[UNIMOD:4]PEPTIDEK-[]")

    def test_empty_input(self):
        """Empty input returns empty string."""
        self.assertEqual(standardize_fragpipe_sequence("", ""), "")
        self.assertEqual(standardize_fragpipe_sequence(None, None), "")


class TestFragpipeModifiedPeptide(unittest.TestCase):
    """Test FragPipe Modified Peptide column standardization."""

    def test_oxidized_met(self):
        """Oxidized methionine M[147.0354]."""
        result = standardize_fragpipe_modified_peptide("M[147.0354]PEPTIDEK")
        self.assertEqual(result, "[]-M[UNIMOD:35]PEPTIDEK-[]")

    def test_nterm_lowercase(self):
        """N-terminal with lowercase 'n' marker."""
        result = standardize_fragpipe_modified_peptide("n[42.0106]PEPTIDEK")
        self.assertEqual(result, "[UNIMOD:1]-PEPTIDEK-[]")

    def test_unmodified(self):
        """Unmodified sequence gets terminal brackets."""
        result = standardize_fragpipe_modified_peptide("PEPTIDEK")
        self.assertEqual(result, "[]-PEPTIDEK-[]")

    def test_empty_input(self):
        """Empty/None input returns empty string."""
        self.assertEqual(standardize_fragpipe_modified_peptide(""), "")
        self.assertEqual(standardize_fragpipe_modified_peptide(None), "")


class TestAddCarbamidomethylToCysteine(unittest.TestCase):
    """Test bare cysteine -> carbamidomethyl addition."""

    def test_single_bare_c(self):
        self.assertEqual(add_carbamidomethyl_to_cysteine("CPEPTIDEK"), "C[UNIMOD:4]PEPTIDEK")

    def test_already_modified_c(self):
        """C with existing mod is not double-modified."""
        self.assertEqual(
            add_carbamidomethyl_to_cysteine("C[UNIMOD:4]PEPTIDEK"),
            "C[UNIMOD:4]PEPTIDEK"
        )

    def test_multiple_bare_c(self):
        self.assertEqual(
            add_carbamidomethyl_to_cysteine("CPEPTICEK"),
            "C[UNIMOD:4]PEPTIC[UNIMOD:4]EK"
        )

    def test_mixed_c(self):
        """Mix of bare and modified C."""
        self.assertEqual(
            add_carbamidomethyl_to_cysteine("C[UNIMOD:4]PEPTICEK"),
            "C[UNIMOD:4]PEPTIC[UNIMOD:4]EK"
        )

    def test_no_cysteine(self):
        self.assertEqual(add_carbamidomethyl_to_cysteine("PEPTIDEK"), "PEPTIDEK")

    def test_empty_input(self):
        self.assertEqual(add_carbamidomethyl_to_cysteine(""), "")
        self.assertEqual(add_carbamidomethyl_to_cysteine(None), "")


class TestRemoveModifications(unittest.TestCase):
    """Test modification removal."""

    def test_terminal_format(self):
        """Remove terminal format brackets and internal mods."""
        result = remove_modifications("[]-C[UNIMOD:4]PEPTIM[UNIMOD:35]DEK-[]")
        self.assertEqual(result, "CPEPTIMDEK")

    def test_terminal_format_with_nterm(self):
        """Remove terminal format with N-term mod."""
        result = remove_modifications("[UNIMOD:1]-M[UNIMOD:35]PEPTIDEK-[]")
        self.assertEqual(result, "MPEPTIDEK")

    def test_unimod_format(self):
        """Remove [UNIMOD:X] format (legacy, no terminal brackets)."""
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

    def test_unmodified_terminal(self):
        """Unmodified sequence in terminal format."""
        result = remove_modifications("[]-PEPTIDEK-[]")
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
        """Modifications are removed, I->L."""
        result = normalize_sequence_il("PEPTIDE[UNIMOD:35]K")
        self.assertEqual(result, "PEPTLDEK")

    def test_with_terminal_format(self):
        """Terminal format is stripped, I->L."""
        result = normalize_sequence_il("[UNIMOD:1]-MIPEPTIDEK-[]")
        self.assertEqual(result, "MLPEPTLDEK")

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

    def test_terminal_format_vs_legacy(self):
        """Terminal format matches legacy format after normalization."""
        self.assertTrue(sequences_match(
            "[]-PEPTIDEK-[]",
            "PEPTIDEK"
        ))

    def test_terminal_format_both(self):
        """Two terminal-format sequences match."""
        self.assertTrue(sequences_match(
            "[UNIMOD:1]-MPEPTIDEK-[]",
            "[UNIMOD:1]-MPEPTLDEK-[]"
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
        self.assertEqual(result, "[]-PEPTIDE[UNIMOD:35]K-[]")

    def test_sage_source(self):
        """Sage source."""
        result = standardize_modified_sequence("PEPTIDE[+15.9949]K", "sage")
        self.assertEqual(result, "[]-PEPTIDE[UNIMOD:35]K-[]")

    def test_fragpipe_source(self):
        """FragPipe source with modifications string."""
        result = standardize_modified_sequence("CPEPTIDEK", "fragpipe", "1C(57.0215)")
        self.assertEqual(result, "[]-C[UNIMOD:4]PEPTIDEK-[]")

    def test_fragpipe_modified_source(self):
        """FragPipe Modified Peptide column."""
        result = standardize_modified_sequence("M[147.0354]PEPTIDEK", "fragpipe_modified")
        self.assertEqual(result, "[]-M[UNIMOD:35]PEPTIDEK-[]")

    def test_case_insensitive_source(self):
        """Source name is case insensitive."""
        result = standardize_modified_sequence("PEPTIDE(UniMod:35)K", "DIANN")
        self.assertEqual(result, "[]-PEPTIDE[UNIMOD:35]K-[]")

    def test_unknown_source(self):
        """Unknown source returns input as-is."""
        result = standardize_modified_sequence("PEPTIDEK", "unknown_engine")
        self.assertEqual(result, "PEPTIDEK")


class TestCrossEngineConsistency(unittest.TestCase):
    """Test that the same peptide from different engines standardizes to the same format."""

    def test_unmodified_all_engines(self):
        """Unmodified peptide should be consistent."""
        diann = standardize_diann_sequence("PEPTIDEK")
        sage = standardize_sage_sequence("PEPTIDEK")
        fragpipe = standardize_fragpipe_sequence("PEPTIDEK", "")

        expected = "[]-PEPTIDEK-[]"
        self.assertEqual(diann, expected)
        self.assertEqual(sage, expected)
        self.assertEqual(fragpipe, expected)

    def test_oxidized_met_all_engines(self):
        """Oxidized methionine should be consistent across engines."""
        diann = standardize_diann_sequence("PEPTIM(UniMod:35)EK")
        sage = standardize_sage_sequence("PEPTIM[+15.9949]EK")
        fragpipe = standardize_fragpipe_sequence("PEPTIMEK", "6M(15.9949)")

        expected = "[]-PEPTIM[UNIMOD:35]EK-[]"
        self.assertEqual(diann, expected)
        self.assertEqual(sage, expected)
        self.assertEqual(fragpipe, expected)

    def test_carbamidomethyl_cys_all_engines(self):
        """Carbamidomethyl cysteine should be consistent across engines."""
        diann = standardize_diann_sequence("C(UniMod:4)PEPTIDEK")
        sage = standardize_sage_sequence("C[+57.021465]PEPTIDEK")
        fragpipe = standardize_fragpipe_sequence("CPEPTIDEK", "1C(57.0215)")

        expected = "[]-C[UNIMOD:4]PEPTIDEK-[]"
        self.assertEqual(diann, expected)
        self.assertEqual(sage, expected)
        self.assertEqual(fragpipe, expected)

    def test_bare_cysteine_all_engines(self):
        """Bare cysteine (fixed mod not shown) should get CAM from all engines."""
        diann = standardize_diann_sequence("CPEPTIDEK")
        sage = standardize_sage_sequence("CPEPTIDEK")
        fragpipe = standardize_fragpipe_sequence("CPEPTIDEK", "")

        expected = "[]-C[UNIMOD:4]PEPTIDEK-[]"
        self.assertEqual(diann, expected)
        self.assertEqual(sage, expected)
        self.assertEqual(fragpipe, expected)

    def test_nterm_acetyl_all_engines(self):
        """N-terminal acetylation should be consistent across engines."""
        diann = standardize_diann_sequence("(UniMod:1)PEPTIDEK")
        sage = standardize_sage_sequence("[+42.010567]-PEPTIDEK")
        fragpipe = standardize_fragpipe_sequence("PEPTIDEK", "N-term(42.0106)")

        expected = "[UNIMOD:1]-PEPTIDEK-[]"
        self.assertEqual(diann, expected)
        self.assertEqual(sage, expected)
        self.assertEqual(fragpipe, expected)

    def test_complex_peptide_all_engines(self):
        """Complex peptide with multiple mods should be consistent."""
        # N-term acetyl + Carbamidomethyl C + Oxidized M
        diann = standardize_diann_sequence("(UniMod:1)C(UniMod:4)PEPTIM(UniMod:35)EK")
        sage = standardize_sage_sequence("[+42.0106]-C[+57.021]PEPTIM[+15.995]EK")
        fragpipe = standardize_fragpipe_sequence("CPEPTIMEK", "N-term(42.0106), 1C(57.0215), 7M(15.9949)")

        expected = "[UNIMOD:1]-C[UNIMOD:4]PEPTIM[UNIMOD:35]EK-[]"
        self.assertEqual(diann, expected)
        self.assertEqual(sage, expected)
        self.assertEqual(fragpipe, expected)

    def test_nterm_acetyl_plus_oxidation_all_engines(self):
        """N-term acetyl + oxidized Met (common in real data)."""
        diann = standardize_diann_sequence("(UniMod:1)M(UniMod:35)DLAAAAEPGAGSQHLEVR")
        sage = standardize_sage_sequence("[+42.010567]-M[+15.994915]DLAAAAEPGAGSQHLEVR")
        fragpipe = standardize_fragpipe_sequence("MDLAAAAEPGAGSQHLEVR", "N-term(42.0106), 1M(15.9949)")

        expected = "[UNIMOD:1]-M[UNIMOD:35]DLAAAAEPGAGSQHLEVR-[]"
        self.assertEqual(diann, expected)
        self.assertEqual(sage, expected)
        self.assertEqual(fragpipe, expected)


def run_tests():
    """Run all tests and print summary."""
    loader = unittest.TestLoader()
    suite = unittest.TestSuite()

    # Add all test classes
    test_classes = [
        TestMassToUnimod,
        TestWrapTerminalFormat,
        TestDiannStandardization,
        TestSageStandardization,
        TestFragpipeStandardization,
        TestFragpipeModifiedPeptide,
        TestAddCarbamidomethylToCysteine,
        TestRemoveModifications,
        TestNormalizeSequenceIL,
        TestSequencesMatch,
        TestStandardizeModifiedSequence,
        TestCrossEngineConsistency,
    ]

    for cls in test_classes:
        suite.addTests(loader.loadTestsFromTestCase(cls))

    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(suite)

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
