"""Round-trip test for the v0.3 acyl modification registration.

The acyl masses must resolve to the correct UNIMOD ids in BOTH converters, or
the corpus's modified_sequence encoding breaks and acyl fragments get skipped at
Tier-3 (TO_V03_EXPANSION.md §2.1, Codex review).
"""
import sequence_utils as su
from fragment_matching import mass_delta_to_unimod

# (name, K-delta mass, correct UNIMOD id)  -- Lactyl is 2114 (NOT 1926)
ACYL = [("Crotonyl", 68.026215, 1363), ("Lactyl", 72.021129, 2114),
        ("Malonyl", 86.000394, 747), ("Succinyl", 100.016044, 64)]
K_RESIDUE = 128.094963


def test_acyl_delta_to_unimod():
    for name, delta, uid in ACYL:
        assert su.mass_to_unimod(delta) == f"[UNIMOD:{uid}]", name


def test_acyl_fragpipe_absolute_to_unimod():
    for name, delta, uid in ACYL:
        assert su.mass_to_unimod(K_RESIDUE + delta) == f"[UNIMOD:{uid}]", name


def test_acyl_tier3_converter():
    for name, delta, uid in ACYL:
        assert mass_delta_to_unimod(delta) == uid, name


def test_no_regression_on_existing_mods():
    assert su.mass_to_unimod(57.021464) == "[UNIMOD:4]"     # Carbamidomethyl
    assert su.mass_to_unimod(79.966331) == "[UNIMOD:21]"    # Phospho
    assert su.mass_to_unimod(15.9949) == "[UNIMOD:35]"      # Oxidation


if __name__ == "__main__":
    for fn in [test_acyl_delta_to_unimod, test_acyl_fragpipe_absolute_to_unimod,
               test_acyl_tier3_converter, test_no_regression_on_existing_mods]:
        fn()
    print("all acyl UNIMOD round-trip tests passed")
