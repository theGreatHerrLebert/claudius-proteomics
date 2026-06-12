#!/usr/bin/env python3
"""
Sequence Standardization Utilities

Converts modified peptide sequences from different search engines to a unified
UNIMOD terminal format:

    [N-term]-SEQUENCE[UNIMOD:X]-[C-term]

Examples:
    []-PEPTIDEK-[]                          # Unmodified
    []-AAC[UNIMOD:4]PEPTIDEK-[]             # Carbamidomethyl (C)
    []-M[UNIMOD:35]PEPTIDEK-[]              # Oxidation (M)
    [UNIMOD:1]-PEPTIDEK-[]                  # N-term Acetyl
    [UNIMOD:1]-M[UNIMOD:35]C[UNIMOD:4]K-[] # N-term + internal mods

Supported input formats:
- FragPipe: M[147.0354] (absolute mass) or 4C(57.0215) (position + delta)
- DIA-NN: (UniMod:35)
- Sage: [+15.9949] (mass delta), [+42.010567]-SEQ (N-term with hyphen)
- MaxQuant: _(ac)PEPTIDE_ with (ox), (ph), etc.
"""

import re
from typing import Optional, Dict, Tuple
import numpy as np


# Mass to UNIMOD mapping (rounded integer mass delta -> UNIMOD ID)
# Based on common proteomics modifications
MASS_TO_UNIMOD: Dict[int, Tuple[int, str]] = {
    # Cysteine modifications
    57: (4, "Carbamidomethyl"),      # C + 57.021464
    58: (4, "Carbamidomethyl"),      # Rounding tolerance

    # Oxidation
    16: (35, "Oxidation"),            # M + 15.9949
    15: (35, "Oxidation"),            # Rounding tolerance

    # N-terminal acetylation
    42: (1, "Acetyl"),                # N-term + 42.0106
    43: (1, "Acetyl"),                # Rounding tolerance

    # Phosphorylation
    80: (21, "Phospho"),              # S/T/Y + 79.9663
    79: (21, "Phospho"),              # Rounding tolerance

    # Methylation
    14: (34, "Methyl"),               # K/R + 14.0157

    # Dimethylation
    28: (36, "Dimethyl"),             # K/R + 28.0314

    # Deamidation
    1: (7, "Deamidated"),             # N/Q + 0.9840

    # Pyro-glutamate (N-terminal)
    -17: (28, "Gln->pyro-Glu"),       # Q - 17.0265 (N-term Q cyclization)
    -18: (27, "Glu->pyro-Glu"),       # E - 18.0106 (N-term E water loss)

    # Ammonia loss
    -16: (385, "Ammonia-loss"),       # -17 rounded differently

    # TMT labels (approximate)
    229: (737, "TMT6plex"),           # TMT6plex + 229.1629
    304: (2016, "TMTpro"),            # TMTpro + 304.2071

    # Lysine acylations (v0.3 expansion) — distinct integers, no collisions
    68: (1363, "Crotonyl"),           # K + 68.026215
    72: (2114, "Lactyl"),             # K + 72.021129  (2114 = L-lactyl; NOT 1926, which is a methylglyoxal artefact)
    86: (747, "Malonyl"),             # K + 86.000394
    100: (64, "Succinyl"),            # K + 100.016044
}

# Absolute mass to UNIMOD (for FragPipe format)
ABSOLUTE_MASS_TO_UNIMOD: Dict[int, Tuple[int, str]] = {
    # Carbamidomethyl-Cysteine: 103.009 (C) + 57.021 = 160.030
    160: (4, "Carbamidomethyl"),

    # Oxidized Methionine: 131.040 (M) + 15.995 = 147.035
    147: (35, "Oxidation"),

    # Lysine acylations (K residue 128.09496 + delta) for the FragPipe absolute format
    196: (1363, "Crotonyl"),          # K 128.095 + 68.026 = 196.121
    200: (2114, "Lactyl"),            # K 128.095 + 72.021 = 200.116
    214: (747, "Malonyl"),            # K 128.095 + 86.000 = 214.095
    228: (64, "Succinyl"),            # K 128.095 + 100.016 = 228.111
}


def mass_to_unimod(mass: float) -> Optional[str]:
    """Convert modification mass (delta) to UNIMOD annotation.

    Args:
        mass: Modification mass delta (e.g., 57.021 for carbamidomethyl)

    Returns:
        UNIMOD string like "[UNIMOD:4]" or None if not found
    """
    rounded = int(round(mass))

    if rounded in MASS_TO_UNIMOD:
        unimod_id, _ = MASS_TO_UNIMOD[rounded]
        return f"[UNIMOD:{unimod_id}]"

    # Try absolute mass lookup
    if rounded in ABSOLUTE_MASS_TO_UNIMOD:
        unimod_id, _ = ABSOLUTE_MASS_TO_UNIMOD[rounded]
        return f"[UNIMOD:{unimod_id}]"

    return None


def wrap_terminal_format(sequence: str) -> str:
    """Wrap a standardized sequence in [nterm]-SEQ-[cterm] format.

    Detects N-terminal modifications (mod before the first amino acid letter)
    and separates them with a hyphen. C-terminal bracket is always [].

    Examples:
        PEPTIDEK                    -> []-PEPTIDEK-[]
        [UNIMOD:1]MPEPTIDEK         -> [UNIMOD:1]-MPEPTIDEK-[]
        [UNIMOD:1]-MPEPTIDEK        -> [UNIMOD:1]-MPEPTIDEK-[]  (Sage already has hyphen)
        MC[UNIMOD:4]PEPTIDEK        -> []-MC[UNIMOD:4]PEPTIDEK-[]
    """
    if not sequence or not isinstance(sequence, str):
        return ""

    nterm = "[]"
    seq = sequence

    # Detect leading modification before first amino acid letter.
    # Matches [UNIMOD:X] or [+/-X.XXX] or [X.XXX] at start, optionally followed by -
    nterm_match = re.match(r'^(\[(?:UNIMOD:\d+|[+-]?\d+\.?\d*)\])-?', seq)
    if nterm_match:
        nterm = nterm_match.group(1)
        seq = seq[nterm_match.end():]

    return f"{nterm}-{seq}-[]"


def standardize_diann_sequence(modified_sequence: str) -> str:
    """Convert DIA-NN modified sequence to standard UNIMOD format.

    DIA-NN format: PEPTIDE(UniMod:35)K
    Standard format: PEPTIDE[UNIMOD:35]K

    Also adds carbamidomethyl [UNIMOD:4] to bare cysteines (fixed mod).

    Args:
        modified_sequence: DIA-NN modified sequence

    Returns:
        Standardized sequence with [UNIMOD:X] format
    """
    if not modified_sequence or not isinstance(modified_sequence, str):
        return ""

    # Replace (UniMod:X) with [UNIMOD:X]
    result = modified_sequence.replace("(", "[").replace(")", "]")
    # Handle case-insensitive UniMod -> UNIMOD
    result = re.sub(r'unimod', 'UNIMOD', result, flags=re.IGNORECASE)

    # Add carbamidomethyl to bare cysteines (fixed mod not shown by DIA-NN)
    result = add_carbamidomethyl_to_cysteine(result)

    return wrap_terminal_format(result)


def standardize_sage_sequence(modified_sequence: str) -> str:
    """Convert Sage modified sequence to standard UNIMOD format.

    Sage format: PEPTIDE[+15.9949]K or PEPTIDE[+57.021465]K
    Standard format: PEPTIDE[UNIMOD:35]K or PEPTIDE[UNIMOD:4]K

    Also ensures any bare cysteines get carbamidomethyl (edge case).

    Args:
        modified_sequence: Sage modified sequence

    Returns:
        Standardized sequence with [UNIMOD:X] format
    """
    if not modified_sequence or not isinstance(modified_sequence, str):
        return ""

    # Find all mass modifications: [+X.XXX] or [-X.XXX]
    pattern = r'\[([+-]?\d+\.?\d*)\]'

    def replace_mass(match):
        mass_str = match.group(1)
        try:
            mass = float(mass_str)
            # First try the exact mass (preserving sign for negative mods like pyro-Glu)
            unimod = mass_to_unimod(mass)
            if unimod:
                return unimod
            # Then try absolute value for positive-only lookups
            unimod = mass_to_unimod(abs(mass))
            if unimod:
                return unimod
            # Keep original if no mapping found
            return match.group(0)
        except ValueError:
            return match.group(0)

    result = re.sub(pattern, replace_mass, modified_sequence)

    # Add carbamidomethyl to any bare cysteines (edge case - Sage usually includes it)
    result = add_carbamidomethyl_to_cysteine(result)

    return wrap_terminal_format(result)


def standardize_fragpipe_sequence(peptide: str, modifications: str) -> str:
    """Convert FragPipe modified sequence to standard UNIMOD format.

    FragPipe format:
        peptide = "PEPTIDEK"
        modifications = "4C(57.0215), 7M(15.9949)" or "N-term(42.0106)"

    Standard format: PEPT[UNIMOD:4]IDE[UNIMOD:35]K

    Args:
        peptide: Plain peptide sequence
        modifications: Modification string from FragPipe

    Returns:
        Standardized sequence with [UNIMOD:X] format
    """
    if not peptide or not isinstance(peptide, str):
        return ""

    if not modifications or not isinstance(modifications, str) or modifications == "":
        seq = add_carbamidomethyl_to_cysteine(peptide)
        return f"[]-{seq}-[]"

    # Parse modifications into position -> UNIMOD mapping
    mods_by_position = {}
    n_term_mod = None

    # Pattern for position-based mods: "4C(57.0215)"
    pos_pattern = r"(\d+)([A-Za-z])\(([\d.]+)\)"

    try:
        # Split on comma (with or without space)
        mod_list = re.split(r',\s*', modifications)

        for mod in mod_list:
            mod = mod.strip()

            # N-terminal modification
            if mod.startswith("N-term"):
                mass_match = re.search(r"\(([\d.]+)\)", mod)
                if mass_match:
                    mass = float(mass_match.group(1))
                    unimod = mass_to_unimod(mass)
                    if unimod:
                        n_term_mod = unimod
            else:
                # Position-based modification
                match = re.match(pos_pattern, mod)
                if match:
                    pos = int(match.group(1)) - 1  # Convert to 0-indexed
                    aa = match.group(2)
                    mass = float(match.group(3))
                    unimod = mass_to_unimod(mass)
                    if unimod:
                        mods_by_position[pos] = (aa, unimod)
    except Exception:
        # If parsing fails, return wrapped plain sequence
        seq = add_carbamidomethyl_to_cysteine(peptide)
        return f"[]-{seq}-[]"

    # Build modified sequence (without N-term, that goes in the terminal wrapper)
    result = []

    for i, aa in enumerate(peptide):
        result.append(aa)
        if i in mods_by_position:
            _, unimod = mods_by_position[i]
            result.append(unimod)

    seq = "".join(result)
    seq = add_carbamidomethyl_to_cysteine(seq)

    # Wrap with terminal format: N-term mod goes in the bracket
    nterm = n_term_mod if n_term_mod else "[]"
    return f"{nterm}-{seq}-[]"


def standardize_fragpipe_modified_peptide(modified_peptide: str) -> str:
    """Convert FragPipe 'Modified Peptide' column to standard UNIMOD format.

    FragPipe Modified Peptide format: M[147.0354]PEPTIDEK or n[42.0106]PEPTIDEK
    Standard format: M[UNIMOD:35]PEPTIDEK or [UNIMOD:1]PEPTIDEK

    Also adds carbamidomethyl [UNIMOD:4] to bare cysteines (fixed mod).

    Args:
        modified_peptide: FragPipe Modified Peptide string

    Returns:
        Standardized sequence with [UNIMOD:X] format
    """
    if not modified_peptide or not isinstance(modified_peptide, str):
        return ""

    # Handle n-terminal lowercase 'n' marker
    result = modified_peptide
    if result.startswith("n["):
        result = result[1:]  # Remove 'n' prefix

    # Find all mass modifications: X[XXX.XXXX] or X[-XXX.XXXX]
    pattern = r'\[([+-]?\d+\.?\d*)\]'

    def replace_mass(match):
        mass_str = match.group(1)
        try:
            mass = float(mass_str)
            rounded = int(round(mass))

            # First try absolute mass lookup (for FragPipe format like [147.0354])
            if rounded in ABSOLUTE_MASS_TO_UNIMOD:
                unimod_id, _ = ABSOLUTE_MASS_TO_UNIMOD[rounded]
                return f"[UNIMOD:{unimod_id}]"

            # Then try delta mass lookup (for masses like [57.0215] or [-17.0265])
            if rounded in MASS_TO_UNIMOD:
                unimod_id, _ = MASS_TO_UNIMOD[rounded]
                return f"[UNIMOD:{unimod_id}]"

            # Keep original if no mapping found
            return match.group(0)
        except ValueError:
            return match.group(0)

    result = re.sub(pattern, replace_mass, result)

    # Add carbamidomethyl to bare cysteines (fixed mod not shown by FragPipe)
    result = add_carbamidomethyl_to_cysteine(result)

    return wrap_terminal_format(result)


def add_carbamidomethyl_to_cysteine(sequence: str) -> str:
    """Add carbamidomethyl modification to all bare cysteines.

    Many search engines (FragPipe, DIA-NN) treat carbamidomethyl as a fixed
    modification and don't include it in the modified sequence output.
    This function adds [UNIMOD:4] to all C that don't already have a modification.

    Args:
        sequence: Modified sequence (may have some C modified, some bare)

    Returns:
        Sequence with all C having [UNIMOD:4]
    """
    if not sequence or not isinstance(sequence, str):
        return ""

    # Find all C that are NOT followed by [ (i.e., bare C)
    # Replace with C[UNIMOD:4]
    result = re.sub(r'C(?!\[)', 'C[UNIMOD:4]', sequence)
    return result


def remove_modifications(sequence: str) -> str:
    """Remove all modification annotations from sequence.

    Handles:
    - Terminal format: [nterm]-SEQ-[cterm]
    - [UNIMOD:X]
    - [+X.XXX]
    - (UniMod:X)
    - [XXX.XXX] (absolute mass)

    Args:
        sequence: Modified sequence (with or without terminal format)

    Returns:
        Plain sequence without modifications
    """
    if not sequence or not isinstance(sequence, str):
        return ""

    result = sequence

    # Strip terminal format: [nterm]-SEQ-[cterm]
    terminal_match = re.match(r'^\[([^\]]*)\]-(.*)-\[([^\]]*)\]$', result)
    if terminal_match:
        result = terminal_match.group(2)  # Extract just the sequence part

    # Remove [UNIMOD:X]
    result = re.sub(r'\[UNIMOD:\d+\]', '', result)
    # Remove [+/-X.XXX]
    result = re.sub(r'\[[+-]?\d+\.?\d*\]', '', result)
    # Remove (UniMod:X)
    result = re.sub(r'\(UniMod:\d+\)', '', result, flags=re.IGNORECASE)

    return result


def normalize_sequence_il(sequence: str) -> str:
    """Normalize sequence by replacing I with L (isoleucine/leucine are isobaric).

    Also removes modifications and converts to uppercase.

    Args:
        sequence: Peptide sequence (may contain modifications)

    Returns:
        Normalized plain sequence
    """
    plain = remove_modifications(sequence)
    return plain.upper().replace("I", "L")


def sequences_match(seq1: str, seq2: str) -> bool:
    """Check if two sequences match after I/L normalization.

    Args:
        seq1: First sequence (may contain modifications)
        seq2: Second sequence (may contain modifications)

    Returns:
        True if sequences match after normalization
    """
    return normalize_sequence_il(seq1) == normalize_sequence_il(seq2)


# Convenience function that handles any format
def standardize_modified_sequence(
    sequence: str,
    source: str,
    modifications: Optional[str] = None
) -> str:
    """Standardize modified sequence from any search engine.

    Args:
        sequence: Modified or plain sequence
        source: Source engine: "diann", "sage", "fragpipe", "fragpipe_modified"
        modifications: For FragPipe, the "Assigned Modifications" column

    Returns:
        Standardized sequence with [UNIMOD:X] format
    """
    source = source.lower()

    if source == "diann":
        return standardize_diann_sequence(sequence)
    elif source == "sage":
        return standardize_sage_sequence(sequence)
    elif source == "fragpipe" and modifications:
        return standardize_fragpipe_sequence(sequence, modifications)
    elif source == "fragpipe_modified":
        return standardize_fragpipe_modified_peptide(sequence)
    else:
        return sequence


if __name__ == "__main__":
    # Test cases
    print("Testing sequence standardization:\n")

    # DIA-NN
    diann_tests = [
        "PEPTIDE(UniMod:35)K",
        "(UniMod:1)MPEPTIDEK",
        "AC(UniMod:4)PEPTIDEK",
        "PEPTIDEK",
    ]
    print("DIA-NN:")
    for seq in diann_tests:
        print(f"  {seq:40s} -> {standardize_diann_sequence(seq)}")

    # Sage
    sage_tests = [
        "PEPTIDE[+15.9949]K",
        "C[+57.021465]PEPTIDEK",
        "[+42.010567]-MDSPWDELALAFSR",
        "[+42.010567]-M[+15.994915]DLAAAAEPGAGSQHLEVR",
        "PEPTIDEK",
    ]
    print("\nSage:")
    for seq in sage_tests:
        print(f"  {seq:50s} -> {standardize_sage_sequence(seq)}")

    # FragPipe with modifications string
    fp_tests = [
        ("CPEPTIDEK", "1C(57.0215)"),
        ("MPEPTIDEK", "N-term(42.0106)"),
        ("MPEPTIDEK", "N-term(42.0106), 1M(15.9949)"),
        ("PEPTIDEK", ""),
    ]
    print("\nFragPipe (Assigned Modifications):")
    for pep, mods in fp_tests:
        print(f"  {pep} + '{mods}' -> {standardize_fragpipe_sequence(pep, mods)}")

    # FragPipe Modified Peptide column
    fp_mod_tests = [
        "M[147.0354]PEPTIDEK",
        "n[42.0106]PEPTIDEK",
        "PEPTIDEK",
    ]
    print("\nFragPipe (Modified Peptide column):")
    for seq in fp_mod_tests:
        print(f"  {seq:30s} -> {standardize_fragpipe_modified_peptide(seq)}")

    # Cross-engine consistency
    print("\nCross-engine consistency (same peptide, different formats):")
    fp = standardize_fragpipe_sequence("CPEPTIDEK", "1C(57.0215)")
    dn = standardize_diann_sequence("C(UniMod:4)PEPTIDEK")
    sg = standardize_sage_sequence("C[+57.021465]PEPTIDEK")
    print(f"  FragPipe: {fp}")
    print(f"  DIA-NN:   {dn}")
    print(f"  Sage:     {sg}")
    print(f"  All match: {fp == dn == sg}")

    # N-term consistency
    print("\nN-term acetyl consistency:")
    fp_nt = standardize_fragpipe_sequence("MPEPTIDEK", "N-term(42.0106)")
    dn_nt = standardize_diann_sequence("(UniMod:1)MPEPTIDEK")
    sg_nt = standardize_sage_sequence("[+42.010567]-MPEPTIDEK")
    print(f"  FragPipe: {fp_nt}")
    print(f"  DIA-NN:   {dn_nt}")
    print(f"  Sage:     {sg_nt}")
    print(f"  All match: {fp_nt == dn_nt == sg_nt}")

    # remove_modifications with terminal format
    print("\nremove_modifications:")
    for seq in ["[]-PEPTIDEK-[]", "[UNIMOD:1]-M[UNIMOD:35]PEPTIDEK-[]", "PEPTIDE[UNIMOD:4]K"]:
        print(f"  {seq:45s} -> '{remove_modifications(seq)}'")

    # normalize_sequence_il with terminal format
    print("\nnormalize_sequence_il:")
    for seq in ["[]-PEPTIDEK-[]", "[UNIMOD:1]-MIPEPTIDEK-[]"]:
        print(f"  {seq:40s} -> '{normalize_sequence_il(seq)}'")
