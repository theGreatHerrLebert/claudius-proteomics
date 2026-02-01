#!/usr/bin/env python3
"""
Sequence Standardization Utilities

Converts modified peptide sequences from different search engines to a unified
UNIMOD format: PEPTIDE[UNIMOD:X]SEQUENCE

Supported formats:
- FragPipe: M[147.0354] (absolute mass) or 4C(57.0215) (position + delta)
- DIA-NN: (UniMod:35)
- Sage: [+15.9949] (mass delta)
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

    # TMT labels (approximate)
    229: (737, "TMT6plex"),           # TMT6plex + 229.1629
    304: (2016, "TMTpro"),            # TMTpro + 304.2071
}

# Absolute mass to UNIMOD (for FragPipe format)
ABSOLUTE_MASS_TO_UNIMOD: Dict[int, Tuple[int, str]] = {
    # Carbamidomethyl-Cysteine: 103.009 (C) + 57.021 = 160.030
    160: (4, "Carbamidomethyl"),

    # Oxidized Methionine: 131.040 (M) + 15.995 = 147.035
    147: (35, "Oxidation"),
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


def standardize_diann_sequence(modified_sequence: str) -> str:
    """Convert DIA-NN modified sequence to standard UNIMOD format.

    DIA-NN format: PEPTIDE(UniMod:35)K
    Standard format: PEPTIDE[UNIMOD:35]K

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

    return result


def standardize_sage_sequence(modified_sequence: str) -> str:
    """Convert Sage modified sequence to standard UNIMOD format.

    Sage format: PEPTIDE[+15.9949]K or PEPTIDE[+57.021465]K
    Standard format: PEPTIDE[UNIMOD:35]K or PEPTIDE[UNIMOD:4]K

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
            unimod = mass_to_unimod(abs(mass))
            if unimod:
                return unimod
            # Keep original if no mapping found
            return match.group(0)
        except ValueError:
            return match.group(0)

    return re.sub(pattern, replace_mass, modified_sequence)


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
        return peptide

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
        # If parsing fails, return plain sequence
        return peptide

    # Build modified sequence
    result = []

    # Add N-terminal modification
    if n_term_mod:
        result.append(n_term_mod)

    for i, aa in enumerate(peptide):
        result.append(aa)
        if i in mods_by_position:
            _, unimod = mods_by_position[i]
            result.append(unimod)

    return "".join(result)


def standardize_fragpipe_modified_peptide(modified_peptide: str) -> str:
    """Convert FragPipe 'Modified Peptide' column to standard UNIMOD format.

    FragPipe Modified Peptide format: M[147.0354]PEPTIDEK or n[42.0106]PEPTIDEK
    Standard format: M[UNIMOD:35]PEPTIDEK or [UNIMOD:1]PEPTIDEK

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

    # Find all mass modifications: X[XXX.XXXX]
    pattern = r'\[(\d+\.?\d*)\]'

    def replace_mass(match):
        mass_str = match.group(1)
        try:
            mass = float(mass_str)
            # This is absolute mass, need to figure out the delta
            # Common absolute masses:
            # 147 = oxidized M (131 + 16)
            # 160 = carbamidomethyl C (103 + 57)
            rounded = int(round(mass))

            if rounded in ABSOLUTE_MASS_TO_UNIMOD:
                unimod_id, _ = ABSOLUTE_MASS_TO_UNIMOD[rounded]
                return f"[UNIMOD:{unimod_id}]"

            # Try delta mass approach for some common ones
            # N-term acetyl is usually shown as [42.XXXX]
            if 41 <= rounded <= 43:
                return "[UNIMOD:1]"

            # Keep original if no mapping found
            return match.group(0)
        except ValueError:
            return match.group(0)

    return re.sub(pattern, replace_mass, result)


def remove_modifications(sequence: str) -> str:
    """Remove all modification annotations from sequence.

    Handles:
    - [UNIMOD:X]
    - [+X.XXX]
    - (UniMod:X)
    - [XXX.XXX] (absolute mass)

    Args:
        sequence: Modified sequence

    Returns:
        Plain sequence without modifications
    """
    if not sequence or not isinstance(sequence, str):
        return ""

    # Remove [UNIMOD:X]
    result = re.sub(r'\[UNIMOD:\d+\]', '', sequence)
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
    diann_seq = "PEPTIDE(UniMod:35)K"
    print(f"DIA-NN: {diann_seq}")
    print(f"  -> {standardize_diann_sequence(diann_seq)}")

    # Sage
    sage_seq = "PEPTIDE[+15.9949]K"
    print(f"\nSage: {sage_seq}")
    print(f"  -> {standardize_sage_sequence(sage_seq)}")

    sage_seq2 = "C[+57.021465]PEPTIDEK"
    print(f"Sage: {sage_seq2}")
    print(f"  -> {standardize_sage_sequence(sage_seq2)}")

    # FragPipe with modifications string
    fp_pep = "CPEPTIDEK"
    fp_mods = "1C(57.0215)"
    print(f"\nFragPipe: {fp_pep} + '{fp_mods}'")
    print(f"  -> {standardize_fragpipe_sequence(fp_pep, fp_mods)}")

    # FragPipe Modified Peptide column
    fp_mod = "M[147.0354]PEPTIDEK"
    print(f"\nFragPipe Modified: {fp_mod}")
    print(f"  -> {standardize_fragpipe_modified_peptide(fp_mod)}")

    # N-term acetyl
    fp_nterm = "n[42.0106]PEPTIDEK"
    print(f"FragPipe N-term: {fp_nterm}")
    print(f"  -> {standardize_fragpipe_modified_peptide(fp_nterm)}")

    # Sequence matching
    seq1 = "PEPTIDE[UNIMOD:35]K"
    seq2 = "PEPTIDE(UniMod:35)K"
    print(f"\nMatching: '{seq1}' vs '{seq2}'")
    print(f"  Plain: {remove_modifications(seq1)} vs {remove_modifications(seq2)}")
    print(f"  Match: {sequences_match(seq1, seq2)}")
