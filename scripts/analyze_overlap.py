#!/usr/bin/env python3
"""
Analyze overlap between search engine results (FragPipe, DIA-NN, Sage).

This script computes consensus statistics for a single PRIDE accession,
comparing peptide identifications from multiple search engines with proper
sequence normalization (I/L equivalence, UNIMOD formatting).

Supports 2-way (FragPipe + DIA-NN) and 3-way (+ Sage) comparisons.

Adapted from timsim-validate parsing and comparison utilities.

Usage:
    python scripts/analyze_overlap.py --accession PXD019086
    python scripts/analyze_overlap.py --accession PXD019086 --include-sage
    python scripts/analyze_overlap.py --accession PXD019086 --html --consensus
"""

import argparse
import base64
import json
import re
from pathlib import Path
from datetime import datetime
from typing import Tuple, Set, Dict, Any, Optional

import pandas as pd
import numpy as np


# =============================================================================
# Sequence Normalization (from timsim-validate)
# =============================================================================

def replace_I_with_L(sequence: str) -> str:
    """
    Replace all occurrences of 'I' with 'L' in the sequence, except inside brackets.

    This normalizes leucine/isoleucine which are indistinguishable by mass spectrometry.

    Args:
        sequence: The input peptide sequence.

    Returns:
        The modified string with 'I' replaced by 'L' outside of modification brackets.
    """
    if pd.isna(sequence):
        return ""

    def replace_match(match):
        if match.group(1):
            return match.group(1)  # Return text inside brackets unchanged
        return match.group(2).replace("I", "L")  # Replace outside brackets

    pattern = r'(\[.*?\])|([^\[\]]+)'
    modified_sequence = re.sub(pattern, replace_match, str(sequence))
    return modified_sequence.replace("UNLMOD", "UNIMOD")


def format_diann_sequence(modified_sequence: str) -> str:
    """
    Convert a DIA-NN modified sequence to standard UNIMOD format.

    DIA-NN uses parentheses and "UniMod", we convert to brackets and "UNIMOD".

    Args:
        modified_sequence: The modified peptide sequence from DIA-NN.

    Returns:
        The formatted peptide sequence with UNIMOD annotations.
    """
    if pd.isna(modified_sequence):
        return ""
    return str(modified_sequence).replace("(", "[").replace(")", "]").replace("UniMod", "UNIMOD")


def remove_unimod_annotation(sequence: str) -> str:
    """
    Remove UNIMOD annotations from a sequence.

    Args:
        sequence: Sequence with UNIMOD annotations like AAC[UNIMOD:4]PK

    Returns:
        Plain sequence without modifications: AACPK
    """
    if pd.isna(sequence):
        return ""
    # Remove [UNIMOD:X] and [+X.XXX] patterns
    clean = re.sub(r'\[UNIMOD:\d+\]', '', str(sequence))
    clean = re.sub(r'\[[+-]?\d+\.?\d*\]', '', clean)
    # Remove (UniMod:X) patterns (DIA-NN style)
    clean = re.sub(r'\(UniMod:\d+\)', '', clean)
    return clean


def mass_to_unimod(mass: float) -> str:
    """
    Convert a modification mass to UNIMOD annotation.

    Args:
        mass: The modification mass in Daltons.

    Returns:
        UNIMOD annotation string.
    """
    mass_map = {
        57: "[UNIMOD:4]",   # Carbamidomethyl (C)
        16: "[UNIMOD:35]",  # Oxidation (M)
        42: "[UNIMOD:1]",   # Acetyl (Protein N-term)
        80: "[UNIMOD:21]",  # Phospho (STY)
    }
    rounded = int(round(mass))
    return mass_map.get(rounded, f"[{mass:.4f}]")


def fragpipe_mods_to_unimod(sequence: str, mods: str) -> str:
    """
    Convert FragPipe modifications to UNIMOD-style sequence.

    Args:
        sequence: The plain peptide sequence.
        mods: The FragPipe modifications string (e.g., "2M(15.9949), 6C(57.0214)").

    Returns:
        Sequence with UNIMOD annotations.
    """
    if pd.isna(mods) or not mods or pd.isna(sequence):
        return str(sequence) if not pd.isna(sequence) else ""

    r_dict = {index: aa for index, aa in enumerate(str(sequence))}

    try:
        mods_list = str(mods).split(", ")
        for mod in mods_list:
            # Handle N-terminal acetylation
            if mod == 'N-term(42.0106)':
                r_dict[0] = f"[UNIMOD:1]{sequence[0]}"
            else:
                # Pattern: position + amino acid + mass in parentheses
                pattern = r"^(\d+)([A-Za-z])\(([\d.]+)\)$"
                match = re.match(pattern, mod)
                if match:
                    index, aa, mass = match.groups()
                    unimod = mass_to_unimod(float(mass))
                    r_dict[int(index) - 1] = aa + unimod
    except Exception:
        pass

    return "".join(r_dict.values())


def create_precursor_id(sequence: str, charge: int) -> str:
    """
    Create a unique precursor identifier from sequence and charge.

    Args:
        sequence: The peptide sequence (with or without modifications).
        charge: The precursor charge state.

    Returns:
        A unique precursor identifier string.
    """
    return f"{sequence}_{charge}"


def normalize_sequence_for_matching(sequence: str) -> str:
    """
    Normalize a sequence for matching between search engines.

    Applies I→L normalization and ensures consistent UNIMOD formatting.

    Args:
        sequence: The peptide sequence to normalize.

    Returns:
        Normalized sequence suitable for matching.
    """
    return replace_I_with_L(sequence)


# =============================================================================
# Data Loading
# =============================================================================

def load_fragpipe(
    accession: str,
    base_dir: Path = Path("data/processed"),
    level: str = "ion",
    pep_threshold: Optional[float] = None,
) -> pd.DataFrame:
    """
    Load FragPipe results with standardized columns.

    Args:
        accession: PRIDE accession
        base_dir: Base directory for processed data
        level: 'peptide' for combined_peptide.tsv, 'ion' for combined_ion.tsv (has charge)

    Returns:
        DataFrame with standardized columns
    """
    if level == "ion":
        path = base_dir / accession / "combined_ion.tsv"
    else:
        path = base_dir / accession / "combined_peptide.tsv"

    if not path.exists():
        raise FileNotFoundError(f"FragPipe results not found: {path}")

    df = pd.read_csv(path, sep="\t")

    # Standardize column names
    col_map = {
        "Peptide": "sequence",
        "Peptide Sequence": "sequence",
        "Modified Sequence": "modified_sequence_raw",
        "Charge": "charge",
        "M/Z": "mz",
        "Calibrated Observed M/Z": "mz",
        "Retention": "rt",
        "Ion Mobility": "mobility",
        "1/K0": "mobility",
        "Assigned Modifications": "modifications",
        "Protein": "protein",
        "Protein ID": "protein_id",
    }

    df = df.rename(columns={k: v for k, v in col_map.items() if k in df.columns})

    # Get plain sequence
    if "sequence" not in df.columns and "Peptide" in df.columns:
        df["sequence"] = df["Peptide"]

    # Convert modifications to UNIMOD format
    if "modifications" in df.columns:
        df["sequence_modified"] = df.apply(
            lambda row: fragpipe_mods_to_unimod(row.get("sequence", ""), row.get("modifications", "")),
            axis=1
        )
    elif "modified_sequence_raw" in df.columns:
        # Already in modified format, just standardize
        df["sequence_modified"] = df["modified_sequence_raw"].apply(format_diann_sequence)
    else:
        df["sequence_modified"] = df["sequence"]

    # Create normalized versions for matching
    df["sequence_clean"] = df["sequence"].apply(remove_unimod_annotation)
    df["sequence_normalized"] = df["sequence_clean"].apply(normalize_sequence_for_matching)
    df["sequence_modified_normalized"] = df["sequence_modified"].apply(normalize_sequence_for_matching)

    # Filter by PEP if threshold provided
    # FragPipe uses "Probability" = 1 - PEP, so filter Probability >= (1 - pep_threshold)
    if pep_threshold is not None:
        if "Probability" in df.columns:
            # Probability = 1 - PEP, so for PEP <= 0.05, need Probability >= 0.95
            prob_threshold = 1.0 - pep_threshold
            before = len(df)
            df = df[df["Probability"] >= prob_threshold]
            print(f"  PEP filter (Probability >= {prob_threshold}): {before:,} -> {len(df):,}")
        elif "PEP" in df.columns or "pep" in df.columns:
            pep_col = "PEP" if "PEP" in df.columns else "pep"
            before = len(df)
            df = df[df[pep_col] <= pep_threshold]
            print(f"  PEP filter ({pep_col} <= {pep_threshold}): {before:,} -> {len(df):,}")
        else:
            # Try to load Probability from PSM files
            psm_dir = base_dir / accession
            psm_files = list(psm_dir.rglob("psm.tsv"))
            if psm_files:
                print(f"  Loading Probability from {len(psm_files)} PSM files...")
                psm_dfs = []
                for psm_file in psm_files:
                    psm_df = pd.read_csv(psm_file, sep="\t", usecols=lambda c: c in ["Peptide", "Charge", "Probability"])
                    psm_dfs.append(psm_df)
                psm_all = pd.concat(psm_dfs, ignore_index=True)
                # Get max probability per peptide+charge (best PSM)
                psm_best = psm_all.groupby(["Peptide", "Charge"])["Probability"].max().reset_index()
                psm_best = psm_best.rename(columns={"Peptide": "sequence", "Charge": "charge"})

                # Merge with main df
                before = len(df)
                df = df.merge(psm_best, on=["sequence", "charge"], how="left")
                prob_threshold = 1.0 - pep_threshold
                df = df[df["Probability"] >= prob_threshold]
                print(f"  PEP filter (Probability >= {prob_threshold}): {before:,} -> {len(df):,}")

    df["source"] = "fragpipe"
    return df


def load_diann(
    accession: str,
    base_dir: Path = Path("data/processed"),
    fdr_threshold: float = 0.01,
    pep_threshold: Optional[float] = None,
) -> pd.DataFrame:
    """
    Load DIA-NN report results with standardized columns.

    Args:
        accession: PRIDE accession
        base_dir: Base directory for processed data
        fdr_threshold: FDR threshold for filtering

    Returns:
        DataFrame with standardized columns
    """
    path = base_dir / accession / "diann" / "report.parquet"
    if not path.exists():
        raise FileNotFoundError(f"DIA-NN results not found: {path}")

    df = pd.read_parquet(path)

    # Standardize column names
    col_map = {
        "Stripped.Sequence": "sequence",
        "Modified.Sequence": "modified_sequence_raw",
        "Precursor.Charge": "charge",
        "Precursor.Mz": "mz",
        "RT": "rt",
        "IM": "mobility",
        "CCS": "ccs",
        "Q.Value": "q_value",
        "PEP": "pep",
        "Protein.Ids": "protein",
        "Proteotypic": "proteotypic",
        "Precursor.Quantity": "intensity",
    }

    df = df.rename(columns={k: v for k, v in col_map.items() if k in df.columns})

    # Filter by FDR
    if "q_value" in df.columns:
        df = df[df["q_value"] <= fdr_threshold]

    # Filter by PEP if threshold provided
    if pep_threshold is not None and "pep" in df.columns:
        before = len(df)
        df = df[df["pep"] <= pep_threshold]
        print(f"  PEP filter (pep <= {pep_threshold}): {before:,} -> {len(df):,}")

    # Convert modified sequence to UNIMOD format
    if "modified_sequence_raw" in df.columns:
        df["sequence_modified"] = df["modified_sequence_raw"].apply(format_diann_sequence)
    else:
        df["sequence_modified"] = df["sequence"]

    # Create normalized versions for matching
    df["sequence_clean"] = df["sequence"].str.upper()
    df["sequence_normalized"] = df["sequence_clean"].apply(normalize_sequence_for_matching)
    df["sequence_modified_normalized"] = df["sequence_modified"].apply(normalize_sequence_for_matching)

    df["source"] = "diann"
    return df


def load_sage(
    accession: str,
    base_dir: Path = Path("data/processed"),
    fdr_threshold: float = 0.01,
    pep_threshold: Optional[float] = None,
) -> pd.DataFrame:
    """
    Load Sage results with standardized columns.

    Args:
        accession: PRIDE accession
        base_dir: Base directory for processed data
        fdr_threshold: FDR threshold for filtering (spectrum_q)

    Returns:
        DataFrame with standardized columns
    """
    path = base_dir / accession / "sage" / "results.sage.parquet"
    if not path.exists():
        raise FileNotFoundError(f"Sage results not found: {path}")

    df = pd.read_parquet(path)

    # Standardize column names
    col_map = {
        "stripped_peptide": "sequence",
        "peptide": "modified_sequence_raw",
        "charge": "charge",
        "expmass": "mz",  # Actually precursor mass, will convert
        "rt": "rt",
        "ion_mobility": "mobility",
        "predicted_mobility": "predicted_mobility",
        "spectrum_q": "q_value",
        "posterior_error": "pep",
        "proteins": "protein",
        "hyperscore": "hyperscore",
    }

    df = df.rename(columns={k: v for k, v in col_map.items() if k in df.columns})

    # Filter by FDR
    if "q_value" in df.columns:
        df = df[df["q_value"] <= fdr_threshold]

    # Filter out decoys
    if "is_decoy" in df.columns:
        df = df[~df["is_decoy"]]

    # Convert posterior_error from log space to linear PEP
    if "pep" in df.columns:
        # Sage stores log(PEP), convert to linear
        df["pep"] = np.exp(df["pep"])

    # Filter by PEP if threshold provided
    if pep_threshold is not None and "pep" in df.columns:
        before = len(df)
        df = df[df["pep"] <= pep_threshold]
        print(f"  PEP filter (pep <= {pep_threshold}): {before:,} -> {len(df):,}")

    # Convert modified sequence format (Sage uses different notation)
    # Sage format: AAC[+57.0215]M[+15.9949]PK
    if "modified_sequence_raw" in df.columns:
        def sage_to_unimod(seq):
            if pd.isna(seq):
                return ""
            # Convert [+57.0215] to [UNIMOD:4], [+15.9949] to [UNIMOD:35], etc.
            result = str(seq)
            # Common modifications
            result = re.sub(r'\[\+57\.02\d*\]', '[UNIMOD:4]', result)  # Carbamidomethyl
            result = re.sub(r'\[\+15\.99\d*\]', '[UNIMOD:35]', result)  # Oxidation
            result = re.sub(r'\[\+42\.01\d*\]', '[UNIMOD:1]', result)  # Acetyl
            result = re.sub(r'\[\-17\.02\d*\]', '[UNIMOD:28]', result)  # Gln->pyro-Glu
            result = re.sub(r'\[\-18\.01\d*\]', '[UNIMOD:27]', result)  # Glu->pyro-Glu
            return result
        df["sequence_modified"] = df["modified_sequence_raw"].apply(sage_to_unimod)
    else:
        df["sequence_modified"] = df["sequence"]

    # Create normalized versions for matching
    df["sequence_clean"] = df["sequence"].str.upper()
    df["sequence_normalized"] = df["sequence_clean"].apply(normalize_sequence_for_matching)
    df["sequence_modified_normalized"] = df["sequence_modified"].apply(normalize_sequence_for_matching)

    df["source"] = "sage"
    return df


# =============================================================================
# Overlap Analysis
# =============================================================================

def create_peptide_sets(
    df: pd.DataFrame,
    use_modifications: bool = False,
) -> Tuple[Set[str], Set[Tuple[str, int]]]:
    """
    Create sets of peptides and precursors (peptide + charge) from a DataFrame.

    Args:
        df: DataFrame with sequence columns and charge
        use_modifications: If True, use modified sequences; otherwise plain sequences

    Returns:
        Tuple of (peptide_set, precursor_set)
    """
    if use_modifications:
        seq_col = "sequence_modified_normalized"
    else:
        seq_col = "sequence_normalized"

    peptide_set = set(df[seq_col].dropna())
    precursor_set = set(zip(df[seq_col].dropna(), df.loc[df[seq_col].notna(), "charge"]))

    return peptide_set, precursor_set


def compute_overlap_stats(
    fp_df: pd.DataFrame,
    dn_df: pd.DataFrame,
    match_charge: bool = True,
    use_modifications: bool = False,
) -> dict:
    """
    Compute overlap statistics between FragPipe and DIA-NN.

    Uses I/L normalization for proper matching.

    Args:
        fp_df: FragPipe peptides
        dn_df: DIA-NN peptides
        match_charge: If True, require charge state match (sequence+charge)
        use_modifications: If True, match on modified sequences

    Returns:
        Dictionary of overlap statistics
    """
    fp_peptides, fp_precursors = create_peptide_sets(fp_df, use_modifications)
    dn_peptides, dn_precursors = create_peptide_sets(dn_df, use_modifications)

    if match_charge:
        fp_set = fp_precursors
        dn_set = dn_precursors
        match_type = "sequence+charge" + ("+mods" if use_modifications else "")
    else:
        fp_set = fp_peptides
        dn_set = dn_peptides
        match_type = "sequence_only" + ("+mods" if use_modifications else "")

    # Set operations
    overlap = fp_set & dn_set
    union = fp_set | dn_set
    fp_only = fp_set - dn_set
    dn_only = dn_set - fp_set

    # Jaccard similarity
    jaccard = len(overlap) / len(union) if union else 0.0

    # Overlap rates
    fp_overlap_rate = len(overlap) / len(fp_set) if fp_set else 0.0
    dn_overlap_rate = len(overlap) / len(dn_set) if dn_set else 0.0

    stats = {
        "match_type": match_type,
        "i_l_normalized": True,
        "n_fragpipe": len(fp_set),
        "n_diann": len(dn_set),
        "n_overlap": len(overlap),
        "n_union": len(union),
        "n_fragpipe_only": len(fp_only),
        "n_diann_only": len(dn_only),
        "jaccard_similarity": round(jaccard, 4),
        "fragpipe_overlap_rate": round(fp_overlap_rate, 4),
        "diann_overlap_rate": round(dn_overlap_rate, 4),
    }

    return stats


def compute_threeway_overlap_stats(
    fp_df: pd.DataFrame,
    dn_df: pd.DataFrame,
    sg_df: pd.DataFrame,
    match_charge: bool = True,
    use_modifications: bool = False,
) -> dict:
    """
    Compute 3-way overlap statistics between FragPipe, DIA-NN, and Sage.

    Uses I/L normalization for proper matching.

    Args:
        fp_df: FragPipe peptides
        dn_df: DIA-NN peptides
        sg_df: Sage peptides
        match_charge: If True, require charge state match
        use_modifications: If True, match on modified sequences

    Returns:
        Dictionary of overlap statistics
    """
    fp_peptides, fp_precursors = create_peptide_sets(fp_df, use_modifications)
    dn_peptides, dn_precursors = create_peptide_sets(dn_df, use_modifications)
    sg_peptides, sg_precursors = create_peptide_sets(sg_df, use_modifications)

    if match_charge:
        fp_set = fp_precursors
        dn_set = dn_precursors
        sg_set = sg_precursors
        match_type = "sequence+charge" + ("+mods" if use_modifications else "")
    else:
        fp_set = fp_peptides
        dn_set = dn_peptides
        sg_set = sg_peptides
        match_type = "sequence_only" + ("+mods" if use_modifications else "")

    # All possible regions in 3-way Venn diagram
    all_three = fp_set & dn_set & sg_set
    fp_dn_only = (fp_set & dn_set) - sg_set
    fp_sg_only = (fp_set & sg_set) - dn_set
    dn_sg_only = (dn_set & sg_set) - fp_set
    fp_only = fp_set - dn_set - sg_set
    dn_only = dn_set - fp_set - sg_set
    sg_only = sg_set - fp_set - dn_set

    # Pairwise overlaps
    fp_dn = fp_set & dn_set
    fp_sg = fp_set & sg_set
    dn_sg = dn_set & sg_set

    # Union
    union = fp_set | dn_set | sg_set

    stats = {
        "match_type": match_type,
        "i_l_normalized": True,
        "n_engines": 3,
        # Individual counts
        "n_fragpipe": len(fp_set),
        "n_diann": len(dn_set),
        "n_sage": len(sg_set),
        # 3-way overlap
        "n_all_three": len(all_three),
        # Pairwise (excluding 3-way)
        "n_fp_dn_only": len(fp_dn_only),
        "n_fp_sg_only": len(fp_sg_only),
        "n_dn_sg_only": len(dn_sg_only),
        # Unique to each
        "n_fragpipe_only": len(fp_only),
        "n_diann_only": len(dn_only),
        "n_sage_only": len(sg_only),
        # Totals
        "n_union": len(union),
        "n_at_least_two": len(fp_dn | fp_sg | dn_sg),
        # Pairwise overlaps (total)
        "n_fp_dn": len(fp_dn),
        "n_fp_sg": len(fp_sg),
        "n_dn_sg": len(dn_sg),
        # Rates
        "three_way_rate": round(len(all_three) / len(union), 4) if union else 0.0,
        "at_least_two_rate": round(len(fp_dn | fp_sg | dn_sg) / len(union), 4) if union else 0.0,
        # Individual overlap rates (% found by at least one other engine)
        "fragpipe_validation_rate": round((len(fp_set) - len(fp_only)) / len(fp_set), 4) if fp_set else 0.0,
        "diann_validation_rate": round((len(dn_set) - len(dn_only)) / len(dn_set), 4) if dn_set else 0.0,
        "sage_validation_rate": round((len(sg_set) - len(sg_only)) / len(sg_set), 4) if sg_set else 0.0,
    }

    return stats


def analyze_differences(
    fp_df: pd.DataFrame,
    dn_df: pd.DataFrame,
    sg_df: Optional[pd.DataFrame] = None,
) -> dict:
    """Analyze characteristics of engine-specific peptides."""

    analysis = {}
    dfs = {"fragpipe": fp_df, "diann": dn_df}
    if sg_df is not None:
        dfs["sage"] = sg_df

    for name, df in dfs.items():
        # Sequence length distribution
        if "sequence_clean" in df.columns:
            lengths = df["sequence_clean"].str.len()
            analysis[f"{name}_seq_length"] = {
                "mean": round(float(lengths.mean()), 1),
                "median": int(lengths.median()),
                "min": int(lengths.min()),
                "max": int(lengths.max()),
            }

        # Charge state distribution
        if "charge" in df.columns:
            charges = df["charge"].value_counts().to_dict()
            analysis[f"{name}_charges"] = {str(k): int(v) for k, v in sorted(charges.items())}

        # Mobility ranges
        if "mobility" in df.columns:
            mob = df["mobility"].dropna()
            if len(mob) > 0:
                analysis[f"{name}_mobility"] = {
                    "mean": round(float(mob.mean()), 4),
                    "min": round(float(mob.min()), 4),
                    "max": round(float(mob.max()), 4),
                }

    return analysis


def create_consensus_dataframes(
    fp_df: pd.DataFrame,
    dn_df: pd.DataFrame,
    match_charge: bool = True,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Create OVERLAP and UNION consensus DataFrames.

    Args:
        fp_df: FragPipe results
        dn_df: DIA-NN results
        match_charge: Match on sequence+charge

    Returns:
        Tuple of (overlap_df, union_df)
    """
    # Create match keys
    fp_df = fp_df.copy()
    dn_df = dn_df.copy()

    if match_charge:
        fp_df["match_key"] = fp_df["sequence_normalized"] + "_" + fp_df["charge"].astype(str)
        dn_df["match_key"] = dn_df["sequence_normalized"] + "_" + dn_df["charge"].astype(str)
    else:
        fp_df["match_key"] = fp_df["sequence_normalized"]
        dn_df["match_key"] = dn_df["sequence_normalized"]

    fp_keys = set(fp_df["match_key"].dropna())
    dn_keys = set(dn_df["match_key"].dropna())

    overlap_keys = fp_keys & dn_keys
    union_keys = fp_keys | dn_keys

    # OVERLAP: Found by both - merge data
    fp_overlap = fp_df[fp_df["match_key"].isin(overlap_keys)].copy()
    dn_overlap = dn_df[dn_df["match_key"].isin(overlap_keys)].copy()

    # Keep best from each (deduplicate)
    fp_overlap = fp_overlap.drop_duplicates(subset=["match_key"], keep="first")
    dn_overlap = dn_overlap.drop_duplicates(subset=["match_key"], keep="first")

    # Merge - only include columns that exist
    fp_cols = ["match_key", "sequence_normalized", "sequence_modified", "charge", "mz", "rt", "mobility", "protein"]
    fp_cols = [c for c in fp_cols if c in fp_overlap.columns]
    overlap_df = fp_overlap[fp_cols].copy()
    # Rename columns if they exist
    rename_map = {"rt": "rt_fragpipe", "mobility": "mobility_fragpipe", "mz": "mz_fragpipe"}
    overlap_df = overlap_df.rename(columns={k: v for k, v in rename_map.items() if k in overlap_df.columns})

    dn_cols = ["match_key", "rt", "mobility", "mz", "q_value", "pep", "ccs"]
    dn_cols = [c for c in dn_cols if c in dn_overlap.columns]
    dn_subset = dn_overlap[dn_cols].copy()
    dn_subset = dn_subset.rename(columns={
        "rt": "rt_diann",
        "mobility": "mobility_diann",
        "mz": "mz_diann",
    })

    overlap_df = overlap_df.merge(dn_subset, on="match_key", how="left")
    overlap_df["source"] = "both"
    overlap_df["confidence_weight"] = 1.0

    # UNION: Found by either
    # Start with overlap
    union_records = []

    for key in union_keys:
        in_fp = key in fp_keys
        in_dn = key in dn_keys

        if in_fp and in_dn:
            # From overlap
            row = fp_overlap[fp_overlap["match_key"] == key].iloc[0]
            dn_row = dn_overlap[dn_overlap["match_key"] == key].iloc[0] if key in dn_overlap["match_key"].values else None
            record = {
                "match_key": key,
                "sequence_normalized": row["sequence_normalized"],
                "sequence_modified": row["sequence_modified"],
                "charge": row["charge"],
                "mz": row.get("mz"),
                "rt": row.get("rt"),
                "mobility": row.get("mobility"),
                "source": "both",
                "confidence_weight": 1.0,
            }
            if dn_row is not None and "ccs" in dn_row:
                record["ccs"] = dn_row.get("ccs")
                record["q_value"] = dn_row.get("q_value")
        elif in_fp:
            row = fp_df[fp_df["match_key"] == key].iloc[0]
            record = {
                "match_key": key,
                "sequence_normalized": row["sequence_normalized"],
                "sequence_modified": row["sequence_modified"],
                "charge": row["charge"],
                "mz": row.get("mz"),
                "rt": row.get("rt"),
                "mobility": row.get("mobility"),
                "source": "fragpipe",
                "confidence_weight": 0.7,
            }
        else:
            row = dn_df[dn_df["match_key"] == key].iloc[0]
            record = {
                "match_key": key,
                "sequence_normalized": row["sequence_normalized"],
                "sequence_modified": row["sequence_modified"],
                "charge": row["charge"],
                "mz": row.get("mz"),
                "rt": row.get("rt"),
                "mobility": row.get("mobility"),
                "ccs": row.get("ccs"),
                "q_value": row.get("q_value"),
                "source": "diann",
                "confidence_weight": 0.7,
            }

        union_records.append(record)

    union_df = pd.DataFrame(union_records)

    return overlap_df, union_df


# =============================================================================
# HTML Report Generation
# =============================================================================

HTML_TEMPLATE = """<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Overlap Analysis: {accession}</title>
    <style>
        :root {{
            --pass-color: #27ae60;
            --diann-color: #3498db;
            --fragpipe-color: #e74c3c;
            --overlap-color: #9b59b6;
            --bg-light: #f8f9fa;
            --border-color: #dee2e6;
        }}
        * {{ box-sizing: border-box; }}
        body {{
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
            line-height: 1.6;
            max-width: 1200px;
            margin: 0 auto;
            padding: 20px;
            background: #fff;
            color: #333;
        }}
        h1 {{
            color: #2c3e50;
            border-bottom: 3px solid var(--overlap-color);
            padding-bottom: 10px;
        }}
        h2 {{
            color: #34495e;
            border-bottom: 2px solid var(--border-color);
            padding-bottom: 8px;
            margin-top: 30px;
        }}
        .header-info {{
            background: var(--bg-light);
            padding: 15px;
            border-radius: 8px;
            margin-bottom: 20px;
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 10px;
        }}
        .header-info .label {{
            font-weight: 600;
            color: #666;
            font-size: 0.85em;
        }}
        .summary-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(250px, 1fr));
            gap: 20px;
            margin: 20px 0;
        }}
        .stat-card {{
            background: var(--bg-light);
            border-radius: 10px;
            padding: 20px;
            text-align: center;
        }}
        .stat-card.fragpipe {{ border-left: 4px solid var(--fragpipe-color); }}
        .stat-card.diann {{ border-left: 4px solid var(--diann-color); }}
        .stat-card.overlap {{ border-left: 4px solid var(--overlap-color); }}
        .stat-card .number {{
            font-size: 2.5em;
            font-weight: bold;
            color: #333;
        }}
        .stat-card .label {{
            font-size: 0.9em;
            color: #666;
            margin-top: 5px;
        }}
        .metric-row {{
            display: flex;
            justify-content: space-between;
            padding: 8px 0;
            border-bottom: 1px solid var(--border-color);
        }}
        .metric-row:last-child {{ border-bottom: none; }}
        .metric-label {{ font-weight: 500; color: #555; }}
        .metric-value {{ font-weight: 600; font-family: monospace; }}
        table {{
            width: 100%;
            border-collapse: collapse;
            margin: 15px 0;
        }}
        th, td {{
            padding: 12px 15px;
            text-align: left;
            border-bottom: 1px solid var(--border-color);
        }}
        th {{ background: var(--bg-light); font-weight: 600; }}
        .bar-container {{
            display: flex;
            height: 30px;
            border-radius: 5px;
            overflow: hidden;
            margin: 20px 0;
        }}
        .bar-fp {{ background: var(--fragpipe-color); }}
        .bar-overlap {{ background: var(--overlap-color); }}
        .bar-dn {{ background: var(--diann-color); }}
        .bar-segment {{
            display: flex;
            align-items: center;
            justify-content: center;
            color: white;
            font-weight: bold;
            font-size: 0.85em;
        }}
        .footer {{
            margin-top: 40px;
            padding-top: 20px;
            border-top: 2px solid var(--border-color);
            text-align: center;
            color: #666;
            font-size: 0.9em;
        }}
    </style>
</head>
<body>
    <h1>Overlap Analysis: {accession}</h1>

    <div class="header-info">
        <div><div class="label">Accession</div><div>{accession}</div></div>
        <div><div class="label">Generated</div><div>{timestamp}</div></div>
        <div><div class="label">Match Type</div><div>{match_type}</div></div>
        <div><div class="label">I/L Normalized</div><div>Yes</div></div>
    </div>

    <h2>Summary</h2>
    <div class="summary-grid">
        <div class="stat-card fragpipe">
            <div class="number">{n_fragpipe:,}</div>
            <div class="label">FragPipe Precursors</div>
        </div>
        <div class="stat-card diann">
            <div class="number">{n_diann:,}</div>
            <div class="label">DIA-NN Precursors</div>
        </div>
        <div class="stat-card overlap">
            <div class="number">{n_overlap:,}</div>
            <div class="label">OVERLAP (both)</div>
        </div>
        <div class="stat-card">
            <div class="number">{n_union:,}</div>
            <div class="label">UNION (either)</div>
        </div>
    </div>

    <h2>Overlap Visualization</h2>
    <div class="bar-container">
        <div class="bar-segment bar-fp" style="width: {fp_only_pct}%">{n_fragpipe_only:,}</div>
        <div class="bar-segment bar-overlap" style="width: {overlap_pct}%">{n_overlap:,}</div>
        <div class="bar-segment bar-dn" style="width: {dn_only_pct}%">{n_diann_only:,}</div>
    </div>
    <div style="display: flex; justify-content: space-between; font-size: 0.85em; color: #666;">
        <span>FragPipe only ({fp_only_pct:.1f}%)</span>
        <span>Overlap ({overlap_pct:.1f}%)</span>
        <span>DIA-NN only ({dn_only_pct:.1f}%)</span>
    </div>

    <h2>Metrics</h2>
    <table>
        <tr><th>Metric</th><th>Value</th><th>Description</th></tr>
        <tr><td>Jaccard Similarity</td><td>{jaccard:.1%}</td><td>Overlap / Union</td></tr>
        <tr><td>FragPipe Overlap Rate</td><td>{fp_overlap_rate:.1%}</td><td>% of FragPipe found by DIA-NN</td></tr>
        <tr><td>DIA-NN Overlap Rate</td><td>{dn_overlap_rate:.1%}</td><td>% of DIA-NN found by FragPipe</td></tr>
    </table>

    <h2>Charge Distribution</h2>
    <table>
        <tr><th>Charge</th><th>FragPipe</th><th>DIA-NN</th></tr>
        {charge_rows}
    </table>

    <h2>Sequence Length</h2>
    <table>
        <tr><th>Statistic</th><th>FragPipe</th><th>DIA-NN</th></tr>
        <tr><td>Min</td><td>{fp_len_min}</td><td>{dn_len_min}</td></tr>
        <tr><td>Median</td><td>{fp_len_median}</td><td>{dn_len_median}</td></tr>
        <tr><td>Max</td><td>{fp_len_max}</td><td>{dn_len_max}</td></tr>
    </table>

    <div class="footer">
        <p>Generated by <strong>San Jos&eacute; Pipeline</strong> | analyze_overlap.py</p>
        <p>Using I/L normalization and UNIMOD standardization from timsim-validate</p>
    </div>
</body>
</html>
"""

HTML_TEMPLATE_3WAY = """<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>3-Way Overlap Analysis: {accession}</title>
    <style>
        :root {{
            --fragpipe-color: #e74c3c;
            --diann-color: #3498db;
            --sage-color: #27ae60;
            --all-three-color: #9b59b6;
            --two-engines-color: #f39c12;
            --bg-light: #f8f9fa;
            --border-color: #dee2e6;
        }}
        * {{ box-sizing: border-box; }}
        body {{
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
            line-height: 1.6;
            max-width: 1200px;
            margin: 0 auto;
            padding: 20px;
            background: #fff;
            color: #333;
        }}
        h1 {{
            color: #2c3e50;
            border-bottom: 3px solid var(--all-three-color);
            padding-bottom: 10px;
        }}
        h2 {{
            color: #34495e;
            border-bottom: 2px solid var(--border-color);
            padding-bottom: 8px;
            margin-top: 30px;
        }}
        .header-info {{
            background: var(--bg-light);
            padding: 15px;
            border-radius: 8px;
            margin-bottom: 20px;
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(180px, 1fr));
            gap: 10px;
        }}
        .header-info .label {{
            font-weight: 600;
            color: #666;
            font-size: 0.85em;
        }}
        .summary-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 15px;
            margin: 20px 0;
        }}
        .stat-card {{
            background: var(--bg-light);
            border-radius: 10px;
            padding: 15px;
            text-align: center;
        }}
        .stat-card.fragpipe {{ border-left: 4px solid var(--fragpipe-color); }}
        .stat-card.diann {{ border-left: 4px solid var(--diann-color); }}
        .stat-card.sage {{ border-left: 4px solid var(--sage-color); }}
        .stat-card.all-three {{ border-left: 4px solid var(--all-three-color); background: #f3e5f5; }}
        .stat-card.two-plus {{ border-left: 4px solid var(--two-engines-color); }}
        .stat-card .number {{
            font-size: 2em;
            font-weight: bold;
            color: #333;
        }}
        .stat-card .label {{
            font-size: 0.85em;
            color: #666;
            margin-top: 5px;
        }}
        .stat-card .pct {{
            font-size: 0.9em;
            color: #888;
        }}
        table {{
            width: 100%;
            border-collapse: collapse;
            margin: 15px 0;
        }}
        th, td {{
            padding: 12px 15px;
            text-align: left;
            border-bottom: 1px solid var(--border-color);
        }}
        th {{ background: var(--bg-light); font-weight: 600; }}
        .highlight {{ background: #e8f5e9; font-weight: 600; }}
        .venn-container {{
            display: flex;
            justify-content: center;
            margin: 30px 0;
        }}
        .venn-svg {{
            max-width: 500px;
        }}
        .bar-stacked {{
            display: flex;
            height: 40px;
            border-radius: 5px;
            overflow: hidden;
            margin: 20px 0;
        }}
        .bar-segment {{
            display: flex;
            align-items: center;
            justify-content: center;
            color: white;
            font-weight: bold;
            font-size: 0.8em;
            min-width: 30px;
        }}
        .bar-fp {{ background: var(--fragpipe-color); }}
        .bar-dn {{ background: var(--diann-color); }}
        .bar-sg {{ background: var(--sage-color); }}
        .bar-all {{ background: var(--all-three-color); }}
        .bar-fp-dn {{ background: #8e44ad; }}
        .bar-fp-sg {{ background: #c0392b; }}
        .bar-dn-sg {{ background: #16a085; }}
        .legend {{
            display: flex;
            flex-wrap: wrap;
            gap: 15px;
            justify-content: center;
            margin: 10px 0;
            font-size: 0.85em;
        }}
        .legend-item {{
            display: flex;
            align-items: center;
            gap: 5px;
        }}
        .legend-color {{
            width: 16px;
            height: 16px;
            border-radius: 3px;
        }}
        .validation-table td:nth-child(3) {{
            font-weight: bold;
        }}
        .footer {{
            margin-top: 40px;
            padding-top: 20px;
            border-top: 2px solid var(--border-color);
            text-align: center;
            color: #666;
            font-size: 0.9em;
        }}
    </style>
</head>
<body>
    <h1>3-Way Orthogonal Validation: {accession}</h1>

    <div class="header-info">
        <div><div class="label">Accession</div><div>{accession}</div></div>
        <div><div class="label">Generated</div><div>{timestamp}</div></div>
        <div><div class="label">Match Type</div><div>{match_type}</div></div>
        <div><div class="label">PEP Filter</div><div>{pep_filter}</div></div>
        <div><div class="label">I/L Normalized</div><div>Yes</div></div>
    </div>

    <h2>Search Engine Results</h2>
    <div class="summary-grid">
        <div class="stat-card fragpipe">
            <div class="number">{n_fragpipe:,}</div>
            <div class="label">FragPipe</div>
            <div class="pct">{fp_validation_rate:.1%} validated</div>
        </div>
        <div class="stat-card diann">
            <div class="number">{n_diann:,}</div>
            <div class="label">DIA-NN</div>
            <div class="pct">{dn_validation_rate:.1%} validated</div>
        </div>
        <div class="stat-card sage">
            <div class="number">{n_sage:,}</div>
            <div class="label">Sage</div>
            <div class="pct">{sg_validation_rate:.1%} validated</div>
        </div>
    </div>

    <h2>Consensus Summary</h2>
    <div class="summary-grid">
        <div class="stat-card all-three">
            <div class="number">{n_all_three:,}</div>
            <div class="label">All 3 Engines</div>
            <div class="pct">{three_way_pct:.1%} of union</div>
        </div>
        <div class="stat-card two-plus">
            <div class="number">{n_at_least_two:,}</div>
            <div class="label">At Least 2 Engines</div>
            <div class="pct">{two_plus_pct:.1%} of union</div>
        </div>
        <div class="stat-card">
            <div class="number">{n_union:,}</div>
            <div class="label">Union (Any Engine)</div>
            <div class="pct">100%</div>
        </div>
    </div>

    <h2>Overlap Breakdown</h2>
    <div class="bar-stacked">
        <div class="bar-segment bar-all" style="width: {all_three_pct}%" title="All 3">{n_all_three:,}</div>
        <div class="bar-segment bar-fp-dn" style="width: {fp_dn_only_pct}%" title="FP+DN only">{n_fp_dn_only:,}</div>
        <div class="bar-segment bar-fp-sg" style="width: {fp_sg_only_pct}%" title="FP+Sage only">{n_fp_sg_only:,}</div>
        <div class="bar-segment bar-dn-sg" style="width: {dn_sg_only_pct}%" title="DN+Sage only">{n_dn_sg_only:,}</div>
        <div class="bar-segment bar-fp" style="width: {fp_only_pct}%" title="FragPipe only">{n_fragpipe_only:,}</div>
        <div class="bar-segment bar-dn" style="width: {dn_only_pct}%" title="DIA-NN only">{n_diann_only:,}</div>
        <div class="bar-segment bar-sg" style="width: {sg_only_pct}%" title="Sage only">{n_sage_only:,}</div>
    </div>
    <div class="legend">
        <div class="legend-item"><div class="legend-color" style="background: var(--all-three-color)"></div>All 3</div>
        <div class="legend-item"><div class="legend-color" style="background: #8e44ad"></div>FP+DN</div>
        <div class="legend-item"><div class="legend-color" style="background: #c0392b"></div>FP+Sage</div>
        <div class="legend-item"><div class="legend-color" style="background: #16a085"></div>DN+Sage</div>
        <div class="legend-item"><div class="legend-color" style="background: var(--fragpipe-color)"></div>FP only</div>
        <div class="legend-item"><div class="legend-color" style="background: var(--diann-color)"></div>DN only</div>
        <div class="legend-item"><div class="legend-color" style="background: var(--sage-color)"></div>Sage only</div>
    </div>

    <h2>Pairwise Overlaps</h2>
    <table>
        <tr><th>Pair</th><th>Total Overlap</th><th>Exclusive (not in 3rd)</th></tr>
        <tr><td>FragPipe ∩ DIA-NN</td><td>{n_fp_dn:,}</td><td>{n_fp_dn_only:,}</td></tr>
        <tr><td>FragPipe ∩ Sage</td><td>{n_fp_sg:,}</td><td>{n_fp_sg_only:,}</td></tr>
        <tr><td>DIA-NN ∩ Sage</td><td>{n_dn_sg:,}</td><td>{n_dn_sg_only:,}</td></tr>
    </table>

    <h2>Validation Quality</h2>
    <table class="validation-table">
        <tr><th>Engine</th><th>Unique (not in others)</th><th>Validation Rate</th></tr>
        <tr><td>FragPipe</td><td>{n_fragpipe_only:,} ({fp_unique_pct:.1%})</td><td>{fp_validation_rate:.1%}</td></tr>
        <tr><td>DIA-NN</td><td>{n_diann_only:,} ({dn_unique_pct:.1%})</td><td>{dn_validation_rate:.1%}</td></tr>
        <tr class="highlight"><td>Sage</td><td>{n_sage_only:,} ({sg_unique_pct:.1%})</td><td>{sg_validation_rate:.1%}</td></tr>
    </table>

    <h2>Charge Distribution</h2>
    <table>
        <tr><th>Charge</th><th>FragPipe</th><th>DIA-NN</th><th>Sage</th></tr>
        {charge_rows}
    </table>

    <h2>Sequence Length</h2>
    <table>
        <tr><th>Statistic</th><th>FragPipe</th><th>DIA-NN</th><th>Sage</th></tr>
        <tr><td>Min</td><td>{fp_len_min}</td><td>{dn_len_min}</td><td>{sg_len_min}</td></tr>
        <tr><td>Median</td><td>{fp_len_median}</td><td>{dn_len_median}</td><td>{sg_len_median}</td></tr>
        <tr><td>Max</td><td>{fp_len_max}</td><td>{dn_len_max}</td><td>{sg_len_max}</td></tr>
    </table>

    <div class="footer">
        <p>Generated by <strong>San Jos&eacute; Pipeline</strong> | Triple Orthogonal Validation</p>
        <p>Using I/L normalization and UNIMOD standardization</p>
    </div>
</body>
</html>
"""


def generate_html_report(
    accession: str,
    stats: dict,
    analysis: dict,
    output_path: Path,
    pep_threshold: Optional[float] = None,
) -> str:
    """Generate a self-contained HTML report."""

    # Check if this is a 3-way comparison
    is_3way = "n_sage" in stats

    if is_3way:
        return generate_html_report_3way(accession, stats, analysis, output_path, pep_threshold)

    # 2-way report (original logic)
    total = stats["n_union"]
    fp_only_pct = (stats["n_fragpipe_only"] / total * 100) if total else 0
    overlap_pct = (stats["n_overlap"] / total * 100) if total else 0
    dn_only_pct = (stats["n_diann_only"] / total * 100) if total else 0

    # Build charge distribution rows
    fp_charges = analysis.get("fragpipe_charges", {})
    dn_charges = analysis.get("diann_charges", {})
    all_charges = sorted(set(fp_charges.keys()) | set(dn_charges.keys()), key=lambda x: int(x))
    charge_rows = "\n".join([
        f"<tr><td>{c}+</td><td>{fp_charges.get(c, 0):,}</td><td>{dn_charges.get(c, 0):,}</td></tr>"
        for c in all_charges
    ])

    # Sequence lengths
    fp_len = analysis.get("fragpipe_seq_length", {})
    dn_len = analysis.get("diann_seq_length", {})

    html = HTML_TEMPLATE.format(
        accession=accession,
        timestamp=datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        match_type=stats["match_type"],
        n_fragpipe=stats["n_fragpipe"],
        n_diann=stats["n_diann"],
        n_overlap=stats["n_overlap"],
        n_union=stats["n_union"],
        n_fragpipe_only=stats["n_fragpipe_only"],
        n_diann_only=stats["n_diann_only"],
        fp_only_pct=fp_only_pct,
        overlap_pct=overlap_pct,
        dn_only_pct=dn_only_pct,
        jaccard=stats["jaccard_similarity"],
        fp_overlap_rate=stats["fragpipe_overlap_rate"],
        dn_overlap_rate=stats["diann_overlap_rate"],
        charge_rows=charge_rows,
        fp_len_min=fp_len.get("min", "N/A"),
        fp_len_median=fp_len.get("median", "N/A"),
        fp_len_max=fp_len.get("max", "N/A"),
        dn_len_min=dn_len.get("min", "N/A"),
        dn_len_median=dn_len.get("median", "N/A"),
        dn_len_max=dn_len.get("max", "N/A"),
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w") as f:
        f.write(html)

    return str(output_path)


def generate_html_report_3way(
    accession: str,
    stats: dict,
    analysis: dict,
    output_path: Path,
    pep_threshold: Optional[float] = None,
) -> str:
    """Generate a 3-way comparison HTML report."""

    total = stats["n_union"]

    # Calculate percentages
    all_three_pct = (stats["n_all_three"] / total * 100) if total else 0
    fp_dn_only_pct = (stats["n_fp_dn_only"] / total * 100) if total else 0
    fp_sg_only_pct = (stats["n_fp_sg_only"] / total * 100) if total else 0
    dn_sg_only_pct = (stats["n_dn_sg_only"] / total * 100) if total else 0
    fp_only_pct = (stats["n_fragpipe_only"] / total * 100) if total else 0
    dn_only_pct = (stats["n_diann_only"] / total * 100) if total else 0
    sg_only_pct = (stats["n_sage_only"] / total * 100) if total else 0

    # Unique percentages (of engine's own count)
    fp_unique_pct = (stats["n_fragpipe_only"] / stats["n_fragpipe"] * 100) if stats["n_fragpipe"] else 0
    dn_unique_pct = (stats["n_diann_only"] / stats["n_diann"] * 100) if stats["n_diann"] else 0
    sg_unique_pct = (stats["n_sage_only"] / stats["n_sage"] * 100) if stats["n_sage"] else 0

    # Build charge distribution rows (3 columns)
    fp_charges = analysis.get("fragpipe_charges", {})
    dn_charges = analysis.get("diann_charges", {})
    sg_charges = analysis.get("sage_charges", {})
    all_charges = sorted(
        set(fp_charges.keys()) | set(dn_charges.keys()) | set(sg_charges.keys()),
        key=lambda x: int(x)
    )
    charge_rows = "\n".join([
        f"<tr><td>{c}+</td><td>{fp_charges.get(c, 0):,}</td><td>{dn_charges.get(c, 0):,}</td><td>{sg_charges.get(c, 0):,}</td></tr>"
        for c in all_charges
    ])

    # Sequence lengths
    fp_len = analysis.get("fragpipe_seq_length", {})
    dn_len = analysis.get("diann_seq_length", {})
    sg_len = analysis.get("sage_seq_length", {})

    # PEP filter description
    pep_filter = f"≤ {pep_threshold} ({(1-pep_threshold)*100:.0f}% conf.)" if pep_threshold else "None (q-value only)"

    html = HTML_TEMPLATE_3WAY.format(
        accession=accession,
        timestamp=datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        match_type=stats["match_type"],
        pep_filter=pep_filter,
        n_fragpipe=stats["n_fragpipe"],
        n_diann=stats["n_diann"],
        n_sage=stats["n_sage"],
        n_all_three=stats["n_all_three"],
        n_at_least_two=stats["n_at_least_two"],
        n_union=stats["n_union"],
        three_way_pct=stats["three_way_rate"],
        two_plus_pct=stats["at_least_two_rate"],
        # Pairwise
        n_fp_dn=stats["n_fp_dn"],
        n_fp_sg=stats["n_fp_sg"],
        n_dn_sg=stats["n_dn_sg"],
        n_fp_dn_only=stats["n_fp_dn_only"],
        n_fp_sg_only=stats["n_fp_sg_only"],
        n_dn_sg_only=stats["n_dn_sg_only"],
        # Unique
        n_fragpipe_only=stats["n_fragpipe_only"],
        n_diann_only=stats["n_diann_only"],
        n_sage_only=stats["n_sage_only"],
        # Percentages for bar
        all_three_pct=all_three_pct,
        fp_dn_only_pct=fp_dn_only_pct,
        fp_sg_only_pct=fp_sg_only_pct,
        dn_sg_only_pct=dn_sg_only_pct,
        fp_only_pct=fp_only_pct,
        dn_only_pct=dn_only_pct,
        sg_only_pct=sg_only_pct,
        # Unique percentages
        fp_unique_pct=fp_unique_pct,
        dn_unique_pct=dn_unique_pct,
        sg_unique_pct=sg_unique_pct,
        # Validation rates
        fp_validation_rate=stats["fragpipe_validation_rate"],
        dn_validation_rate=stats["diann_validation_rate"],
        sg_validation_rate=stats["sage_validation_rate"],
        # Charge rows
        charge_rows=charge_rows,
        # Sequence lengths
        fp_len_min=fp_len.get("min", "N/A"),
        fp_len_median=fp_len.get("median", "N/A"),
        fp_len_max=fp_len.get("max", "N/A"),
        dn_len_min=dn_len.get("min", "N/A"),
        dn_len_median=dn_len.get("median", "N/A"),
        dn_len_max=dn_len.get("max", "N/A"),
        sg_len_min=sg_len.get("min", "N/A"),
        sg_len_median=sg_len.get("median", "N/A"),
        sg_len_max=sg_len.get("max", "N/A"),
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w") as f:
        f.write(html)

    return str(output_path)


# =============================================================================
# CLI
# =============================================================================

def print_summary(accession: str, stats: dict, analysis: dict):
    """Print a formatted summary to console."""
    print(f"\n{'='*60}")
    print(f"  Overlap Analysis: {accession}")
    print(f"{'='*60}\n")

    print(f"  Match type: {stats['match_type']}")
    print(f"  I/L normalized: {stats.get('i_l_normalized', True)}\n")

    print("  Precursor Counts:")
    print(f"    FragPipe:      {stats['n_fragpipe']:>8,}")
    print(f"    DIA-NN:        {stats['n_diann']:>8,}")
    print(f"    Overlap:       {stats['n_overlap']:>8,}  (found by both)")
    print(f"    Union:         {stats['n_union']:>8,}  (found by either)")
    print(f"    FragPipe only: {stats['n_fragpipe_only']:>8,}")
    print(f"    DIA-NN only:   {stats['n_diann_only']:>8,}")

    print(f"\n  Similarity Metrics:")
    print(f"    Jaccard:           {stats['jaccard_similarity']:.1%}")
    print(f"    FragPipe overlap:  {stats['fragpipe_overlap_rate']:.1%}  (% of FP found by DN)")
    print(f"    DIA-NN overlap:    {stats['diann_overlap_rate']:.1%}  (% of DN found by FP)")

    if "fragpipe_seq_length" in analysis:
        print(f"\n  Sequence Length:")
        fp = analysis["fragpipe_seq_length"]
        dn = analysis.get("diann_seq_length", {})
        print(f"    FragPipe: {fp['min']}-{fp['max']} aa (median: {fp['median']})")
        if dn:
            print(f"    DIA-NN:   {dn['min']}-{dn['max']} aa (median: {dn['median']})")

    if "fragpipe_charges" in analysis:
        print(f"\n  Charge Distribution:")
        print(f"    FragPipe: {analysis['fragpipe_charges']}")
        if "diann_charges" in analysis:
            print(f"    DIA-NN:   {analysis['diann_charges']}")

    print(f"\n{'='*60}\n")


def print_threeway_summary(accession: str, stats: dict, analysis: dict):
    """Print a formatted 3-way summary to console."""
    print(f"\n{'='*70}")
    print(f"  3-Way Overlap Analysis: {accession}")
    print(f"{'='*70}\n")

    print(f"  Match type: {stats['match_type']}")
    print(f"  I/L normalized: {stats.get('i_l_normalized', True)}\n")

    print("  Precursor Counts:")
    print(f"    FragPipe:        {stats['n_fragpipe']:>8,}")
    print(f"    DIA-NN:          {stats['n_diann']:>8,}")
    print(f"    Sage:            {stats['n_sage']:>8,}")
    print(f"    ---")
    print(f"    All three:       {stats['n_all_three']:>8,}  (found by all engines)")
    print(f"    At least two:    {stats['n_at_least_two']:>8,}  (found by 2+ engines)")
    print(f"    Union:           {stats['n_union']:>8,}  (found by any)")

    print(f"\n  Pairwise Overlaps:")
    print(f"    FP ∩ DN:         {stats['n_fp_dn']:>8,}")
    print(f"    FP ∩ Sage:       {stats['n_fp_sg']:>8,}")
    print(f"    DN ∩ Sage:       {stats['n_dn_sg']:>8,}")

    print(f"\n  Unique to Each Engine:")
    print(f"    FragPipe only:   {stats['n_fragpipe_only']:>8,}")
    print(f"    DIA-NN only:     {stats['n_diann_only']:>8,}")
    print(f"    Sage only:       {stats['n_sage_only']:>8,}")

    print(f"\n  Validation Rates (% found by at least one other engine):")
    print(f"    FragPipe:        {stats['fragpipe_validation_rate']:.1%}")
    print(f"    DIA-NN:          {stats['diann_validation_rate']:.1%}")
    print(f"    Sage:            {stats['sage_validation_rate']:.1%}")

    print(f"\n  Consensus Quality:")
    print(f"    3-way agreement: {stats['three_way_rate']:.1%}  (% of union)")
    print(f"    2+ agreement:    {stats['at_least_two_rate']:.1%}  (% of union)")

    print(f"\n{'='*70}\n")


def main():
    parser = argparse.ArgumentParser(
        description="Analyze overlap between search engine results"
    )
    parser.add_argument(
        "--accession", "-a",
        required=True,
        help="PRIDE accession to analyze (e.g., PXD019086)"
    )
    parser.add_argument(
        "--output", "-o",
        type=Path,
        help="Output JSON file for statistics"
    )
    parser.add_argument(
        "--base-dir", "-d",
        type=Path,
        default=Path("data/processed"),
        help="Base directory for processed data"
    )
    parser.add_argument(
        "--include-sage",
        action="store_true",
        help="Include Sage in 3-way comparison"
    )
    parser.add_argument(
        "--pep-threshold",
        type=float,
        default=None,
        help="PEP threshold filter (e.g., 0.05 for 95%% confidence)"
    )
    parser.add_argument(
        "--no-charge",
        action="store_true",
        help="Match on sequence only, ignore charge state"
    )
    parser.add_argument(
        "--with-mods",
        action="store_true",
        help="Include modifications in matching"
    )
    parser.add_argument(
        "--html",
        action="store_true",
        help="Generate HTML report"
    )
    parser.add_argument(
        "--consensus",
        action="store_true",
        help="Generate consensus parquet files (overlap.parquet, union.parquet)"
    )
    parser.add_argument(
        "--quiet", "-q",
        action="store_true",
        help="Suppress console output"
    )

    args = parser.parse_args()

    # Load data
    pep_thresh = args.pep_threshold
    if pep_thresh is not None and not args.quiet:
        print(f"Applying PEP threshold: {pep_thresh} ({(1-pep_thresh)*100:.0f}% confidence)")

    if not args.quiet:
        print(f"Loading FragPipe results for {args.accession}...")
    fp_df = load_fragpipe(args.accession, args.base_dir, pep_threshold=pep_thresh)

    if not args.quiet:
        print(f"Loading DIA-NN results for {args.accession}...")
    dn_df = load_diann(args.accession, args.base_dir, pep_threshold=pep_thresh)

    sg_df = None
    if args.include_sage:
        sage_path = args.base_dir / args.accession / "sage" / "results.sage.parquet"
        if sage_path.exists():
            if not args.quiet:
                print(f"Loading Sage results for {args.accession}...")
            sg_df = load_sage(args.accession, args.base_dir, pep_threshold=pep_thresh)
        else:
            if not args.quiet:
                print(f"Warning: Sage results not found at {sage_path}, skipping 3-way comparison")

    # Compute statistics
    if sg_df is not None:
        stats = compute_threeway_overlap_stats(
            fp_df, dn_df, sg_df,
            match_charge=not args.no_charge,
            use_modifications=args.with_mods
        )
        analysis = analyze_differences(fp_df, dn_df, sg_df)
        if not args.quiet:
            print_threeway_summary(args.accession, stats, analysis)
    else:
        stats = compute_overlap_stats(
            fp_df, dn_df,
            match_charge=not args.no_charge,
            use_modifications=args.with_mods
        )
        analysis = analyze_differences(fp_df, dn_df)
        if not args.quiet:
            print_summary(args.accession, stats, analysis)

    # Save JSON output
    output = {
        "accession": args.accession,
        "timestamp": datetime.now().isoformat(),
        "stats": stats,
        "analysis": analysis,
    }

    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        with open(args.output, "w") as f:
            json.dump(output, f, indent=2, default=lambda x: float(x) if hasattr(x, 'item') else x)
        if not args.quiet:
            print(f"Results saved to: {args.output}")

    # Generate HTML report
    if args.html:
        html_path = args.base_dir / args.accession / "consensus" / "overlap_report.html"
        generate_html_report(args.accession, stats, analysis, html_path, pep_threshold=pep_thresh)
        if not args.quiet:
            print(f"HTML report saved to: {html_path}")

    # Generate consensus parquet files
    if args.consensus:
        if not args.quiet:
            print("Generating consensus DataFrames...")

        overlap_df, union_df = create_consensus_dataframes(
            fp_df, dn_df,
            match_charge=not args.no_charge
        )

        consensus_dir = args.base_dir / args.accession / "consensus"
        consensus_dir.mkdir(parents=True, exist_ok=True)

        overlap_path = consensus_dir / "overlap.parquet"
        union_path = consensus_dir / "union.parquet"

        overlap_df.to_parquet(overlap_path, index=False)
        union_df.to_parquet(union_path, index=False)

        if not args.quiet:
            print(f"Overlap saved to: {overlap_path} ({len(overlap_df):,} precursors)")
            print(f"Union saved to: {union_path} ({len(union_df):,} precursors)")

    return output


if __name__ == "__main__":
    main()
