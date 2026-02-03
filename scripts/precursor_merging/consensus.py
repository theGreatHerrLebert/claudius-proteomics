"""
Consensus Calculation

Calculates consensus metrics across search engines:
- n_engines: Number of engines that identified a precursor
- consensus_peptide: Best peptide from engine agreement
- confidence_weight: Weight based on engine agreement
"""

from typing import Dict, List, Tuple, Optional
import pandas as pd

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))
from sequence_utils import normalize_sequence_il


def count_engines(row: pd.Series) -> int:
    """Count number of engines that identified this precursor.

    Args:
        row: DataFrame row with engine peptide columns

    Returns:
        Number of engines (0-3)
    """
    n = 0
    if pd.notna(row.get("fragpipe_peptide")) and row.get("fragpipe_peptide"):
        n += 1
    if pd.notna(row.get("diann_peptide")) and row.get("diann_peptide"):
        n += 1
    if pd.notna(row.get("sage_peptide")) and row.get("sage_peptide"):
        n += 1
    return n


def get_consensus_peptide(row: pd.Series) -> Tuple[str, float]:
    """Get consensus peptide and confidence weight.

    Strategy:
    1. Collect all identified peptides with their confidence scores
    2. If all engines agree -> weight = 1.0
    3. If disagreement -> weight = n_engines / (n_unique + 1)

    Args:
        row: DataFrame row with engine peptide and score columns

    Returns:
        Tuple of (consensus_peptide, confidence_weight)
    """
    peptides: List[Tuple[str, str, float]] = []

    # Collect peptides with normalized sequences and confidence scores
    if pd.notna(row.get("fragpipe_peptide")) and row.get("fragpipe_peptide"):
        norm_seq = normalize_sequence_il(str(row["fragpipe_peptide"]))
        confidence = row.get("fragpipe_probability", 0)
        if pd.isna(confidence):
            confidence = 0
        peptides.append(("fragpipe", norm_seq, float(confidence)))

    if pd.notna(row.get("diann_peptide")) and row.get("diann_peptide"):
        norm_seq = normalize_sequence_il(str(row["diann_peptide"]))
        qvalue = row.get("diann_qvalue", 1)
        if pd.isna(qvalue):
            qvalue = 1
        confidence = 1 - float(qvalue)  # Convert q-value to confidence
        peptides.append(("diann", norm_seq, confidence))

    if pd.notna(row.get("sage_peptide")) and row.get("sage_peptide"):
        norm_seq = normalize_sequence_il(str(row["sage_peptide"]))
        qvalue = row.get("sage_qvalue", 1)
        if pd.isna(qvalue):
            qvalue = 1
        confidence = 1 - float(qvalue)
        peptides.append(("sage", norm_seq, confidence))

    if not peptides:
        return "", 0.0

    # Check agreement
    unique_seqs = set(p[1] for p in peptides)

    if len(unique_seqs) == 1:
        # All engines agree
        return peptides[0][1], 1.0
    else:
        # Disagreement - return most confident
        best = max(peptides, key=lambda x: x[2])
        weight = len(peptides) / (len(unique_seqs) + 1)
        return best[1], weight


def calculate_consensus(df: pd.DataFrame) -> pd.DataFrame:
    """Calculate consensus metrics for all precursors.

    Adds columns:
    - n_engines: Number of engines identifying precursor
    - consensus_peptide: Best peptide from agreement
    - confidence_weight: Weight based on agreement

    Args:
        df: DataFrame with engine peptide columns

    Returns:
        DataFrame with consensus columns added
    """
    result = df.copy()

    # Count engines
    result["n_engines"] = result.apply(count_engines, axis=1)

    # Get consensus peptide and weight
    consensus = result.apply(get_consensus_peptide, axis=1)
    result["consensus_peptide"] = consensus.apply(lambda x: x[0])
    result["confidence_weight"] = consensus.apply(lambda x: x[1])

    return result


def normalize_peptide_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Add normalized peptide columns for matching.

    Adds {engine}_peptide_norm and {engine}_modified_norm columns
    with I/L normalization applied.

    Args:
        df: DataFrame with engine peptide columns

    Returns:
        DataFrame with normalized columns added
    """
    result = df.copy()

    # Normalize plain peptides
    for col in ["fragpipe_peptide", "sage_peptide", "diann_peptide"]:
        if col in result.columns:
            result[f"{col}_norm"] = result[col].apply(
                lambda x: normalize_sequence_il(str(x)) if pd.notna(x) else ""
            )

    # Normalize modified sequences
    for col in ["fragpipe_modified", "sage_modified", "diann_modified"]:
        if col in result.columns:
            result[f"{col}_norm"] = result[col].apply(
                lambda x: normalize_sequence_il(str(x)) if pd.notna(x) else ""
            )

    return result
