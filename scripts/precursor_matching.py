#!/usr/bin/env python3
"""
Precursor Matching Rules

Defines the rule-set for merging precursor identifications across search engines.

Strategy:
1. PRIMARY: Match by standardized modified sequence + charge (highest confidence)
2. SECONDARY: Match by m/z + charge + RT + IM (for unidentified or disagreeing precursors)

Matching Tiers:
- Tier 1: Exact sequence match (I/L normalized) + charge + raw_file
- Tier 2: m/z match (within tolerance) + charge + RT + IM + raw_file
- Tier 3: m/z match (within tolerance) + charge + raw_file (no RT/IM available)

Output includes match quality scores and disagreement flags.
"""

from dataclasses import dataclass
from typing import Optional, Dict, List, Tuple, Set
import pandas as pd
import numpy as np
from enum import Enum

from sequence_utils import normalize_sequence_il, remove_modifications


class MatchTier(Enum):
    """Match quality tiers."""
    SEQUENCE_EXACT = 1      # Exact sequence + charge match
    SEQUENCE_IL_NORM = 2    # I/L normalized sequence + charge match
    COORDINATE_FULL = 3     # m/z + charge + RT + IM match
    COORDINATE_PARTIAL = 4  # m/z + charge match (missing RT/IM)
    NO_MATCH = 5            # No match found


@dataclass
class MatchConfig:
    """Configuration for precursor matching."""
    # m/z tolerance
    mz_tol_ppm: float = 20.0

    # RT tolerance (seconds)
    rt_tol_sec: float = 30.0

    # Ion mobility tolerance (1/K0)
    im_tol: float = 0.05

    # Charge must match exactly
    charge_must_match: bool = True

    # Whether to normalize I/L in sequences
    normalize_il: bool = True


@dataclass
class MatchResult:
    """Result of a precursor match attempt."""
    matched: bool
    tier: MatchTier
    score: float  # Lower is better

    # Differences (for diagnostics)
    mz_diff_ppm: Optional[float] = None
    rt_diff_sec: Optional[float] = None
    im_diff: Optional[float] = None
    sequence_match: bool = False
    charge_match: bool = False


def normalize_sequence_for_matching(sequence: str, normalize_il: bool = True) -> str:
    """Normalize sequence for matching purposes.

    1. Remove modifications
    2. Uppercase
    3. Optionally normalize I/L
    """
    if pd.isna(sequence) or not sequence:
        return ""

    plain = remove_modifications(str(sequence))
    upper = plain.upper()

    if normalize_il:
        return upper.replace("I", "L")
    return upper


def compute_mz_diff_ppm(mz1: float, mz2: float) -> float:
    """Compute m/z difference in ppm."""
    if mz1 == 0 or mz2 == 0:
        return float('inf')
    return abs(mz1 - mz2) / mz1 * 1e6


def match_precursors(
    source_row: Dict,
    target_row: Dict,
    config: MatchConfig,
    source_prefix: str = "",
    target_prefix: str = "",
) -> MatchResult:
    """
    Attempt to match two precursor records.

    Args:
        source_row: Source precursor data (dict-like)
        target_row: Target precursor data (dict-like)
        config: Matching configuration
        source_prefix: Prefix for source columns (e.g., "fragpipe_")
        target_prefix: Prefix for target columns (e.g., "sage_")

    Returns:
        MatchResult with match status and quality metrics
    """

    def get_val(row, key, prefix=""):
        """Get value from row with optional prefix."""
        full_key = f"{prefix}{key}" if prefix else key
        val = row.get(full_key)
        if pd.isna(val):
            return None
        return val

    # Extract values
    source_seq = get_val(source_row, "modified", source_prefix) or get_val(source_row, "peptide", source_prefix)
    target_seq = get_val(target_row, "modified", target_prefix) or get_val(target_row, "peptide", target_prefix)

    source_charge = get_val(source_row, "charge", source_prefix)
    target_charge = get_val(target_row, "charge", target_prefix)

    source_mz = get_val(source_row, "mz", source_prefix)
    target_mz = get_val(target_row, "mz", target_prefix)

    source_rt = get_val(source_row, "rt", source_prefix)
    target_rt = get_val(target_row, "rt", target_prefix)

    source_im = get_val(source_row, "mobility", source_prefix)
    target_im = get_val(target_row, "mobility", target_prefix)

    # Check charge match
    charge_match = False
    if source_charge is not None and target_charge is not None:
        charge_match = int(source_charge) == int(target_charge)

    if config.charge_must_match and not charge_match:
        return MatchResult(
            matched=False,
            tier=MatchTier.NO_MATCH,
            score=float('inf'),
            charge_match=False
        )

    # Normalize sequences
    source_seq_norm = normalize_sequence_for_matching(source_seq, config.normalize_il)
    target_seq_norm = normalize_sequence_for_matching(target_seq, config.normalize_il)

    # Check sequence match
    sequence_match = False
    if source_seq_norm and target_seq_norm:
        sequence_match = source_seq_norm == target_seq_norm

    # Compute coordinate differences
    mz_diff_ppm = None
    rt_diff_sec = None
    im_diff = None

    if source_mz is not None and target_mz is not None:
        mz_diff_ppm = compute_mz_diff_ppm(float(source_mz), float(target_mz))

    if source_rt is not None and target_rt is not None:
        rt_diff_sec = abs(float(source_rt) - float(target_rt))

    if source_im is not None and target_im is not None:
        im_diff = abs(float(source_im) - float(target_im))

    # Determine match tier and score
    if sequence_match and charge_match:
        # Tier 1/2: Sequence match
        tier = MatchTier.SEQUENCE_EXACT if not config.normalize_il else MatchTier.SEQUENCE_IL_NORM
        score = 0.0  # Perfect match

        return MatchResult(
            matched=True,
            tier=tier,
            score=score,
            mz_diff_ppm=mz_diff_ppm,
            rt_diff_sec=rt_diff_sec,
            im_diff=im_diff,
            sequence_match=True,
            charge_match=True
        )

    # Check coordinate-based match
    mz_within_tol = mz_diff_ppm is not None and mz_diff_ppm <= config.mz_tol_ppm
    rt_within_tol = rt_diff_sec is not None and rt_diff_sec <= config.rt_tol_sec
    im_within_tol = im_diff is not None and im_diff <= config.im_tol

    if mz_within_tol and charge_match:
        if rt_within_tol and im_within_tol:
            # Tier 3: Full coordinate match
            score = (mz_diff_ppm / config.mz_tol_ppm +
                    rt_diff_sec / config.rt_tol_sec +
                    im_diff / config.im_tol) / 3

            return MatchResult(
                matched=True,
                tier=MatchTier.COORDINATE_FULL,
                score=score,
                mz_diff_ppm=mz_diff_ppm,
                rt_diff_sec=rt_diff_sec,
                im_diff=im_diff,
                sequence_match=sequence_match,
                charge_match=True
            )
        elif rt_within_tol or im_within_tol:
            # Partial coordinate match
            score = mz_diff_ppm / config.mz_tol_ppm
            if rt_within_tol:
                score += rt_diff_sec / config.rt_tol_sec
            if im_within_tol:
                score += im_diff / config.im_tol
            score /= (1 + int(rt_within_tol) + int(im_within_tol))

            return MatchResult(
                matched=True,
                tier=MatchTier.COORDINATE_PARTIAL,
                score=score + 1.0,  # Penalty for partial
                mz_diff_ppm=mz_diff_ppm,
                rt_diff_sec=rt_diff_sec,
                im_diff=im_diff,
                sequence_match=sequence_match,
                charge_match=True
            )
        else:
            # m/z + charge only
            score = mz_diff_ppm / config.mz_tol_ppm + 2.0  # Penalty

            return MatchResult(
                matched=True,
                tier=MatchTier.COORDINATE_PARTIAL,
                score=score,
                mz_diff_ppm=mz_diff_ppm,
                rt_diff_sec=rt_diff_sec,
                im_diff=im_diff,
                sequence_match=sequence_match,
                charge_match=True
            )

    # No match
    return MatchResult(
        matched=False,
        tier=MatchTier.NO_MATCH,
        score=float('inf'),
        mz_diff_ppm=mz_diff_ppm,
        rt_diff_sec=rt_diff_sec,
        im_diff=im_diff,
        sequence_match=sequence_match,
        charge_match=charge_match
    )


class PrecursorMatcher:
    """
    Matches precursors across search engines using a tiered strategy.

    Strategy:
    1. Build hash index on (normalized_sequence, charge, raw_file)
    2. For sequence matches, directly link precursors
    3. For non-matches, fall back to coordinate-based matching
    """

    def __init__(self, config: Optional[MatchConfig] = None):
        self.config = config or MatchConfig()

    def create_sequence_index(
        self,
        df: pd.DataFrame,
        sequence_col: str,
        charge_col: str,
        raw_file_col: str = "raw_file",
    ) -> Dict[Tuple[str, int, str], List[int]]:
        """
        Create hash index for sequence-based matching.

        Returns:
            Dict mapping (norm_sequence, charge, raw_file) -> list of row indices
        """
        index = {}

        for idx, row in df.iterrows():
            seq = row.get(sequence_col)
            charge = row.get(charge_col)
            raw_file = row.get(raw_file_col, "")

            if pd.isna(seq) or pd.isna(charge):
                continue

            norm_seq = normalize_sequence_for_matching(seq, self.config.normalize_il)
            if not norm_seq:
                continue

            key = (norm_seq, int(charge), str(raw_file))
            if key not in index:
                index[key] = []
            index[key].append(idx)

        return index

    def create_coordinate_index(
        self,
        df: pd.DataFrame,
        mz_col: str,
        charge_col: str,
        raw_file_col: str = "raw_file",
        bin_size: float = 0.01,  # Da
    ) -> Dict[Tuple[int, int, str], List[Tuple[int, float]]]:
        """
        Create hash index for coordinate-based matching.

        Returns:
            Dict mapping (mz_bin, charge, raw_file) -> list of (row_idx, mz)
        """
        index = {}

        for idx, row in df.iterrows():
            mz = row.get(mz_col)
            charge = row.get(charge_col)
            raw_file = row.get(raw_file_col, "")

            if pd.isna(mz) or pd.isna(charge):
                continue

            mz = float(mz)
            mz_bin = int(mz / bin_size)

            # Add to current and neighbor bins
            for offset in [-1, 0, 1]:
                key = (mz_bin + offset, int(charge), str(raw_file))
                if key not in index:
                    index[key] = []
                index[key].append((idx, mz))

        return index

    def match_dataframes(
        self,
        source_df: pd.DataFrame,
        target_df: pd.DataFrame,
        source_cols: Dict[str, str],
        target_cols: Dict[str, str],
        raw_file_col: str = "raw_file",
    ) -> pd.DataFrame:
        """
        Match precursors between two dataframes.

        Args:
            source_df: Source dataframe (e.g., FragPipe)
            target_df: Target dataframe (e.g., Sage)
            source_cols: Column mapping for source {
                'sequence': col_name,
                'modified': col_name,
                'charge': col_name,
                'mz': col_name,
                'rt': col_name,
                'mobility': col_name,
            }
            target_cols: Column mapping for target (same structure)
            raw_file_col: Column name for raw file

        Returns:
            DataFrame with match results
        """
        results = []

        # Build indices
        target_seq_index = self.create_sequence_index(
            target_df,
            sequence_col=target_cols.get('modified', target_cols.get('sequence')),
            charge_col=target_cols['charge'],
            raw_file_col=raw_file_col,
        )

        target_coord_index = self.create_coordinate_index(
            target_df,
            mz_col=target_cols['mz'],
            charge_col=target_cols['charge'],
            raw_file_col=raw_file_col,
        )

        matched_target_indices = set()

        # Match each source precursor
        for source_idx, source_row in source_df.iterrows():
            raw_file = source_row.get(raw_file_col, "")

            # Try sequence-based match first
            source_seq = source_row.get(source_cols.get('modified', source_cols.get('sequence')))
            source_charge = source_row.get(source_cols['charge'])

            best_match = None
            best_target_idx = None

            if not pd.isna(source_seq) and not pd.isna(source_charge):
                norm_seq = normalize_sequence_for_matching(source_seq, self.config.normalize_il)
                key = (norm_seq, int(source_charge), str(raw_file))

                if key in target_seq_index:
                    # Found sequence match - pick best one
                    for target_idx in target_seq_index[key]:
                        if target_idx in matched_target_indices:
                            continue

                        target_row = target_df.loc[target_idx]
                        match_result = match_precursors(
                            dict(source_row),
                            dict(target_row),
                            self.config,
                            source_prefix="",
                            target_prefix="",
                        )

                        if match_result.matched:
                            if best_match is None or match_result.score < best_match.score:
                                best_match = match_result
                                best_target_idx = target_idx

            # If no sequence match, try coordinate-based
            if best_match is None:
                source_mz = source_row.get(source_cols['mz'])

                if not pd.isna(source_mz) and not pd.isna(source_charge):
                    mz_bin = int(float(source_mz) / 0.01)
                    key = (mz_bin, int(source_charge), str(raw_file))

                    if key in target_coord_index:
                        for target_idx, target_mz in target_coord_index[key]:
                            if target_idx in matched_target_indices:
                                continue

                            target_row = target_df.loc[target_idx]
                            match_result = match_precursors(
                                dict(source_row),
                                dict(target_row),
                                self.config,
                            )

                            if match_result.matched:
                                if best_match is None or match_result.score < best_match.score:
                                    best_match = match_result
                                    best_target_idx = target_idx

            # Record result
            result = {
                'source_idx': source_idx,
                'target_idx': best_target_idx,
                'matched': best_match is not None and best_match.matched,
                'match_tier': best_match.tier.name if best_match else 'NO_MATCH',
                'match_score': best_match.score if best_match else float('inf'),
                'sequence_match': best_match.sequence_match if best_match else False,
                'mz_diff_ppm': best_match.mz_diff_ppm if best_match else None,
                'rt_diff_sec': best_match.rt_diff_sec if best_match else None,
                'im_diff': best_match.im_diff if best_match else None,
            }
            results.append(result)

            if best_target_idx is not None:
                matched_target_indices.add(best_target_idx)

        return pd.DataFrame(results)


def merge_engine_results(
    base_df: pd.DataFrame,
    engine_dfs: Dict[str, pd.DataFrame],
    engine_cols: Dict[str, Dict[str, str]],
    config: Optional[MatchConfig] = None,
    raw_file_col: str = "raw_file",
    precursor_id_col: str = "precursor_id",
) -> pd.DataFrame:
    """
    Merge results from multiple search engines into a unified precursor index.

    Args:
        base_df: Base dataframe with raw precursor info
        engine_dfs: Dict of engine name -> results dataframe
        engine_cols: Dict of engine name -> column mapping
        config: Matching configuration
        raw_file_col: Column name for raw file
        precursor_id_col: Column name for precursor ID

    Returns:
        Merged dataframe with all engine results
    """
    config = config or MatchConfig()
    matcher = PrecursorMatcher(config)

    result_df = base_df.copy()

    for engine_name, engine_df in engine_dfs.items():
        cols = engine_cols[engine_name]
        print(f"\nMatching {engine_name}...")

        # Match engine results to base precursors
        match_results = matcher.match_dataframes(
            source_df=result_df,
            target_df=engine_df,
            source_cols={
                'sequence': 'consensus_peptide',
                'modified': 'fragpipe_modified',  # Use first engine's modified as reference
                'charge': 'raw_charge',
                'mz': 'raw_mz',
                'rt': 'fragpipe_rt',
                'mobility': 'raw_mobility',
            },
            target_cols=cols,
            raw_file_col=raw_file_col,
        )

        # Report match statistics
        tier_counts = match_results['match_tier'].value_counts()
        print(f"  Match tiers:")
        for tier, count in tier_counts.items():
            print(f"    {tier}: {count}")

        # Merge matched results
        matched = match_results[match_results['matched']]

        for _, match_row in matched.iterrows():
            source_idx = match_row['source_idx']
            target_idx = match_row['target_idx']

            if target_idx is not None:
                target_row = engine_df.loc[target_idx]

                # Add engine columns to result
                for col in engine_df.columns:
                    if col.startswith(f"{engine_name}_"):
                        result_df.loc[source_idx, col] = target_row[col]

                # Add match quality info
                result_df.loc[source_idx, f"{engine_name}_match_tier"] = match_row['match_tier']
                result_df.loc[source_idx, f"{engine_name}_match_score"] = match_row['match_score']

    return result_df


if __name__ == "__main__":
    # Test the matching logic
    print("Testing precursor matching rules...")

    config = MatchConfig()

    # Test case 1: Exact sequence match
    source = {
        'fragpipe_modified': 'PEPTIDE[UNIMOD:4]K',
        'fragpipe_charge': 2,
        'fragpipe_mz': 500.0,
        'fragpipe_rt': 100.0,
        'fragpipe_mobility': 1.0,
    }
    target = {
        'sage_modified': 'PEPTIDE[UNIMOD:4]K',
        'sage_charge': 2,
        'sage_mz': 500.001,
        'sage_rt': 100.5,
        'sage_mobility': 1.01,
    }

    result = match_precursors(source, target, config, 'fragpipe_', 'sage_')
    print(f"\nTest 1 - Exact sequence match:")
    print(f"  Matched: {result.matched}, Tier: {result.tier.name}, Score: {result.score:.3f}")

    # Test case 2: I/L difference
    source2 = {
        'fragpipe_modified': 'PEPTIDEK',
        'fragpipe_charge': 2,
        'fragpipe_mz': 500.0,
    }
    target2 = {
        'sage_modified': 'PEPTLDEK',
        'sage_charge': 2,
        'sage_mz': 500.0,
    }

    result2 = match_precursors(source2, target2, config, 'fragpipe_', 'sage_')
    print(f"\nTest 2 - I/L normalized match:")
    print(f"  Matched: {result2.matched}, Tier: {result2.tier.name}")

    # Test case 3: Coordinate-only match
    source3 = {
        'fragpipe_modified': 'AAAAAAA',
        'fragpipe_charge': 2,
        'fragpipe_mz': 500.0,
        'fragpipe_rt': 100.0,
        'fragpipe_mobility': 1.0,
    }
    target3 = {
        'sage_modified': 'BBBBBBB',
        'sage_charge': 2,
        'sage_mz': 500.005,
        'sage_rt': 100.2,
        'sage_mobility': 1.02,
    }

    result3 = match_precursors(source3, target3, config, 'fragpipe_', 'sage_')
    print(f"\nTest 3 - Coordinate match (different sequences):")
    print(f"  Matched: {result3.matched}, Tier: {result3.tier.name}, Score: {result3.score:.3f}")
    print(f"  m/z diff: {result3.mz_diff_ppm:.1f} ppm, RT diff: {result3.rt_diff_sec:.1f} sec")

    # Test case 4: No match (charge mismatch)
    source4 = {
        'fragpipe_modified': 'PEPTIDEK',
        'fragpipe_charge': 2,
        'fragpipe_mz': 500.0,
    }
    target4 = {
        'sage_modified': 'PEPTIDEK',
        'sage_charge': 3,
        'sage_mz': 500.0,
    }

    result4 = match_precursors(source4, target4, config, 'fragpipe_', 'sage_')
    print(f"\nTest 4 - Charge mismatch:")
    print(f"  Matched: {result4.matched}, Tier: {result4.tier.name}")

    print("\nAll tests complete!")
