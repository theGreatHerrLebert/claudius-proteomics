"""
Sequence-Based Precursor Matching

Matches precursors by normalized sequence + charge.
This is the highest-confidence matching strategy (Tier 1/2).
"""

from typing import Dict, List, Tuple, Set, Optional
import pandas as pd

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))
from sequence_utils import remove_modifications

from .config import MatchConfig, MatchTier


def normalize_sequence_for_matching(sequence: str, normalize_il: bool = True) -> str:
    """Normalize sequence for matching purposes.

    1. Remove modifications
    2. Uppercase
    3. Optionally normalize I/L (isoleucine/leucine are isobaric)

    Args:
        sequence: Peptide sequence (may contain modifications)
        normalize_il: Whether to replace I with L

    Returns:
        Normalized plain sequence
    """
    if pd.isna(sequence) or not sequence:
        return ""

    plain = remove_modifications(str(sequence))
    upper = plain.upper()

    if normalize_il:
        return upper.replace("I", "L")
    return upper


class SequenceMatcher:
    """Matches precursors by sequence + charge using hash-based lookup.

    This provides O(1) lookup for sequence-based matching, which is the
    highest confidence matching tier.
    """

    def __init__(self, config: Optional[MatchConfig] = None):
        self.config = config or MatchConfig()

    def create_index(
        self,
        df: pd.DataFrame,
        sequence_col: str,
        charge_col: str,
        raw_file_col: str = "raw_file",
    ) -> Dict[Tuple[str, int, str], List[int]]:
        """Create hash index for sequence-based matching.

        The index maps (normalized_sequence, charge, raw_file) -> list of row indices.
        This enables O(1) lookup for sequence matching.

        Args:
            df: DataFrame to index
            sequence_col: Column containing sequence (may have modifications)
            charge_col: Column containing charge
            raw_file_col: Column containing raw file name

        Returns:
            Dict mapping (norm_sequence, charge, raw_file) -> list of row indices
        """
        index: Dict[Tuple[str, int, str], List[int]] = {}

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

    def match(
        self,
        source_df: pd.DataFrame,
        target_df: pd.DataFrame,
        source_sequence_col: str,
        source_charge_col: str,
        target_sequence_col: str,
        target_charge_col: str,
        raw_file_col: str = "raw_file",
    ) -> List[Tuple[int, int, MatchTier]]:
        """Match source precursors to target by sequence + charge.

        Args:
            source_df: Source DataFrame (e.g., FragPipe precursors)
            target_df: Target DataFrame (e.g., DIA-NN results)
            source_sequence_col: Sequence column in source
            source_charge_col: Charge column in source
            target_sequence_col: Sequence column in target
            target_charge_col: Charge column in target
            raw_file_col: Raw file column (must exist in both)

        Returns:
            List of (source_idx, target_idx, tier) tuples for successful matches
        """
        # Build target index
        target_index = self.create_index(
            target_df,
            sequence_col=target_sequence_col,
            charge_col=target_charge_col,
            raw_file_col=raw_file_col,
        )

        matches = []
        matched_target_indices: Set[int] = set()

        for source_idx, source_row in source_df.iterrows():
            source_seq = source_row.get(source_sequence_col)
            source_charge = source_row.get(source_charge_col)
            raw_file = source_row.get(raw_file_col, "")

            if pd.isna(source_seq) or pd.isna(source_charge):
                continue

            norm_seq = normalize_sequence_for_matching(source_seq, self.config.normalize_il)
            if not norm_seq:
                continue

            key = (norm_seq, int(source_charge), str(raw_file))

            if key in target_index:
                # Find first unmatched target
                for target_idx in target_index[key]:
                    if target_idx not in matched_target_indices:
                        tier = (MatchTier.SEQUENCE_IL_NORM if self.config.normalize_il
                                else MatchTier.SEQUENCE_EXACT)
                        matches.append((source_idx, target_idx, tier))
                        matched_target_indices.add(target_idx)
                        break  # One target per source

        return matches

    def match_to_index(
        self,
        source_row: pd.Series,
        target_index: Dict[Tuple[str, int, str], List[int]],
        source_sequence_col: str,
        source_charge_col: str,
        raw_file_col: str = "raw_file",
        matched_indices: Optional[Set[int]] = None,
    ) -> Optional[Tuple[int, MatchTier]]:
        """Match a single source row against a pre-built target index.

        Useful when matching incrementally or when the target index
        is reused across multiple sources.

        Args:
            source_row: Single row from source DataFrame
            target_index: Pre-built target index
            source_sequence_col: Sequence column in source
            source_charge_col: Charge column in source
            raw_file_col: Raw file column
            matched_indices: Set of already-matched target indices (optional)

        Returns:
            (target_idx, tier) if matched, None otherwise
        """
        matched_indices = matched_indices or set()

        source_seq = source_row.get(source_sequence_col)
        source_charge = source_row.get(source_charge_col)
        raw_file = source_row.get(raw_file_col, "")

        if pd.isna(source_seq) or pd.isna(source_charge):
            return None

        norm_seq = normalize_sequence_for_matching(source_seq, self.config.normalize_il)
        if not norm_seq:
            return None

        key = (norm_seq, int(source_charge), str(raw_file))

        if key in target_index:
            for target_idx in target_index[key]:
                if target_idx not in matched_indices:
                    tier = (MatchTier.SEQUENCE_IL_NORM if self.config.normalize_il
                            else MatchTier.SEQUENCE_EXACT)
                    return target_idx, tier

        return None
