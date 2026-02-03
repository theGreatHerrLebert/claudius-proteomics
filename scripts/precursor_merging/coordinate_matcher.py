"""
Coordinate-Based Precursor Matching

Matches precursors by m/z + charge + RT + IM coordinates.
This is the fallback matching strategy (Tier 3/4) when sequence matching fails.
"""

from typing import Dict, List, Tuple, Set, Optional
import pandas as pd

from .config import MatchConfig, MatchTier, MatchResult


class CoordinateMatcher:
    """Matches precursors by coordinate (m/z, RT, IM) using hash-binned lookup.

    Uses m/z bins with neighbor lookup to enable efficient O(1) candidate
    retrieval, followed by precise tolerance checking.
    """

    def __init__(self, config: Optional[MatchConfig] = None):
        self.config = config or MatchConfig()
        # Bin size: at m/z=500 with 20ppm, tolerance=0.01. Use bin_size=0.01
        self.bin_size = 0.01  # Da

    def create_index(
        self,
        df: pd.DataFrame,
        mz_col: str,
        charge_col: str,
        raw_file_col: str = "raw_file",
    ) -> Dict[Tuple[int, int, str], List[Tuple[int, float]]]:
        """Create hash index for coordinate-based matching.

        The index maps (mz_bin, charge, raw_file) -> list of (row_idx, mz) tuples.
        Includes neighbor bins (offset -1, 0, +1) for edge cases.

        Args:
            df: DataFrame to index
            mz_col: Column containing m/z values
            charge_col: Column containing charge
            raw_file_col: Column containing raw file name

        Returns:
            Dict mapping (mz_bin, charge, raw_file) -> list of (row_idx, mz)
        """
        index: Dict[Tuple[int, int, str], List[Tuple[int, float]]] = {}

        for idx, row in df.iterrows():
            mz = row.get(mz_col)
            charge = row.get(charge_col)
            raw_file = row.get(raw_file_col, "")

            if pd.isna(mz) or pd.isna(charge):
                continue

            mz = float(mz)
            mz_bin = int(mz / self.bin_size)

            # Add to current and neighbor bins for edge-case handling
            for offset in [-1, 0, 1]:
                key = (mz_bin + offset, int(charge), str(raw_file))
                if key not in index:
                    index[key] = []
                index[key].append((idx, mz))

        return index

    def match_single(
        self,
        source_mz: float,
        source_charge: int,
        source_rt: Optional[float],
        source_im: Optional[float],
        raw_file: str,
        target_index: Dict[Tuple[int, int, str], List[Tuple[int, float]]],
        target_df: pd.DataFrame,
        target_rt_col: str,
        target_im_col: str,
        matched_indices: Optional[Set[int]] = None,
    ) -> Optional[Tuple[int, MatchResult]]:
        """Match a single precursor against a pre-built target index.

        Args:
            source_mz: Source precursor m/z
            source_charge: Source precursor charge
            source_rt: Source RT in seconds (optional)
            source_im: Source ion mobility (optional)
            raw_file: Raw file name
            target_index: Pre-built coordinate index
            target_df: Target DataFrame for retrieving RT/IM values
            target_rt_col: RT column in target (should be in seconds)
            target_im_col: IM column in target
            matched_indices: Set of already-matched target indices

        Returns:
            (target_idx, MatchResult) if matched, None otherwise
        """
        matched_indices = matched_indices or set()
        config = self.config

        mz_bin = int(source_mz / self.bin_size)
        key = (mz_bin, source_charge, raw_file)

        if key not in target_index:
            return None

        mz_tol = source_mz * config.mz_tol_ppm / 1e6

        best_match: Optional[Tuple[int, MatchResult]] = None
        best_score = float('inf')

        for target_idx, target_mz in target_index[key]:
            if target_idx in matched_indices:
                continue

            mz_diff = abs(source_mz - target_mz)
            if mz_diff > mz_tol:
                continue

            mz_diff_ppm = mz_diff / source_mz * 1e6
            score = mz_diff_ppm / config.mz_tol_ppm

            # Check RT/IM for tier determination
            has_rt_match = False
            has_im_match = False
            rt_diff_sec = None
            im_diff = None

            target_row = target_df.loc[target_idx]

            target_rt = target_row.get(target_rt_col)
            if source_rt is not None and pd.notna(target_rt):
                rt_diff_sec = abs(float(source_rt) - float(target_rt))
                if rt_diff_sec <= config.rt_tol_sec:
                    has_rt_match = True
                    score += rt_diff_sec / config.rt_tol_sec
                else:
                    # RT mismatch - skip this candidate
                    continue

            target_im = target_row.get(target_im_col)
            if source_im is not None and pd.notna(target_im):
                im_diff = abs(float(source_im) - float(target_im))
                if im_diff <= config.im_tol:
                    has_im_match = True
                    score += im_diff / config.im_tol
                else:
                    # IM mismatch - skip this candidate
                    continue

            if score < best_score:
                best_score = score
                full_match = has_rt_match and has_im_match
                result = MatchResult.coordinate_match(
                    score=score,
                    full=full_match,
                    mz_diff_ppm=mz_diff_ppm,
                    rt_diff_sec=rt_diff_sec,
                    im_diff=im_diff,
                )
                best_match = (target_idx, result)

        return best_match

    def match(
        self,
        source_df: pd.DataFrame,
        target_df: pd.DataFrame,
        source_mz_col: str,
        source_charge_col: str,
        source_rt_col: str,
        source_im_col: str,
        target_mz_col: str,
        target_charge_col: str,
        target_rt_col: str,
        target_im_col: str,
        raw_file_col: str = "raw_file",
        skip_source_indices: Optional[Set[int]] = None,
    ) -> List[Tuple[int, int, MatchResult]]:
        """Match source precursors to target by coordinates.

        Args:
            source_df: Source DataFrame
            target_df: Target DataFrame
            source_mz_col: m/z column in source
            source_charge_col: Charge column in source
            source_rt_col: RT column in source (seconds)
            source_im_col: IM column in source
            target_mz_col: m/z column in target
            target_charge_col: Charge column in target
            target_rt_col: RT column in target (seconds)
            target_im_col: IM column in target
            raw_file_col: Raw file column (must exist in both)
            skip_source_indices: Source indices to skip (already matched)

        Returns:
            List of (source_idx, target_idx, MatchResult) tuples
        """
        skip_source_indices = skip_source_indices or set()

        # Build target index
        target_index = self.create_index(
            target_df,
            mz_col=target_mz_col,
            charge_col=target_charge_col,
            raw_file_col=raw_file_col,
        )

        matches = []
        matched_target_indices: Set[int] = set()

        for source_idx, source_row in source_df.iterrows():
            if source_idx in skip_source_indices:
                continue

            source_mz = source_row.get(source_mz_col)
            source_charge = source_row.get(source_charge_col)
            raw_file = source_row.get(raw_file_col, "")

            if pd.isna(source_mz) or pd.isna(source_charge):
                continue

            source_rt = source_row.get(source_rt_col)
            if pd.isna(source_rt):
                source_rt = None
            else:
                source_rt = float(source_rt)

            source_im = source_row.get(source_im_col)
            if pd.isna(source_im):
                source_im = None
            else:
                source_im = float(source_im)

            match_result = self.match_single(
                source_mz=float(source_mz),
                source_charge=int(source_charge),
                source_rt=source_rt,
                source_im=source_im,
                raw_file=raw_file,
                target_index=target_index,
                target_df=target_df,
                target_rt_col=target_rt_col,
                target_im_col=target_im_col,
                matched_indices=matched_target_indices,
            )

            if match_result is not None:
                target_idx, result = match_result
                matches.append((source_idx, target_idx, result))
                matched_target_indices.add(target_idx)

        return matches
