"""
FragPipe-Anchored Precursor Merger

Implements the 3-step merging strategy:
1. Direct FragPipe → Raw join (by precursor_id)
2. Sequence match DIA-NN/Sage to FragPipe
3. Coordinate fallback for unmatched entries

FragPipe is the anchor because it provides direct precursor_id mapping
from the Spectrum column, enabling exact joins to raw timsTOF data.
"""

from pathlib import Path
from typing import Optional, Dict, Set, List
import numpy as np
import pandas as pd

from .config import MatchConfig, MatchTier
from .sequence_matcher import SequenceMatcher, normalize_sequence_for_matching
from .coordinate_matcher import CoordinateMatcher
from .consensus import calculate_consensus, normalize_peptide_columns


class FragPipeAnchoredMerger:
    """Orchestrates merging of search engine results with raw precursor data.

    Strategy:
    1. Start from raw precursors (all fragmented)
    2. Join FragPipe directly by precursor_id (exact match)
    3. Match DIA-NN to merged data by sequence+charge, then coordinates
    4. Match Sage to merged data by sequence+charge, then coordinates
    5. Calculate consensus metrics
    """

    def __init__(self, config: Optional[MatchConfig] = None):
        self.config = config or MatchConfig()
        self.sequence_matcher = SequenceMatcher(self.config)
        self.coordinate_matcher = CoordinateMatcher(self.config)

        # Track match statistics
        self.stats: Dict[str, Dict[str, int]] = {}

    def merge(
        self,
        raw_precursors: pd.DataFrame,
        fragpipe_df: pd.DataFrame,
        diann_df: pd.DataFrame,
        sage_df: pd.DataFrame,
    ) -> pd.DataFrame:
        """Merge all sources into unified precursor index.

        Args:
            raw_precursors: Raw precursor data from timsTOF
            fragpipe_df: FragPipe PSM results
            diann_df: DIA-NN results
            sage_df: Sage results

        Returns:
            Merged DataFrame with all engine results
        """
        print(f"  Merge config: {self.config}")
        self.stats = {}

        # Step 1: Start from raw precursors, join FragPipe by precursor_id
        merged = self._join_fragpipe_direct(raw_precursors, fragpipe_df)
        print(f"  After FragPipe join: {(merged['fragpipe_peptide'].notna()).sum()} identified")

        # Step 2: Join Sage directly by precursor_id (scannr = precursor_id)
        if not sage_df.empty:
            merged = self._join_sage_direct(merged, sage_df)
            print(f"  After Sage join: {(merged['sage_peptide'].notna()).sum()} identified")

        # Step 3: Match DIA-NN using tiered strategy (no direct precursor_id)
        if not diann_df.empty:
            merged = self._match_engine(
                merged, diann_df,
                engine_name="diann",
                peptide_col="diann_peptide",
                modified_col="diann_modified",
                charge_col="diann_charge",
                mz_col="diann_mz",
                rt_col="diann_rt",
                im_col="diann_mobility",
            )

        # Step 4: Calculate consensus
        print("\n  Calculating consensus metrics...")
        merged = normalize_peptide_columns(merged)
        merged = calculate_consensus(merged)

        # Print statistics
        self._print_stats()

        return merged

    def _join_fragpipe_direct(
        self,
        raw_df: pd.DataFrame,
        fragpipe_df: pd.DataFrame,
    ) -> pd.DataFrame:
        """Join FragPipe results directly by precursor_id.

        FragPipe provides precursor_id from the Spectrum column, enabling
        exact join to raw timsTOF data.

        Args:
            raw_df: Raw precursor data
            fragpipe_df: FragPipe PSM results

        Returns:
            Merged DataFrame
        """
        if raw_df.empty:
            # Fall back to FragPipe-only
            return fragpipe_df.copy()

        if fragpipe_df.empty:
            return raw_df.copy()

        # Direct join on (raw_file, precursor_id)
        merged = raw_df.merge(
            fragpipe_df,
            on=["raw_file", "precursor_id"],
            how="left"
        )

        # Track stats
        self.stats["fragpipe"] = {
            "direct_join": (merged["fragpipe_peptide"].notna()).sum(),
        }

        return merged

    def _join_sage_direct(
        self,
        merged_df: pd.DataFrame,
        sage_df: pd.DataFrame,
    ) -> pd.DataFrame:
        """Join Sage results directly by precursor_id.

        Sage's scannr field corresponds to timsTOF precursor_id, enabling
        direct join to raw data (confirmed: 100% intersection with
        pasef_meta_data.precursor_id).

        Args:
            merged_df: Current merged DataFrame (with raw + FragPipe)
            sage_df: Sage PSM results (with precursor_id from scannr)

        Returns:
            Merged DataFrame with Sage columns added
        """
        if sage_df.empty:
            return merged_df

        # Direct join on (raw_file, precursor_id)
        # Sage df has precursor_id column (from scannr)
        sage_cols = [c for c in sage_df.columns if c.startswith("sage_")]
        join_cols = ["raw_file", "precursor_id"]
        sage_subset = sage_df[join_cols + sage_cols].copy()

        merged = merged_df.merge(
            sage_subset,
            on=join_cols,
            how="left"
        )

        # Track stats
        self.stats["sage"] = {
            "direct_join": (merged["sage_peptide"].notna()).sum(),
        }

        return merged

    def _match_engine(
        self,
        merged_df: pd.DataFrame,
        engine_df: pd.DataFrame,
        engine_name: str,
        peptide_col: str,
        modified_col: str,
        charge_col: str,
        mz_col: str,
        rt_col: str,
        im_col: str,
    ) -> pd.DataFrame:
        """Match engine results using tiered strategy.

        Strategy:
        1. Match precursors WITH existing IDs by sequence + charge (highest confidence)
        2. Match remaining by coordinates (m/z + RT + IM)

        Args:
            merged_df: Current merged DataFrame
            engine_df: Engine results to match
            engine_name: Engine name (for column prefixes and stats)
            peptide_col: Plain peptide column
            modified_col: Modified sequence column
            charge_col: Charge column
            mz_col: m/z column
            rt_col: RT column (in seconds)
            im_col: Ion mobility column

        Returns:
            Merged DataFrame with engine results added
        """
        print(f"\n  Matching {engine_name} ({len(engine_df)} entries)...")

        result = merged_df.copy()
        engine_stats = {"sequence": 0, "coordinate": 0}
        coord_rt_diffs: List[float] = []  # Track RT diffs for diagnostics

        # Add raw_file column to engine_df if not present
        if "raw_file" not in engine_df.columns:
            # Extract from run/filename column
            run_col = f"{engine_name}_run" if f"{engine_name}_run" in engine_df.columns else None
            if run_col is None and "diann_run" in engine_df.columns:
                run_col = "diann_run"

            if run_col:
                engine_df = engine_df.copy()
                engine_df["raw_file"] = engine_df[run_col].apply(
                    lambda x: Path(x).stem.replace(".d", "") if pd.notna(x) else ""
                )
            elif f"{engine_name}_raw_file" in engine_df.columns:
                engine_df = engine_df.copy()
                engine_df["raw_file"] = engine_df[f"{engine_name}_raw_file"]

        # Process per raw file for efficiency
        all_matches = []

        for raw_file in result["raw_file"].unique():
            file_mask = result["raw_file"] == raw_file
            file_df = result[file_mask]
            engine_file_df = engine_df[engine_df["raw_file"] == raw_file]

            if engine_file_df.empty:
                continue

            matched_engine_indices: Set[int] = set()
            file_matches = []

            # === Pass 1: Sequence matching for identified precursors ===
            identified_mask = file_df["fragpipe_peptide"].notna()
            identified_df = file_df[identified_mask]

            # Build sequence index on engine results
            seq_index = self.sequence_matcher.create_index(
                engine_file_df,
                sequence_col=modified_col,
                charge_col=charge_col,
                raw_file_col="raw_file",
            )

            for idx, row in identified_df.iterrows():
                # Get source sequence (prefer fragpipe_modified)
                source_seq = row.get("fragpipe_modified")
                if pd.isna(source_seq):
                    source_seq = row.get("fragpipe_peptide")

                # Get charge (prefer fragpipe, fallback to raw)
                source_charge = row.get("fragpipe_charge")
                if pd.isna(source_charge):
                    source_charge = row.get("raw_charge")

                if pd.isna(source_seq) or pd.isna(source_charge):
                    continue

                # Try sequence match
                match = self.sequence_matcher.match_to_index(
                    source_row=pd.Series({
                        modified_col: source_seq,
                        charge_col: source_charge,
                        "raw_file": raw_file,
                    }),
                    target_index=seq_index,
                    source_sequence_col=modified_col,
                    source_charge_col=charge_col,
                    raw_file_col="raw_file",
                    matched_indices=matched_engine_indices,
                )

                if match is not None:
                    engine_idx, tier = match
                    file_matches.append((idx, engine_idx, tier.name))
                    matched_engine_indices.add(engine_idx)
                    engine_stats["sequence"] += 1

            # === Pass 2: Coordinate matching for unmatched precursors ===
            # Build coordinate index on remaining engine results
            unmatched_engine_df = engine_file_df[~engine_file_df.index.isin(matched_engine_indices)]

            if not unmatched_engine_df.empty:
                coord_index = self.coordinate_matcher.create_index(
                    unmatched_engine_df,
                    mz_col=mz_col,
                    charge_col=charge_col,
                    raw_file_col="raw_file",
                )

                for idx, row in file_df.iterrows():
                    # Skip if already matched via sequence
                    if any(m[0] == idx for m in file_matches):
                        continue

                    # Get source coordinates
                    source_mz = row.get("raw_mz")
                    if pd.isna(source_mz):
                        source_mz = row.get("fragpipe_mz")

                    source_charge = row.get("raw_charge")
                    if pd.isna(source_charge):
                        source_charge = row.get("fragpipe_charge")

                    source_rt = row.get("raw_rt_seconds")
                    if pd.isna(source_rt):
                        source_rt = row.get("fragpipe_rt")

                    source_im = row.get("raw_mobility")
                    if pd.isna(source_im):
                        source_im = row.get("fragpipe_mobility")

                    if pd.isna(source_mz) or pd.isna(source_charge):
                        continue

                    # Try coordinate match
                    match = self.coordinate_matcher.match_single(
                        source_mz=float(source_mz),
                        source_charge=int(source_charge),
                        source_rt=float(source_rt) if pd.notna(source_rt) else None,
                        source_im=float(source_im) if pd.notna(source_im) else None,
                        raw_file=raw_file,
                        target_index=coord_index,
                        target_df=unmatched_engine_df,
                        target_rt_col=rt_col,
                        target_im_col=im_col,
                        matched_indices=matched_engine_indices,
                    )

                    if match is not None:
                        engine_idx, match_result = match
                        file_matches.append((idx, engine_idx, match_result.tier.name))
                        matched_engine_indices.add(engine_idx)
                        engine_stats["coordinate"] += 1
                        if match_result.rt_diff_sec is not None:
                            coord_rt_diffs.append(match_result.rt_diff_sec)

            all_matches.extend(file_matches)

        # Apply matches to result DataFrame
        for source_idx, engine_idx, tier_name in all_matches:
            engine_row = engine_df.loc[engine_idx]
            for col in engine_df.columns:
                if col.startswith(f"{engine_name}_"):
                    result.loc[source_idx, col] = engine_row[col]
            result.loc[source_idx, f"{engine_name}_match_tier"] = tier_name

        # Track stats
        self.stats[engine_name] = engine_stats

        total_matched = engine_stats["sequence"] + engine_stats["coordinate"]
        print(f"    Sequence matches: {engine_stats['sequence']}")
        print(f"    Coordinate matches: {engine_stats['coordinate']}")
        print(f"    Total: {total_matched}")

        # RT diagnostics for coordinate matches
        if coord_rt_diffs:
            rt_arr = np.array(coord_rt_diffs)
            print(f"    RT diff (coordinate matches): "
                  f"median={np.median(rt_arr):.4f}s, "
                  f"mean={np.mean(rt_arr):.4f}s, "
                  f"p90={np.percentile(rt_arr, 90):.4f}s, "
                  f"max={np.max(rt_arr):.4f}s")

        return result

    def _print_stats(self):
        """Print merge statistics summary."""
        print("\n  Merge statistics:")
        for engine, stats in self.stats.items():
            print(f"    {engine}:")
            for key, value in stats.items():
                print(f"      {key}: {value}")
