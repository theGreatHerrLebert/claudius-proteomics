#!/usr/bin/env python3
"""
Step 3: Stratify and Merge Third-Party Output

UNIMOD-standardizes sequences, computes consensus, stratifies by engine agreement.

Uses FragPipe-anchored merging strategy:
1. Start from raw precursors (from Step 4 raw_features.parquet)
2. Direct join FragPipe by (raw_file, precursor_id) - FragPipe has this!
3. Match DIA-NN/Sage by sequence+charge, fallback to coordinates
4. Output precursor_index with precursor_id enabling clean Step 5 merge

Input: FragPipe/DIA-NN/Sage outputs from Step 2

Outputs:
- data/processed/{accession}/precursor_index.parquet
- data/processed/{accession}/consensus/
    - overlap_stats.json
    - overlap_report.html
    - stratified/
        - all_three.parquet
        - two_plus.parquet
        - fragpipe_only.parquet
        - diann_only.parquet
        - sage_only.parquet
- step3_summary.json
"""

import json
import sys
from pathlib import Path
from typing import Dict, Any, Optional, Set, Tuple

import pandas as pd
import numpy as np

# Add parent to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from lib.sequence_utils import (
    standardize_diann_sequence,
    standardize_sage_sequence,
    standardize_fragpipe_modified_peptide,
    normalize_sequence_il,
)
from lib.precursor_matching import MatchConfig, PrecursorMatcher, MatchTier
from runner.summary import StepSummary, write_step_summary
from scripts.engine_parsers.fragpipe_parser import FragPipeParser
from scripts.engine_parsers.diann_parser import DiannParser
from scripts.engine_parsers.sage_parser import SageParser


def run_step3_stratify(
    accession: str,
    config: Dict[str, Any],
    output_base_dir: Path,
    generate_html: bool = True,
) -> StepSummary:
    """
    Execute Step 3: Stratify and merge search engine results.

    Uses FragPipe-anchored merging:
    1. Load raw precursors from Step 4 raw_features.parquet
    2. Direct join FragPipe by (raw_file, precursor_id)
    3. Match DIA-NN/Sage by sequence+charge, fallback to coordinates
    4. Output precursor_index with precursor_id for Step 5 merge

    Args:
        accession: PRIDE accession
        config: Pipeline configuration dict
        output_base_dir: Base directory for outputs
        generate_html: Whether to generate HTML report

    Returns:
        StepSummary with results
    """
    summary = StepSummary(
        step_name="step3",
        accession=accession,
    )

    processed_dir = output_base_dir / "processed" / accession
    extracted_dir = output_base_dir / "extracted" / accession
    consensus_dir = processed_dir / "consensus"
    stratified_dir = consensus_dir / "stratified"

    try:
        # Load raw precursors from Step 4 (provides the anchor)
        print("  Loading raw precursors from Step 4...")
        raw_features_path = extracted_dir / "raw_features.parquet"
        if raw_features_path.exists():
            raw_precursors = pd.read_parquet(raw_features_path)
            # Keep only columns needed for anchoring
            anchor_cols = ['precursor_id', 'raw_file', 'mz', 'charge', 'rt_seconds', 'mobility']
            available_cols = [c for c in anchor_cols if c in raw_precursors.columns]
            raw_precursors = raw_precursors[available_cols].copy()
            print(f"    {len(raw_precursors)} raw precursors loaded")
        else:
            print("    No raw_features.parquet found - will build index from engine results only")
            raw_precursors = pd.DataFrame()

        # Load search engine results using parsers (preserves precursor_id for FragPipe)
        print("  Loading search engine results...")
        fp_df = _load_fragpipe_with_parser(processed_dir)
        dn_df = _load_diann_with_parser(processed_dir)
        sg_df = _load_sage_with_parser(processed_dir)

        print(f"    FragPipe: {len(fp_df) if fp_df is not None and not fp_df.empty else 0} PSMs")
        print(f"    DIA-NN: {len(dn_df) if dn_df is not None and not dn_df.empty else 0} precursors")
        print(f"    Sage: {len(sg_df) if sg_df is not None and not sg_df.empty else 0} PSMs")

        # Compute overlap statistics (at sequence+charge level for reporting)
        print("  Computing overlap statistics...")
        overlap_stats = _compute_overlap_stats_from_parsers(fp_df, dn_df, sg_df)

        # Build unified precursor index with FragPipe anchoring
        print("  Building unified precursor index (FragPipe-anchored)...")
        precursor_index = _build_precursor_index_anchored(
            raw_precursors, fp_df, dn_df, sg_df, config
        )

        # Stratify by engine agreement
        print("  Stratifying by engine agreement...")
        stratified_counts = _stratify_precursors(precursor_index, stratified_dir)

        # Compute match tier statistics
        match_tiers = _compute_match_tiers(precursor_index)

        # Save precursor index
        index_path = processed_dir / "precursor_index.parquet"
        precursor_index.to_parquet(index_path, index=False)
        print(f"  Saved precursor index: {index_path}")
        print(f"    Columns: {list(precursor_index.columns)}")

        # Save overlap stats
        consensus_dir.mkdir(parents=True, exist_ok=True)
        stats_path = consensus_dir / "overlap_stats.json"
        with open(stats_path, "w") as f:
            json.dump(overlap_stats, f, indent=2)

        # Generate HTML report
        if generate_html:
            html_path = consensus_dir / "overlap_report.html"
            _generate_html_report(accession, overlap_stats, html_path)
            print(f"  Generated HTML report: {html_path}")

        # Update summary
        summary.data = {
            "overlap_stats": overlap_stats,
            "stratified_counts": stratified_counts,
            "match_tiers": match_tiers,
            "n_total_precursors": len(precursor_index),
        }
        summary.outputs = [
            str(index_path),
            str(consensus_dir),
        ]
        summary.complete(success=True)

    except Exception as e:
        summary.complete(success=False, error_message=str(e))
        raise

    # Write summary file
    write_step_summary(summary, processed_dir)

    return summary


def _load_fragpipe_with_parser(processed_dir: Path) -> Optional[pd.DataFrame]:
    """Load FragPipe results using parser that extracts precursor_id.

    FragPipe provides precursor_id from the Spectrum column format:
    rawfile.scannum.scannum.charge - the scannum IS the timsTOF precursor_id.
    """
    parser = FragPipeParser()
    # processed_dir is like data/processed/PXD019086, parser expects base_dir + accession
    base_dir = processed_dir.parent
    accession = processed_dir.name
    return parser.parse(base_dir, accession)


def _load_diann_with_parser(processed_dir: Path) -> Optional[pd.DataFrame]:
    """Load DIA-NN results using parser.

    DIA-NN does not have precursor_id - requires sequence/coordinate matching.
    RT is converted to seconds at parse time.
    """
    parser = DiannParser()
    base_dir = processed_dir.parent
    accession = processed_dir.name
    return parser.parse(base_dir, accession)


def _load_sage_with_parser(processed_dir: Path) -> Optional[pd.DataFrame]:
    """Load Sage results using parser.

    Sage scannr+1 = raw precursor_id (raw extraction uses 1-based indexing).
    m/z is calculated from experimental mass at parse time.
    """
    parser = SageParser()
    base_dir = processed_dir.parent
    accession = processed_dir.name
    return parser.parse(base_dir, accession)


# Keep legacy loaders for backward compatibility (used by _compute_overlap_stats)
def _load_fragpipe(processed_dir: Path) -> Optional[pd.DataFrame]:
    """Load FragPipe results (legacy loader for overlap stats)."""
    # Try combined_ion.tsv first
    combined_ion = processed_dir / "combined_ion.tsv"
    if combined_ion.exists():
        return pd.read_csv(combined_ion, sep="\t")

    # Try psm.tsv files
    psm_files = list(processed_dir.rglob("psm.tsv"))
    if psm_files:
        dfs = [pd.read_csv(f, sep="\t") for f in psm_files]
        return pd.concat(dfs, ignore_index=True)

    return None


def _load_diann(processed_dir: Path) -> Optional[pd.DataFrame]:
    """Load DIA-NN results (legacy loader for overlap stats)."""
    report_path = processed_dir / "diann" / "report.parquet"
    if report_path.exists():
        return pd.read_parquet(report_path)

    report_tsv = processed_dir / "diann" / "report.tsv"
    if report_tsv.exists():
        return pd.read_csv(report_tsv, sep="\t")

    return None


def _load_sage(processed_dir: Path) -> Optional[pd.DataFrame]:
    """Load Sage results (legacy loader for overlap stats)."""
    results_path = processed_dir / "sage" / "results.sage.parquet"
    if results_path.exists():
        df = pd.read_parquet(results_path)
        # Filter decoys
        if "is_decoy" in df.columns:
            df = df[~df["is_decoy"]]
        return df

    return None


def _standardize_fragpipe(df: pd.DataFrame) -> pd.DataFrame:
    """Standardize FragPipe sequences."""
    df = df.copy()

    # Try different column names for modified peptide
    mod_col = None
    for col in ["Modified Sequence", "Modified Peptide", "modified_sequence", "modified_peptide"]:
        if col in df.columns:
            mod_col = col
            break

    # Get bare peptide column for fallback
    bare_col = None
    for col in ["Peptide", "Peptide Sequence"]:
        if col in df.columns:
            bare_col = col
            break

    if mod_col and bare_col:
        # Use Modified Peptide when available, fall back to bare Peptide when NaN/empty
        def standardize_with_fallback(row):
            mod_seq = row[mod_col]
            if pd.notna(mod_seq) and mod_seq:
                result = standardize_fragpipe_modified_peptide(mod_seq)
                if result:
                    return result
            # Fallback to bare peptide
            return row[bare_col] if pd.notna(row[bare_col]) else ""
        df["modified_unimod"] = df.apply(standardize_with_fallback, axis=1)
    elif mod_col:
        df["modified_unimod"] = df[mod_col].apply(standardize_fragpipe_modified_peptide)
    elif bare_col:
        df["modified_unimod"] = df[bare_col]
    else:
        raise ValueError(f"No peptide column found in FragPipe output. Columns: {list(df.columns[:10])}")

    # Create normalized sequence for matching
    df["sequence_normalized"] = df["modified_unimod"].apply(
        lambda x: normalize_sequence_il(x) if pd.notna(x) else ""
    )

    return df


def _standardize_diann(df: pd.DataFrame) -> pd.DataFrame:
    """Standardize DIA-NN sequences."""
    df = df.copy()

    mod_col = "Modified.Sequence" if "Modified.Sequence" in df.columns else "Precursor.Id"
    df["modified_unimod"] = df[mod_col].apply(standardize_diann_sequence)

    df["sequence_normalized"] = df["modified_unimod"].apply(
        lambda x: normalize_sequence_il(x) if pd.notna(x) else ""
    )

    return df


def _standardize_sage(df: pd.DataFrame) -> pd.DataFrame:
    """Standardize Sage sequences."""
    df = df.copy()

    if "peptide" in df.columns:
        df["modified_unimod"] = df["peptide"].apply(standardize_sage_sequence)
    else:
        df["modified_unimod"] = df.get("stripped_peptide", "")

    df["sequence_normalized"] = df["modified_unimod"].apply(
        lambda x: normalize_sequence_il(x) if pd.notna(x) else ""
    )

    return df


def _compute_overlap_stats_from_parsers(
    fp_df: Optional[pd.DataFrame],
    dn_df: Optional[pd.DataFrame],
    sg_df: Optional[pd.DataFrame],
) -> Dict[str, Any]:
    """Compute 3-way overlap statistics from parser output.

    Uses normalized sequence for fair comparison across engines.
    This is computed at (sequence, charge) level for reporting purposes.
    """

    def get_precursor_set(df: Optional[pd.DataFrame], mod_col: str, charge_col: str) -> Set[Tuple[str, int]]:
        if df is None or df.empty:
            return set()
        # Get normalized sequence for comparison
        sequences = df[mod_col].apply(
            lambda x: normalize_sequence_il(x) if pd.notna(x) else ""
        )
        charges = df[charge_col].astype(int)
        return set(zip(sequences, charges))

    fp_set = get_precursor_set(fp_df, "fragpipe_modified", "fragpipe_charge") if fp_df is not None else set()
    dn_set = get_precursor_set(dn_df, "diann_modified", "diann_charge") if dn_df is not None else set()
    sg_set = get_precursor_set(sg_df, "sage_modified", "sage_charge") if sg_df is not None else set()

    # 3-way Venn regions
    all_three = fp_set & dn_set & sg_set
    fp_dn_only = (fp_set & dn_set) - sg_set
    fp_sg_only = (fp_set & sg_set) - dn_set
    dn_sg_only = (dn_set & sg_set) - fp_set
    fp_only = fp_set - dn_set - sg_set
    dn_only = dn_set - fp_set - sg_set
    sg_only = sg_set - fp_set - dn_set

    union = fp_set | dn_set | sg_set
    at_least_two = (fp_set & dn_set) | (fp_set & sg_set) | (dn_set & sg_set)

    return {
        "n_fragpipe": len(fp_set),
        "n_diann": len(dn_set),
        "n_sage": len(sg_set),
        "n_all_three": len(all_three),
        "n_fp_dn_only": len(fp_dn_only),
        "n_fp_sg_only": len(fp_sg_only),
        "n_dn_sg_only": len(dn_sg_only),
        "n_fragpipe_only": len(fp_only),
        "n_diann_only": len(dn_only),
        "n_sage_only": len(sg_only),
        "n_union": len(union),
        "n_at_least_two": len(at_least_two),
        "three_way_rate": len(all_three) / len(union) if union else 0.0,
        "at_least_two_rate": len(at_least_two) / len(union) if union else 0.0,
    }


def _compute_overlap_stats(
    fp_df: Optional[pd.DataFrame],
    dn_df: Optional[pd.DataFrame],
    sg_df: Optional[pd.DataFrame],
) -> Dict[str, Any]:
    """Compute 3-way overlap statistics (legacy version for non-parser input)."""

    def get_precursor_set(df: Optional[pd.DataFrame], charge_col: str) -> Set[Tuple[str, int]]:
        if df is None or df.empty:
            return set()
        return set(zip(
            df["sequence_normalized"],
            df[charge_col].astype(int)
        ))

    fp_set = get_precursor_set(fp_df, "Charge") if fp_df is not None else set()
    dn_set = get_precursor_set(dn_df, "Precursor.Charge") if dn_df is not None else set()
    sg_set = get_precursor_set(sg_df, "charge") if sg_df is not None else set()

    # 3-way Venn regions
    all_three = fp_set & dn_set & sg_set
    fp_dn_only = (fp_set & dn_set) - sg_set
    fp_sg_only = (fp_set & sg_set) - dn_set
    dn_sg_only = (dn_set & sg_set) - fp_set
    fp_only = fp_set - dn_set - sg_set
    dn_only = dn_set - fp_set - sg_set
    sg_only = sg_set - fp_set - dn_set

    union = fp_set | dn_set | sg_set
    at_least_two = (fp_set & dn_set) | (fp_set & sg_set) | (dn_set & sg_set)

    return {
        "n_fragpipe": len(fp_set),
        "n_diann": len(dn_set),
        "n_sage": len(sg_set),
        "n_all_three": len(all_three),
        "n_fp_dn_only": len(fp_dn_only),
        "n_fp_sg_only": len(fp_sg_only),
        "n_dn_sg_only": len(dn_sg_only),
        "n_fragpipe_only": len(fp_only),
        "n_diann_only": len(dn_only),
        "n_sage_only": len(sg_only),
        "n_union": len(union),
        "n_at_least_two": len(at_least_two),
        "three_way_rate": len(all_three) / len(union) if union else 0.0,
        "at_least_two_rate": len(at_least_two) / len(union) if union else 0.0,
    }


def _build_precursor_index_anchored(
    raw_precursors: pd.DataFrame,
    fp_df: Optional[pd.DataFrame],
    dn_df: Optional[pd.DataFrame],
    sg_df: Optional[pd.DataFrame],
    config: Dict[str, Any],
) -> pd.DataFrame:
    """Build unified precursor index anchored to raw precursor_ids.

    Strategy:
    1. Start from raw precursors (have precursor_id, raw_file, coordinates)
    2. Direct join FragPipe by (raw_file, precursor_id) - FragPipe has this!
    3. Match DIA-NN/Sage by sequence+charge, fallback to coordinates
    4. For engines that don't match raw precursors, add as separate rows

    This ensures the output has precursor_id for Step 5 merge with raw features.

    Schema:
    - Anchor columns: precursor_id, raw_file, mz, charge, rt_seconds, mobility
    - FragPipe columns: fragpipe_* (joined by precursor_id)
    - DIA-NN columns: diann_* (matched by sequence or coordinates)
    - Sage columns: sage_* (matched by sequence or coordinates)
    - Match quality: diann_match_tier, sage_match_tier
    - Derived: n_engines, sequence_normalized
    """
    match_config = MatchConfig()

    # Handle case with no raw precursors - fall back to engine-only index
    if raw_precursors.empty:
        print("    No raw precursors - building engine-only index")
        return _build_precursor_index_engines_only(fp_df, dn_df, sg_df)

    # === Step 1: Start from raw precursors ===
    index = raw_precursors.copy()

    # Ensure required columns exist
    if 'precursor_id' not in index.columns:
        raise ValueError("raw_precursors must have precursor_id column")
    if 'raw_file' not in index.columns:
        raise ValueError("raw_precursors must have raw_file column")

    # Normalize raw_file names (remove .d extension if present)
    index['raw_file'] = index['raw_file'].apply(
        lambda x: str(x).replace('.d', '') if pd.notna(x) else ''
    )

    print(f"    Starting with {len(index)} raw precursors")

    # === Step 2: Direct join FragPipe by (raw_file, precursor_id) ===
    if fp_df is not None and not fp_df.empty:
        # FragPipe parser output already has raw_file and precursor_id
        fp_join = fp_df.copy()

        # Normalize raw_file in FragPipe
        fp_join['raw_file'] = fp_join['raw_file'].apply(
            lambda x: str(x).replace('.d', '') if pd.notna(x) else ''
        )

        # Keep only columns we need
        fp_cols = [c for c in fp_join.columns if c.startswith('fragpipe_')]
        fp_join_cols = ['raw_file', 'precursor_id'] + fp_cols
        fp_join = fp_join[[c for c in fp_join_cols if c in fp_join.columns]]

        # Direct join on (raw_file, precursor_id)
        index = index.merge(
            fp_join,
            on=['raw_file', 'precursor_id'],
            how='left'
        )

        n_fp_joined = index['fragpipe_peptide'].notna().sum() if 'fragpipe_peptide' in index.columns else 0
        print(f"    FragPipe direct join: {n_fp_joined} matches")

        # Add normalized sequence from FragPipe where available
        if 'fragpipe_modified' in index.columns:
            index['sequence_normalized'] = index['fragpipe_modified'].apply(
                lambda x: normalize_sequence_il(x) if pd.notna(x) else ''
            )
        else:
            index['sequence_normalized'] = ''

    else:
        index['sequence_normalized'] = ''

    # === Step 3: Match DIA-NN by sequence+charge, then coordinates ===
    if dn_df is not None and not dn_df.empty:
        index, dn_stats = _match_engine_to_index(
            index=index,
            engine_df=dn_df,
            engine_name='diann',
            sequence_col='diann_modified',
            charge_col='diann_charge',
            mz_col='diann_mz',
            rt_col='diann_rt',
            im_col='diann_mobility',
            config=match_config,
        )
        print(f"    DIA-NN matches: {dn_stats['sequence']} sequence, {dn_stats['coordinate']} coordinate")

    # === Step 4: Direct join Sage by precursor_id (like FragPipe) ===
    # Sage parser converts scannr+1 to precursor_id (raw extraction is 1-indexed)
    if sg_df is not None and not sg_df.empty:
        sg_join = sg_df.copy()

        # Normalize raw_file in Sage
        sg_join['raw_file'] = sg_join['raw_file'].apply(
            lambda x: str(x).replace('.d', '') if pd.notna(x) else ''
        )

        # Keep only columns we need
        sg_cols = [c for c in sg_join.columns if c.startswith('sage_')]
        sg_join_cols = ['raw_file', 'precursor_id'] + sg_cols
        sg_join = sg_join[[c for c in sg_join_cols if c in sg_join.columns]]

        # Direct join on (raw_file, precursor_id)
        index = index.merge(
            sg_join,
            on=['raw_file', 'precursor_id'],
            how='left'
        )

        n_sg_joined = index['sage_peptide'].notna().sum() if 'sage_peptide' in index.columns else 0
        print(f"    Sage direct join: {n_sg_joined} matches")

    # === Step 5: Count engines ===
    fp_present = index.get('fragpipe_peptide', pd.Series(dtype=str)).notna()
    dn_present = index.get('diann_peptide', pd.Series(dtype=str)).notna()
    sg_present = index.get('sage_peptide', pd.Series(dtype=str)).notna()

    index['n_engines'] = fp_present.astype(int) + dn_present.astype(int) + sg_present.astype(int)

    print(f"    Engine distribution: 3={int((index['n_engines']==3).sum())}, "
          f"2={int((index['n_engines']==2).sum())}, 1={int((index['n_engines']==1).sum())}, "
          f"0={int((index['n_engines']==0).sum())}")

    return index


def _match_engine_to_index(
    index: pd.DataFrame,
    engine_df: pd.DataFrame,
    engine_name: str,
    sequence_col: str,
    charge_col: str,
    mz_col: str,
    rt_col: str,
    im_col: str,
    config: MatchConfig,
) -> Tuple[pd.DataFrame, Dict[str, int]]:
    """Match engine results to existing index using tiered strategy.

    Strategy:
    1. For rows WITH FragPipe ID: match by sequence+charge
    2. For rows WITHOUT: match by coordinates (m/z + RT + IM)

    Args:
        index: Current precursor index (has precursor_id, raw_file, optionally sequence_normalized)
        engine_df: Engine results to match
        engine_name: Engine name for column prefixes
        sequence_col: Column with modified sequence in engine_df
        charge_col: Column with charge in engine_df
        mz_col: Column with m/z in engine_df
        rt_col: Column with RT in engine_df (seconds)
        im_col: Column with ion mobility in engine_df
        config: Match configuration

    Returns:
        (updated_index, stats_dict)
    """
    stats = {'sequence': 0, 'coordinate': 0}

    # Initialize engine columns
    engine_cols = [c for c in engine_df.columns if c.startswith(f'{engine_name}_')]
    for col in engine_cols:
        if col not in index.columns:
            index[col] = None
    index[f'{engine_name}_match_tier'] = None

    # Normalize raw_file in engine_df
    engine_df = engine_df.copy()
    engine_df['raw_file'] = engine_df['raw_file'].apply(
        lambda x: str(x).replace('.d', '') if pd.notna(x) else ''
    )

    # Create normalized sequence column for matching
    engine_df['_norm_seq'] = engine_df[sequence_col].apply(
        lambda x: normalize_sequence_il(x) if pd.notna(x) else ''
    )

    matched_engine_indices: Set[int] = set()

    # Process per raw file
    for raw_file in index['raw_file'].unique():
        file_mask = index['raw_file'] == raw_file
        engine_file_df = engine_df[engine_df['raw_file'] == raw_file]

        if engine_file_df.empty:
            continue

        # === Pass 1: Sequence matching for rows with FragPipe ID ===
        has_fp = file_mask & index['sequence_normalized'].notna() & (index['sequence_normalized'] != '')

        # Build sequence index on engine data
        seq_index: Dict[Tuple[str, int], List[int]] = {}
        for idx, row in engine_file_df.iterrows():
            norm_seq = row['_norm_seq']
            charge = row.get(charge_col)
            if not norm_seq or pd.isna(charge):
                continue
            key = (norm_seq, int(charge))
            if key not in seq_index:
                seq_index[key] = []
            seq_index[key].append(idx)

        for idx in index[has_fp].index:
            row = index.loc[idx]
            norm_seq = row['sequence_normalized']
            charge = row.get('charge')

            if pd.isna(charge) or not norm_seq:
                continue

            key = (norm_seq, int(charge))
            if key in seq_index:
                for engine_idx in seq_index[key]:
                    if engine_idx not in matched_engine_indices:
                        # Copy engine columns
                        engine_row = engine_df.loc[engine_idx]
                        for col in engine_cols:
                            index.loc[idx, col] = engine_row.get(col)
                        index.loc[idx, f'{engine_name}_match_tier'] = 'SEQUENCE_IL_NORM'
                        matched_engine_indices.add(engine_idx)
                        stats['sequence'] += 1
                        break

        # === Pass 2: Coordinate matching for unmatched rows ===
        unmatched_engine = engine_file_df[~engine_file_df.index.isin(matched_engine_indices)]
        if unmatched_engine.empty:
            continue

        # Build coordinate index (m/z bins)
        bin_size = 0.01
        coord_index: Dict[Tuple[int, int], List[Tuple[int, float]]] = {}
        for idx, row in unmatched_engine.iterrows():
            mz = row.get(mz_col)
            charge = row.get(charge_col)
            if pd.isna(mz) or pd.isna(charge):
                continue
            mz = float(mz)
            mz_bin = int(mz / bin_size)
            charge = int(charge)
            for offset in [-1, 0, 1]:
                key = (mz_bin + offset, charge)
                if key not in coord_index:
                    coord_index[key] = []
                coord_index[key].append((idx, mz))

        # Match unassigned index rows
        no_engine = file_mask & index[f'{engine_name}_match_tier'].isna()

        for idx in index[no_engine].index:
            row = index.loc[idx]
            source_mz = row.get('mz')
            source_charge = row.get('charge')

            if pd.isna(source_mz) or pd.isna(source_charge):
                continue

            source_mz = float(source_mz)
            source_charge = int(source_charge)
            mz_bin = int(source_mz / bin_size)
            key = (mz_bin, source_charge)

            if key not in coord_index:
                continue

            mz_tol = source_mz * config.mz_tol_ppm / 1e6
            source_rt = row.get('rt_seconds')
            source_im = row.get('mobility')

            best_match = None
            best_score = float('inf')

            for engine_idx, target_mz in coord_index[key]:
                if engine_idx in matched_engine_indices:
                    continue

                mz_diff = abs(source_mz - target_mz)
                if mz_diff > mz_tol:
                    continue

                engine_row = engine_df.loc[engine_idx]
                score = mz_diff / mz_tol

                # Check RT if available
                target_rt = engine_row.get(rt_col)
                if pd.notna(source_rt) and pd.notna(target_rt):
                    rt_diff = abs(float(source_rt) - float(target_rt))
                    if rt_diff > config.rt_tol_sec:
                        continue
                    score += rt_diff / config.rt_tol_sec

                # Check IM if available
                target_im = engine_row.get(im_col)
                if pd.notna(source_im) and pd.notna(target_im):
                    im_diff = abs(float(source_im) - float(target_im))
                    if im_diff > config.im_tol:
                        continue
                    score += im_diff / config.im_tol

                if score < best_score:
                    best_score = score
                    best_match = engine_idx

            if best_match is not None:
                engine_row = engine_df.loc[best_match]
                for col in engine_cols:
                    index.loc[idx, col] = engine_row.get(col)
                index.loc[idx, f'{engine_name}_match_tier'] = 'COORDINATE_FULL'

                # Update sequence_normalized if not set
                if not index.loc[idx, 'sequence_normalized']:
                    index.loc[idx, 'sequence_normalized'] = engine_row['_norm_seq']

                matched_engine_indices.add(best_match)
                stats['coordinate'] += 1

    return index, stats


def _build_precursor_index_engines_only(
    fp_df: Optional[pd.DataFrame],
    dn_df: Optional[pd.DataFrame],
    sg_df: Optional[pd.DataFrame],
) -> pd.DataFrame:
    """Build precursor index from engine results only (no raw precursors).

    Used as fallback when raw_features.parquet is not available.
    Groups by sequence+charge and merges across engines.
    """
    dfs_to_merge = []

    if fp_df is not None and not fp_df.empty:
        # Add normalized sequence
        fp_df = fp_df.copy()
        fp_df['sequence_normalized'] = fp_df['fragpipe_modified'].apply(
            lambda x: normalize_sequence_il(x) if pd.notna(x) else ''
        )
        fp_df['charge'] = fp_df['fragpipe_charge']

        # Keep best per (raw_file, precursor_id) if available, else (sequence, charge)
        if 'precursor_id' in fp_df.columns:
            fp_summary = fp_df.groupby(['raw_file', 'precursor_id']).first().reset_index()
        else:
            fp_summary = fp_df.groupby(['sequence_normalized', 'charge']).first().reset_index()

        dfs_to_merge.append(fp_summary)

    if dn_df is not None and not dn_df.empty:
        dn_df = dn_df.copy()
        dn_df['sequence_normalized'] = dn_df['diann_modified'].apply(
            lambda x: normalize_sequence_il(x) if pd.notna(x) else ''
        )
        dn_df['charge'] = dn_df['diann_charge']
        dn_summary = dn_df.groupby(['sequence_normalized', 'charge']).first().reset_index()
        dfs_to_merge.append(dn_summary)

    if sg_df is not None and not sg_df.empty:
        sg_df = sg_df.copy()
        sg_df['sequence_normalized'] = sg_df['sage_modified'].apply(
            lambda x: normalize_sequence_il(x) if pd.notna(x) else ''
        )
        sg_df['charge'] = sg_df['sage_charge']
        sg_summary = sg_df.groupby(['sequence_normalized', 'charge']).first().reset_index()
        dfs_to_merge.append(sg_summary)

    if not dfs_to_merge:
        return pd.DataFrame(columns=['sequence_normalized', 'charge', 'n_engines'])

    # Merge all on sequence_normalized + charge
    merged = dfs_to_merge[0]
    for df in dfs_to_merge[1:]:
        merged = merged.merge(df, on=['sequence_normalized', 'charge'], how='outer', suffixes=('', '_dup'))
        # Remove duplicate columns
        merged = merged[[c for c in merged.columns if not c.endswith('_dup')]]

    # Count engines
    fp_present = merged.get('fragpipe_peptide', pd.Series(dtype=str)).notna()
    dn_present = merged.get('diann_peptide', pd.Series(dtype=str)).notna()
    sg_present = merged.get('sage_peptide', pd.Series(dtype=str)).notna()

    merged['n_engines'] = fp_present.astype(int) + dn_present.astype(int) + sg_present.astype(int)

    return merged


def _build_precursor_index(
    fp_df: Optional[pd.DataFrame],
    dn_df: Optional[pd.DataFrame],
    sg_df: Optional[pd.DataFrame],
    processed_dir: Path,
) -> pd.DataFrame:
    """Build unified precursor index (legacy version).

    DEPRECATED: Use _build_precursor_index_anchored instead.
    Kept for backward compatibility.
    """
    # This is the old implementation that groups by (sequence, charge)
    # and loses precursor_id linkage
    dfs_to_merge = []

    if fp_df is not None and not fp_df.empty:
        fp_summary = fp_df.groupby(["sequence_normalized", "Charge"]).first().reset_index()
        fp_summary = fp_summary.rename(columns={"Charge": "charge"})

        fp_summary["fragpipe_modified"] = fp_summary.get("modified_unimod", "")
        fp_summary["fragpipe_peptide"] = fp_summary.get("Peptide Sequence", fp_summary.get("Peptide", ""))
        fp_summary["fragpipe_protein"] = fp_summary.get("Protein", "")
        fp_summary["fragpipe_probability"] = fp_summary.get("Probability")
        fp_summary["fragpipe_pep"] = 1.0 - fp_summary.get("Probability", 0.0)
        fp_summary["fragpipe_hyperscore"] = fp_summary.get("Hyperscore")
        fp_summary["fragpipe_qvalue"] = fp_summary.get("Qvalue")
        fp_summary["fragpipe_rt"] = fp_summary.get("Retention")
        fp_summary["fragpipe_mz"] = fp_summary.get("Calibrated Observed M/Z", fp_summary.get("Observed M/Z"))
        fp_summary["fragpipe_mobility"] = fp_summary.get("Ion Mobility")

        fp_cols = [
            "sequence_normalized", "charge",
            "fragpipe_modified", "fragpipe_peptide", "fragpipe_protein",
            "fragpipe_probability", "fragpipe_pep", "fragpipe_hyperscore", "fragpipe_qvalue",
            "fragpipe_rt", "fragpipe_mz", "fragpipe_mobility",
        ]
        fp_summary = fp_summary[[c for c in fp_cols if c in fp_summary.columns]]
        dfs_to_merge.append(("fragpipe", fp_summary))

    if dn_df is not None and not dn_df.empty:
        dn_summary = dn_df.groupby(["sequence_normalized", "Precursor.Charge"]).first().reset_index()
        dn_summary = dn_summary.rename(columns={"Precursor.Charge": "charge"})

        dn_summary["diann_modified"] = dn_summary.get("modified_unimod", "")
        dn_summary["diann_peptide"] = dn_summary.get("Stripped.Sequence", "")
        dn_summary["diann_protein"] = dn_summary.get("Protein.Ids", "")
        dn_summary["diann_qvalue"] = dn_summary.get("Q.Value")
        dn_summary["diann_pep"] = dn_summary.get("PEP")
        dn_summary["diann_global_qvalue"] = dn_summary.get("Global.Q.Value")
        dn_summary["diann_pg_qvalue"] = dn_summary.get("PG.Q.Value")

        rt_minutes = dn_summary.get("RT")
        if rt_minutes is not None:
            dn_summary["diann_rt"] = rt_minutes * 60.0
        else:
            dn_summary["diann_rt"] = None
        dn_summary["diann_mz"] = dn_summary.get("Precursor.Mz")
        dn_summary["diann_mobility"] = dn_summary.get("IM")
        dn_summary["diann_ccs"] = dn_summary.get("CCS")

        dn_cols = [
            "sequence_normalized", "charge",
            "diann_modified", "diann_peptide", "diann_protein",
            "diann_qvalue", "diann_pep", "diann_global_qvalue", "diann_pg_qvalue",
            "diann_rt", "diann_mz", "diann_mobility", "diann_ccs",
        ]
        dn_summary = dn_summary[[c for c in dn_cols if c in dn_summary.columns]]
        dfs_to_merge.append(("diann", dn_summary))

    if sg_df is not None and not sg_df.empty:
        sg_summary = sg_df.groupby(["sequence_normalized", "charge"]).first().reset_index()

        sg_summary["sage_modified"] = sg_summary.get("modified_unimod", "")
        import re
        def strip_mods(seq):
            if pd.isna(seq) or not seq:
                return ""
            return re.sub(r'\[[^\]]+\]', '', str(seq))
        if "peptide" in sg_summary.columns:
            sg_summary["sage_peptide"] = sg_summary["peptide"].apply(strip_mods)
        else:
            sg_summary["sage_peptide"] = ""
        sg_summary["sage_protein"] = sg_summary.get("proteins", "")
        sg_summary["sage_qvalue"] = sg_summary.get("spectrum_q")
        sg_summary["sage_pep"] = sg_summary.get("posterior_error").apply(
            lambda x: np.exp(x) if pd.notna(x) else np.nan
        ) if "posterior_error" in sg_summary.columns else np.nan
        sg_summary["sage_hyperscore"] = sg_summary.get("hyperscore")
        sg_summary["sage_peptide_qvalue"] = sg_summary.get("peptide_q")
        sg_summary["sage_protein_qvalue"] = sg_summary.get("protein_q")
        sg_summary["sage_rt"] = sg_summary.get("rt")
        PROTON_MASS = 1.007276
        if "expmass" in sg_summary.columns:
            sg_summary["sage_mz"] = (sg_summary["expmass"] + sg_summary["charge"] * PROTON_MASS) / sg_summary["charge"]
        else:
            sg_summary["sage_mz"] = None
        sg_summary["sage_mobility"] = sg_summary.get("ion_mobility")

        sg_cols = [
            "sequence_normalized", "charge",
            "sage_modified", "sage_peptide", "sage_protein",
            "sage_qvalue", "sage_pep", "sage_hyperscore", "sage_peptide_qvalue", "sage_protein_qvalue",
            "sage_rt", "sage_mz", "sage_mobility",
        ]
        sg_summary = sg_summary[[c for c in sg_cols if c in sg_summary.columns]]
        dfs_to_merge.append(("sage", sg_summary))

    if not dfs_to_merge:
        return pd.DataFrame(columns=["sequence_normalized", "charge", "n_engines"])

    _, merged = dfs_to_merge[0]
    for name, df in dfs_to_merge[1:]:
        merged = merged.merge(df, on=["sequence_normalized", "charge"], how="outer")

    fp_present = merged.get("fragpipe_peptide", pd.Series(dtype=str)).notna() & (merged.get("fragpipe_peptide", "") != "")
    dn_present = merged.get("diann_peptide", pd.Series(dtype=str)).notna() & (merged.get("diann_peptide", "") != "")
    sg_present = merged.get("sage_peptide", pd.Series(dtype=str)).notna() & (merged.get("sage_peptide", "") != "")

    merged["n_engines"] = fp_present.astype(int) + dn_present.astype(int) + sg_present.astype(int)

    return merged


def _stratify_precursors(
    index_df: pd.DataFrame,
    stratified_dir: Path,
) -> Dict[str, int]:
    """Stratify precursors by engine agreement."""
    stratified_dir.mkdir(parents=True, exist_ok=True)

    # Define conditions - handle missing columns for engines with 0 results
    has_fp = index_df["fragpipe_peptide"].notna() if "fragpipe_peptide" in index_df.columns else pd.Series(False, index=index_df.index)
    has_dn = index_df["diann_peptide"].notna() if "diann_peptide" in index_df.columns else pd.Series(False, index=index_df.index)
    has_sg = index_df["sage_peptide"].notna() if "sage_peptide" in index_df.columns else pd.Series(False, index=index_df.index)

    strata = {
        "all_three": has_fp & has_dn & has_sg,
        "two_plus": (has_fp & has_dn) | (has_fp & has_sg) | (has_dn & has_sg),
        "fragpipe_only": has_fp & ~has_dn & ~has_sg,
        "diann_only": ~has_fp & has_dn & ~has_sg,
        "sage_only": ~has_fp & ~has_dn & has_sg,
    }

    counts = {}
    for name, mask in strata.items():
        subset = index_df[mask]
        if not subset.empty:
            output_path = stratified_dir / f"{name}.parquet"
            subset.to_parquet(output_path, index=False)
        counts[name] = int(mask.sum())

    return counts


def _compute_match_tiers(index_df: pd.DataFrame) -> Dict[str, int]:
    """Compute match tier distribution."""
    # For now, return engine agreement tiers
    return {
        "all_three": int((index_df["n_engines"] == 3).sum()),
        "two_engines": int((index_df["n_engines"] == 2).sum()),
        "one_engine": int((index_df["n_engines"] == 1).sum()),
    }


def _generate_html_report(
    accession: str,
    stats: Dict[str, Any],
    output_path: Path,
) -> None:
    """Generate HTML overlap report."""
    from scripts.analyze_overlap import HTML_TEMPLATE_3WAY

    # Calculate percentages
    total = stats["n_union"]
    if total == 0:
        total = 1  # Avoid division by zero

    html = HTML_TEMPLATE_3WAY.format(
        accession=accession,
        timestamp="Generated by San José Runner",
        match_type="sequence+charge",
        pep_filter="None (all results)",
        n_fragpipe=stats["n_fragpipe"],
        n_diann=stats["n_diann"],
        n_sage=stats["n_sage"],
        n_all_three=stats["n_all_three"],
        n_at_least_two=stats["n_at_least_two"],
        n_union=stats["n_union"],
        three_way_pct=stats["three_way_rate"],
        two_plus_pct=stats["at_least_two_rate"],
        n_fp_dn=stats["n_all_three"] + stats["n_fp_dn_only"],
        n_fp_sg=stats["n_all_three"] + stats["n_fp_sg_only"],
        n_dn_sg=stats["n_all_three"] + stats["n_dn_sg_only"],
        n_fp_dn_only=stats["n_fp_dn_only"],
        n_fp_sg_only=stats["n_fp_sg_only"],
        n_dn_sg_only=stats["n_dn_sg_only"],
        n_fragpipe_only=stats["n_fragpipe_only"],
        n_diann_only=stats["n_diann_only"],
        n_sage_only=stats["n_sage_only"],
        all_three_pct=stats["n_all_three"] / total * 100,
        fp_dn_only_pct=stats["n_fp_dn_only"] / total * 100,
        fp_sg_only_pct=stats["n_fp_sg_only"] / total * 100,
        dn_sg_only_pct=stats["n_dn_sg_only"] / total * 100,
        fp_only_pct=stats["n_fragpipe_only"] / total * 100,
        dn_only_pct=stats["n_diann_only"] / total * 100,
        sg_only_pct=stats["n_sage_only"] / total * 100,
        fp_unique_pct=stats["n_fragpipe_only"] / max(stats["n_fragpipe"], 1) * 100,
        dn_unique_pct=stats["n_diann_only"] / max(stats["n_diann"], 1) * 100,
        sg_unique_pct=stats["n_sage_only"] / max(stats["n_sage"], 1) * 100,
        fp_validation_rate=1 - stats["n_fragpipe_only"] / max(stats["n_fragpipe"], 1),
        dn_validation_rate=1 - stats["n_diann_only"] / max(stats["n_diann"], 1),
        sg_validation_rate=1 - stats["n_sage_only"] / max(stats["n_sage"], 1),
        charge_rows="<tr><td colspan='4'>N/A</td></tr>",
        fp_len_min="N/A",
        fp_len_median="N/A",
        fp_len_max="N/A",
        dn_len_min="N/A",
        dn_len_median="N/A",
        dn_len_max="N/A",
        sg_len_min="N/A",
        sg_len_median="N/A",
        sg_len_max="N/A",
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w") as f:
        f.write(html)


if __name__ == "__main__":
    import argparse
    import yaml

    parser = argparse.ArgumentParser(description="Step 3: Stratify search results")
    parser.add_argument("--accession", "-a", required=True, help="PRIDE accession")
    parser.add_argument("--config", "-c", default="config/config.yaml", help="Config file")
    parser.add_argument("--output-dir", "-o", default="data", help="Output base directory")
    parser.add_argument("--no-html", action="store_true", help="Skip HTML report")

    args = parser.parse_args()

    # Load config
    with open(args.config) as f:
        config = yaml.safe_load(f)

    # Run step
    summary = run_step3_stratify(
        accession=args.accession,
        config=config,
        output_base_dir=Path(args.output_dir),
        generate_html=not args.no_html,
    )

    print(f"\nStep 3 completed: {summary.status}")
    print(f"  Total precursors: {summary.data['n_total_precursors']}")
    print(f"  Stratified counts: {summary.data['stratified_counts']}")
