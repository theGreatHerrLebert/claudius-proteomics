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

Per-group mode (when sample_groups.yaml exists from Step 1):
- Iterates over sample groups, produces one index + consensus per group
- Writes cross-group QC manifest at accession level

Outputs (per-group mode):
- data/processed/{accession}/{group_id}/precursor_index.parquet
- data/processed/{accession}/{group_id}/consensus/
    - overlap_stats.json
    - overlap_report.html
    - stratified/
        - all_three.parquet
        - two_plus.parquet
        - fragpipe_only.parquet
        - diann_only.parquet
        - sage_only.parquet
- data/processed/{accession}/step3_qc_manifest.json
- step3_summary.json

Outputs (legacy single-group mode):
- data/processed/{accession}/precursor_index.parquet
- data/processed/{accession}/consensus/
- step3_summary.json
"""

import json
import sys
from datetime import datetime
from pathlib import Path
from typing import Dict, Any, Optional, Set, Tuple, List

import pandas as pd

# Add parent to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from lib.sequence_utils import normalize_sequence_il
from lib.precursor_matching import MatchConfig
from runner.summary import StepSummary, write_step_summary
from scripts.engine_parsers.fragpipe_parser import FragPipeParser
from scripts.engine_parsers.diann_parser import DiannParser
from scripts.engine_parsers.sage_parser import SageParser
from scripts.sample_group_resolver import SampleGroupManifest


def run_step3_stratify(
    accession: str,
    config: Dict[str, Any],
    output_base_dir: Path,
    generate_html: bool = True,
) -> StepSummary:
    """
    Execute Step 3: Stratify and merge search engine results.

    If sample_groups.yaml exists (written by step 1), runs stratification per
    sample group with group-specific processed dirs. Otherwise falls back to
    single-group legacy behavior.

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
    metadata_dir = output_base_dir / "metadata" / accession

    try:
        # Try to load sample groups from step 1
        sg_path = metadata_dir / "sample_groups.yaml"
        if sg_path.exists():
            manifest = SampleGroupManifest.from_yaml(sg_path)
            print(f"  Loaded {len(manifest.groups)} sample groups from {sg_path}")
        else:
            manifest = None
            print(f"  No sample_groups.yaml found — using single-group mode")

        if manifest and len(manifest.groups) > 0:
            results = _run_per_group(
                manifest=manifest,
                accession=accession,
                config=config,
                processed_dir=processed_dir,
                extracted_dir=extracted_dir,
                generate_html=generate_html,
            )
        else:
            results = _run_single_group(
                accession=accession,
                config=config,
                processed_dir=processed_dir,
                extracted_dir=extracted_dir,
                generate_html=generate_html,
            )

        summary.data = results
        summary.outputs = [str(processed_dir)]
        summary.complete(success=True)

    except Exception as e:
        summary.complete(success=False, error_message=str(e))
        raise

    # Write summary file
    write_step_summary(summary, processed_dir)

    return summary


def _run_per_group(
    manifest: SampleGroupManifest,
    accession: str,
    config: Dict[str, Any],
    processed_dir: Path,
    extracted_dir: Path,
    generate_html: bool,
) -> Dict[str, Any]:
    """Run stratification for each sample group.

    Args:
        manifest: Sample group manifest from step 1
        accession: PRIDE accession
        config: Pipeline configuration dict
        processed_dir: e.g. data/processed/PXD019086
        extracted_dir: e.g. data/extracted/PXD019086
        generate_html: Whether to generate HTML reports

    Returns:
        Aggregated results dict with per-group breakdown
    """
    group_results = {}
    n_skipped = 0

    for group in manifest.groups:
        print(f"\n  === Sample group: {group.group_id} ===")
        print(f"      Organism: {group.organism_name} ({group.organism_key})")
        print(f"      Enzyme:   {group.enzyme}")
        print(f"      Runs:     {group.n_runs}")

        if group.n_runs == 0:
            print(f"      Skipping (no runs)")
            group_results[group.group_id] = {"status": "skipped", "reason": "no runs"}
            n_skipped += 1
            continue

        group_processed_dir = processed_dir / group.group_id
        group_extracted_dir = extracted_dir / group.group_id

        try:
            group_stats = _stratify_single_group(
                accession=accession,
                config=config,
                processed_dir=group_processed_dir,
                extracted_dir=group_extracted_dir,
                generate_html=generate_html,
                report_label=f"{accession}/{group.group_id}",
            )
            group_results[group.group_id] = {
                "status": "success",
                "organism": group.organism_key,
                "enzyme": group.enzyme,
                **group_stats,
            }
        except Exception as e:
            print(f"      ERROR: {e}")
            group_results[group.group_id] = {
                "status": "error",
                "organism": group.organism_key,
                "enzyme": group.enzyme,
                "error": str(e),
            }

    # Write cross-group QC manifest (always, for diagnostics)
    _write_qc_manifest(accession, group_results, processed_dir)

    # Fail if any group had errors
    failed_groups = [
        gid for gid, r in group_results.items() if r.get("status") == "error"
    ]
    if failed_groups:
        errors = "; ".join(
            f"{gid}: {group_results[gid].get('error', 'unknown')}"
            for gid in failed_groups
        )
        raise RuntimeError(
            f"Stratification failed for {len(failed_groups)} group(s): {errors}"
        )

    # Compute totals across groups
    total_precursors = sum(
        r.get("n_total_precursors", 0)
        for r in group_results.values()
        if r.get("status") == "success"
    )

    return {
        "mode": "per_group",
        "n_groups": len(manifest.groups),
        "n_groups_processed": len(manifest.groups) - n_skipped,
        "n_groups_skipped": n_skipped,
        "n_total_precursors": total_precursors,
        "groups": group_results,
    }


def _run_single_group(
    accession: str,
    config: Dict[str, Any],
    processed_dir: Path,
    extracted_dir: Path,
    generate_html: bool,
) -> Dict[str, Any]:
    """Legacy single-group mode: flat processed_dir, no sample groups."""
    stats = _stratify_single_group(
        accession=accession,
        config=config,
        processed_dir=processed_dir,
        extracted_dir=extracted_dir,
        generate_html=generate_html,
        report_label=accession,
    )
    return {
        "mode": "single_group",
        **stats,
    }


def _stratify_single_group(
    accession: str,
    config: Dict[str, Any],
    processed_dir: Path,
    extracted_dir: Path,
    generate_html: bool,
    report_label: str,
) -> Dict[str, Any]:
    """Stratify a single group's search engine results.

    This is the core stratification logic, operating on a single processed_dir
    (either the flat accession dir in legacy mode, or a group subdir).

    Args:
        accession: PRIDE accession (for parser fallback)
        config: Pipeline configuration dict
        processed_dir: Directory with engine results (canonical parquets)
        extracted_dir: Directory with raw features (may not exist)
        generate_html: Whether to generate HTML report
        report_label: Label for HTML report (accession or accession/group_id)

    Returns:
        Dict with overlap_stats, stratified_counts, match_tiers, n_total_precursors
    """
    consensus_dir = processed_dir / "consensus"
    stratified_dir = consensus_dir / "stratified"

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

    # Load search engine results (canonical parquet preferred, parser fallback)
    print("  Loading search engine results...")
    fp_df = _load_engine(processed_dir, "fragpipe", accession=accession)
    dn_df = _load_engine(processed_dir, "diann", accession=accession)
    sg_df = _load_engine(processed_dir, "sage", accession=accession)

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

    # Compute QC stats from engine DataFrames
    print("  Computing QC statistics...")
    qc_stats = _compute_qc_stats(fp_df, dn_df, sg_df)

    # Generate QC plots and HTML report
    if generate_html:
        print("  Generating QC plots...")
        qc_plots = _generate_qc_plots(precursor_index)
        print(f"    Generated {len(qc_plots)} plots")

        html_path = consensus_dir / "overlap_report.html"
        _generate_html_report(report_label, overlap_stats, html_path, qc_stats=qc_stats, plots=qc_plots)
        print(f"  Generated HTML report: {html_path}")

    return {
        "overlap_stats": overlap_stats,
        "stratified_counts": stratified_counts,
        "match_tiers": match_tiers,
        "n_total_precursors": len(precursor_index),
    }


def _write_qc_manifest(
    accession: str,
    group_results: Dict[str, Dict[str, Any]],
    processed_dir: Path,
) -> None:
    """Write cross-group QC manifest at the accession level.

    Args:
        accession: PRIDE accession
        group_results: Per-group results from _run_per_group
        processed_dir: e.g. data/processed/PXD019086
    """
    groups_qc = {}
    total_precursors = 0
    n_processed = 0
    n_skipped = 0

    for group_id, result in group_results.items():
        if result.get("status") != "success":
            n_skipped += 1
            continue

        n_processed += 1
        n_prec = result.get("n_total_precursors", 0)
        total_precursors += n_prec

        overlap = result.get("overlap_stats", {})
        match_tiers = result.get("match_tiers", {})

        groups_qc[group_id] = {
            "organism": result.get("organism", ""),
            "enzyme": result.get("enzyme", ""),
            "n_precursors": n_prec,
            "n_engines": {
                "3": match_tiers.get("all_three", 0),
                "2": match_tiers.get("two_engines", 0),
                "1": match_tiers.get("one_engine", 0),
            },
            "overlap": {
                "three_way_rate": overlap.get("three_way_rate", 0.0),
                "at_least_two_rate": overlap.get("at_least_two_rate", 0.0),
                "n_union": overlap.get("n_union", 0),
            },
        }

    qc_manifest = {
        "accession": accession,
        "generated_at": datetime.now().isoformat(),
        "n_groups": len(group_results),
        "groups": groups_qc,
        "totals": {
            "n_precursors": total_precursors,
            "n_groups_processed": n_processed,
            "n_groups_skipped": n_skipped,
        },
    }

    processed_dir.mkdir(parents=True, exist_ok=True)
    qc_path = processed_dir / "step3_qc_manifest.json"
    with open(qc_path, "w") as f:
        json.dump(qc_manifest, f, indent=2)
    print(f"\n  Wrote QC manifest: {qc_path}")


def _load_engine(
    processed_dir: Path,
    engine_name: str,
    accession: Optional[str] = None,
) -> Optional[pd.DataFrame]:
    """Load engine results, preferring canonical parquet from step2.

    Lookup order:
    1. {engine}_canonical.parquet (produced by runner/engines/ jobs)
    2. Fallback: re-parse native output using the engine parser
       (only in legacy single-group mode where processed_dir.name == accession)

    Args:
        processed_dir: e.g. data/processed/PXD019086 or data/processed/PXD019086/human_trypsin
        engine_name: "fragpipe", "diann", or "sage"
        accession: PRIDE accession for parser fallback. If provided and
            processed_dir.name != accession, parser fallback is skipped
            (per-group mode always has canonical parquets from step 2).

    Returns:
        DataFrame with {engine}_* prefixed columns, or None/empty DataFrame
    """
    # Prefer canonical parquet written by the engine job
    canonical = processed_dir / f"{engine_name}_canonical.parquet"
    if canonical.exists():
        return pd.read_parquet(canonical)

    # In per-group mode, skip parser fallback (canonical parquet is guaranteed)
    if accession and processed_dir.name != accession:
        return None

    # Fallback: re-parse native output (legacy single-group mode only)
    parsers = {"fragpipe": FragPipeParser, "diann": DiannParser, "sage": SageParser}
    parser_cls = parsers.get(engine_name)
    if parser_cls is None:
        return None
    base_dir = processed_dir.parent
    fallback_accession = accession if accession else processed_dir.name
    return parser_cls().parse(base_dir, fallback_accession)


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
    Groups by (raw_file, sequence, charge) and merges across engines to
    preserve per-run identity and avoid inflating consensus counts.
    """
    dfs_to_merge = []

    if fp_df is not None and not fp_df.empty:
        # Add normalized sequence
        fp_df = fp_df.copy()
        fp_df['sequence_normalized'] = fp_df['fragpipe_modified'].apply(
            lambda x: normalize_sequence_il(x) if pd.notna(x) else ''
        )
        fp_df['charge'] = fp_df['fragpipe_charge']

        # Keep best per (raw_file, precursor_id) if available, else (raw_file, sequence, charge)
        if 'precursor_id' in fp_df.columns:
            fp_summary = fp_df.groupby(['raw_file', 'precursor_id']).first().reset_index()
        else:
            fp_summary = fp_df.groupby(['raw_file', 'sequence_normalized', 'charge']).first().reset_index()

        dfs_to_merge.append(fp_summary)

    if dn_df is not None and not dn_df.empty:
        dn_df = dn_df.copy()
        dn_df['sequence_normalized'] = dn_df['diann_modified'].apply(
            lambda x: normalize_sequence_il(x) if pd.notna(x) else ''
        )
        dn_df['charge'] = dn_df['diann_charge']
        dn_summary = dn_df.groupby(['raw_file', 'sequence_normalized', 'charge']).first().reset_index()
        dfs_to_merge.append(dn_summary)

    if sg_df is not None and not sg_df.empty:
        sg_df = sg_df.copy()
        sg_df['sequence_normalized'] = sg_df['sage_modified'].apply(
            lambda x: normalize_sequence_il(x) if pd.notna(x) else ''
        )
        sg_df['charge'] = sg_df['sage_charge']
        sg_summary = sg_df.groupby(['raw_file', 'sequence_normalized', 'charge']).first().reset_index()
        dfs_to_merge.append(sg_summary)

    if not dfs_to_merge:
        return pd.DataFrame(columns=['raw_file', 'sequence_normalized', 'charge', 'n_engines'])

    # Merge all on (raw_file, sequence_normalized, charge) to preserve per-run identity
    merge_keys = ['raw_file', 'sequence_normalized', 'charge']
    merged = dfs_to_merge[0]
    for df in dfs_to_merge[1:]:
        merged = merged.merge(df, on=merge_keys, how='outer', suffixes=('', '_dup'))
        # Remove duplicate columns
        merged = merged[[c for c in merged.columns if not c.endswith('_dup')]]

    # Count engines
    fp_present = merged.get('fragpipe_peptide', pd.Series(dtype=str)).notna()
    dn_present = merged.get('diann_peptide', pd.Series(dtype=str)).notna()
    sg_present = merged.get('sage_peptide', pd.Series(dtype=str)).notna()

    merged['n_engines'] = fp_present.astype(int) + dn_present.astype(int) + sg_present.astype(int)

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


def _compute_qc_stats(
    fp_df: Optional[pd.DataFrame],
    dn_df: Optional[pd.DataFrame],
    sg_df: Optional[pd.DataFrame],
) -> Dict[str, Any]:
    """Compute exhaustive QC statistics from the three engine DataFrames.

    Returns a dict with per-engine charge distributions, sequence length stats,
    ion mobility stats, PEP score distributions, per-raw-file PSM counts,
    and top modifications.
    """
    import numpy as np

    engines = {
        "fragpipe": fp_df,
        "diann": dn_df,
        "sage": sg_df,
    }

    qc: Dict[str, Any] = {}

    # --- Charge distribution per engine ---
    charge_dist: Dict[str, Dict[int, int]] = {}
    for name, df in engines.items():
        col = f"{name}_charge"
        if df is not None and not df.empty and col in df.columns:
            counts = df[col].dropna().astype(int).value_counts().sort_index()
            charge_dist[name] = {int(k): int(v) for k, v in counts.items()}
        else:
            charge_dist[name] = {}
    qc["charge_dist"] = charge_dist

    # --- Sequence length stats per engine ---
    seq_length: Dict[str, Dict[str, Any]] = {}
    for name, df in engines.items():
        col = f"{name}_peptide"
        if df is not None and not df.empty and col in df.columns:
            lengths = df[col].dropna().str.len()
            if len(lengths) > 0:
                seq_length[name] = {
                    "min": int(lengths.min()),
                    "median": int(lengths.median()),
                    "max": int(lengths.max()),
                    "mean": round(float(lengths.mean()), 1),
                }
            else:
                seq_length[name] = {}
        else:
            seq_length[name] = {}
    qc["seq_length"] = seq_length

    # --- Ion mobility stats per engine ---
    mobility_stats: Dict[str, Dict[str, Any]] = {}
    for name, df in engines.items():
        col = f"{name}_mobility"
        if df is not None and not df.empty and col in df.columns:
            vals = df[col].dropna()
            if len(vals) > 0:
                mobility_stats[name] = {
                    "min": round(float(vals.min()), 4),
                    "median": round(float(vals.median()), 4),
                    "max": round(float(vals.max()), 4),
                    "mean": round(float(vals.mean()), 4),
                }
            else:
                mobility_stats[name] = {}
        else:
            mobility_stats[name] = {}
    qc["mobility"] = mobility_stats

    # --- PEP score distribution per engine ---
    pep_stats: Dict[str, Dict[str, Any]] = {}
    for name, df in engines.items():
        col = f"{name}_pep"
        if df is not None and not df.empty and col in df.columns:
            vals = df[col].dropna()
            if len(vals) > 0:
                pep_stats[name] = {
                    "median": f"{float(vals.median()):.2e}",
                    "p95": f"{float(np.percentile(vals, 95)):.2e}",
                    "min": f"{float(vals.min()):.2e}",
                    "max": f"{float(vals.max()):.2e}",
                    "pct_below_001": round(float((vals < 0.01).mean()) * 100, 1),
                }
            else:
                pep_stats[name] = {}
        else:
            pep_stats[name] = {}
    qc["pep_scores"] = pep_stats

    # --- Per-raw-file PSM counts per engine ---
    per_file: Dict[str, Dict[str, int]] = {}
    for name, df in engines.items():
        if df is not None and not df.empty and "raw_file" in df.columns:
            counts = df.groupby("raw_file").size()
            per_file[name] = {str(k): int(v) for k, v in counts.items()}
        else:
            per_file[name] = {}
    qc["per_file_counts"] = per_file

    # --- Top modifications per engine ---
    mod_freq: Dict[str, List[Tuple[str, int]]] = {}
    import re
    for name, df in engines.items():
        col = f"{name}_modified"
        if df is not None and not df.empty and col in df.columns:
            seqs = df[col].dropna()
            mod_counts: Dict[str, int] = {}
            for seq in seqs:
                mods = re.findall(r'\[UNIMOD:\d+\]', str(seq))
                for m in mods:
                    mod_counts[m] = mod_counts.get(m, 0) + 1
            # Sort by frequency, take top 5
            top = sorted(mod_counts.items(), key=lambda x: -x[1])[:5]
            mod_freq[name] = top
        else:
            mod_freq[name] = []
    qc["top_mods"] = mod_freq

    return qc


def _generate_qc_plots(
    precursor_index: pd.DataFrame,
    max_points: int = 20_000,
) -> Dict[str, str]:
    """Generate QC scatter plots as base64-encoded PNGs.

    Produces:
    - 3 RT vs IM scatter plots (one per engine)
    - 3 pairwise IM alignment plots (engine pairs with y=x line and Pearson R)

    Args:
        precursor_index: Unified precursor index with per-engine columns
        max_points: Maximum points to plot (subsampled for performance)

    Returns:
        Dict mapping plot name to base64 PNG string
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import numpy as np
    from io import BytesIO
    import base64

    plots: Dict[str, str] = {}

    def _fig_to_base64(fig) -> str:
        buf = BytesIO()
        fig.savefig(buf, format='png', dpi=100, bbox_inches='tight')
        plt.close(fig)
        buf.seek(0)
        return base64.b64encode(buf.read()).decode('ascii')

    def _subsample(series_list, n):
        """Subsample multiple aligned series to at most n rows."""
        length = len(series_list[0])
        if length <= n:
            return series_list
        idx = np.random.default_rng(42).choice(length, size=n, replace=False)
        return [s.iloc[idx] for s in series_list]

    # --- RT vs IM scatter plots (one per engine) ---
    engines = {
        'FragPipe': ('fragpipe_rt', 'fragpipe_mobility'),
        'DIA-NN': ('diann_rt', 'diann_mobility'),
        'Sage': ('sage_rt', 'sage_mobility'),
    }
    colors = {'FragPipe': '#e74c3c', 'DIA-NN': '#3498db', 'Sage': '#27ae60'}

    for engine_label, (rt_col, im_col) in engines.items():
        if rt_col not in precursor_index.columns or im_col not in precursor_index.columns:
            continue
        mask = precursor_index[rt_col].notna() & precursor_index[im_col].notna()
        if mask.sum() == 0:
            continue

        rt = precursor_index.loc[mask, rt_col].astype(float)
        im = precursor_index.loc[mask, im_col].astype(float)
        rt, im = _subsample([rt, im], max_points)

        fig, ax = plt.subplots(figsize=(5, 4))
        ax.scatter(rt, im, s=1, alpha=0.15, c=colors[engine_label], rasterized=True)
        ax.set_xlabel('RT (seconds)')
        ax.set_ylabel('Ion Mobility (1/K0)')
        ax.set_title(f'{engine_label}: RT vs IM (n={mask.sum():,})')
        fig.tight_layout()
        plots[f'rt_im_{engine_label.lower().replace("-", "")}'] = _fig_to_base64(fig)

    # --- Pairwise IM alignment plots ---
    pairs = [
        ('FragPipe', 'DIA-NN', 'fragpipe_mobility', 'diann_mobility'),
        ('FragPipe', 'Sage', 'fragpipe_mobility', 'sage_mobility'),
        ('DIA-NN', 'Sage', 'diann_mobility', 'sage_mobility'),
    ]

    for engine_a, engine_b, col_a, col_b in pairs:
        if col_a not in precursor_index.columns or col_b not in precursor_index.columns:
            continue
        mask = precursor_index[col_a].notna() & precursor_index[col_b].notna()
        if mask.sum() < 2:
            continue

        im_a = precursor_index.loc[mask, col_a].astype(float)
        im_b = precursor_index.loc[mask, col_b].astype(float)

        # Compute Pearson R before subsampling (on full data)
        r = float(np.corrcoef(im_a, im_b)[0, 1])

        im_a_plot, im_b_plot = _subsample([im_a, im_b], max_points)

        fig, ax = plt.subplots(figsize=(5, 4.5))
        ax.scatter(im_a_plot, im_b_plot, s=1, alpha=0.15, c='#34495e', rasterized=True)

        # y=x reference line
        lo = min(float(im_a_plot.min()), float(im_b_plot.min()))
        hi = max(float(im_a_plot.max()), float(im_b_plot.max()))
        ax.plot([lo, hi], [lo, hi], 'r--', linewidth=1, label='y = x')

        ax.set_xlabel(f'{engine_a} IM (1/K0)')
        ax.set_ylabel(f'{engine_b} IM (1/K0)')
        ax.set_title(f'{engine_a} vs {engine_b} IM  (R={r:.4f}, n={mask.sum():,})')
        ax.legend(loc='lower right', fontsize=8)
        fig.tight_layout()

        key = f'im_align_{engine_a.lower().replace("-", "")}_{engine_b.lower().replace("-", "")}'
        plots[key] = _fig_to_base64(fig)

    return plots


def _generate_html_report(
    accession: str,
    stats: Dict[str, Any],
    output_path: Path,
    qc_stats: Optional[Dict[str, Any]] = None,
    plots: Optional[Dict[str, str]] = None,
) -> None:
    """Generate self-contained HTML overlap + QC report.

    Uses an inline template (no external import) with computed QC stats
    for charge distribution, sequence length, ion mobility, PEP scores,
    per-run PSM counts, and top modifications.
    """
    if qc_stats is None:
        qc_stats = {}
    if plots is None:
        plots = {}

    total = stats["n_union"]
    if total == 0:
        total = 1

    # --- Precompute derived values ---
    all_three_pct = stats["n_all_three"] / total * 100
    fp_dn_only_pct = stats["n_fp_dn_only"] / total * 100
    fp_sg_only_pct = stats["n_fp_sg_only"] / total * 100
    dn_sg_only_pct = stats["n_dn_sg_only"] / total * 100
    fp_only_pct = stats["n_fragpipe_only"] / total * 100
    dn_only_pct = stats["n_diann_only"] / total * 100
    sg_only_pct = stats["n_sage_only"] / total * 100

    fp_unique_pct = stats["n_fragpipe_only"] / max(stats["n_fragpipe"], 1) * 100
    dn_unique_pct = stats["n_diann_only"] / max(stats["n_diann"], 1) * 100
    sg_unique_pct = stats["n_sage_only"] / max(stats["n_sage"], 1) * 100

    fp_val = 1 - stats["n_fragpipe_only"] / max(stats["n_fragpipe"], 1)
    dn_val = 1 - stats["n_diann_only"] / max(stats["n_diann"], 1)
    sg_val = 1 - stats["n_sage_only"] / max(stats["n_sage"], 1)

    n_fp_dn = stats["n_all_three"] + stats["n_fp_dn_only"]
    n_fp_sg = stats["n_all_three"] + stats["n_fp_sg_only"]
    n_dn_sg = stats["n_all_three"] + stats["n_dn_sg_only"]

    # --- Build charge distribution rows ---
    charge_dist = qc_stats.get("charge_dist", {})
    all_charges = sorted(
        set(charge_dist.get("fragpipe", {}).keys())
        | set(charge_dist.get("diann", {}).keys())
        | set(charge_dist.get("sage", {}).keys())
    )
    if all_charges:
        charge_rows = "\n".join(
            f"        <tr><td>{c}+</td>"
            f"<td>{charge_dist.get('fragpipe', {}).get(c, 0):,}</td>"
            f"<td>{charge_dist.get('diann', {}).get(c, 0):,}</td>"
            f"<td>{charge_dist.get('sage', {}).get(c, 0):,}</td></tr>"
            for c in all_charges
        )
    else:
        charge_rows = "        <tr><td colspan='4'>No data</td></tr>"

    # --- Build sequence length rows ---
    sl = qc_stats.get("seq_length", {})

    def _sl(engine, stat):
        v = sl.get(engine, {}).get(stat)
        return str(v) if v is not None else "-"

    # --- Build ion mobility rows ---
    mob = qc_stats.get("mobility", {})

    def _mob(engine, stat):
        v = mob.get(engine, {}).get(stat)
        return str(v) if v is not None else "-"

    # --- Build PEP score rows ---
    pep = qc_stats.get("pep_scores", {})

    def _pep(engine, stat):
        v = pep.get(engine, {}).get(stat)
        return str(v) if v is not None else "-"

    # --- Build per-file PSM count table ---
    pf = qc_stats.get("per_file_counts", {})
    all_files = sorted(
        set(pf.get("fragpipe", {}).keys())
        | set(pf.get("diann", {}).keys())
        | set(pf.get("sage", {}).keys())
    )
    if all_files:
        per_file_rows = "\n".join(
            f"        <tr><td>{f}</td>"
            f"<td>{pf.get('fragpipe', {}).get(f, 0):,}</td>"
            f"<td>{pf.get('diann', {}).get(f, 0):,}</td>"
            f"<td>{pf.get('sage', {}).get(f, 0):,}</td></tr>"
            for f in all_files
        )
    else:
        per_file_rows = "        <tr><td colspan='4'>No data</td></tr>"

    # --- Build top modifications rows ---
    mods = qc_stats.get("top_mods", {})
    mod_names = {
        "[UNIMOD:4]": "Carbamidomethyl (C)",
        "[UNIMOD:35]": "Oxidation (M)",
        "[UNIMOD:1]": "Acetyl (N-term)",
        "[UNIMOD:21]": "Phospho (STY)",
        "[UNIMOD:7]": "Deamidation (NQ)",
        "[UNIMOD:34]": "Methylation",
        "[UNIMOD:28]": "Gln->pyro-Glu",
        "[UNIMOD:27]": "Glu->pyro-Glu",
        "[UNIMOD:5]": "Carbamyl",
        "[UNIMOD:122]": "Formyl (N-term)",
    }

    def _mod_table(engine):
        items = mods.get(engine, [])
        if not items:
            return "<em>-</em>"
        parts = []
        for m, count in items:
            label = mod_names.get(m, m)
            parts.append(f"{label}: {count:,}")
        return "<br>".join(parts)

    # --- Build plot sections ---
    rt_im_keys = ['rt_im_fragpipe', 'rt_im_diann', 'rt_im_sage']
    rt_im_imgs = [
        f'<img src="data:image/png;base64,{plots[k]}" alt="{k}">'
        for k in rt_im_keys if k in plots
    ]
    rt_im_section = ""
    if rt_im_imgs:
        rt_im_section = f"""
    <h2>RT vs Ion Mobility</h2>
    <div class="plot-row">
        {"".join(rt_im_imgs)}
    </div>"""

    im_align_keys = [
        'im_align_fragpipe_diann',
        'im_align_fragpipe_sage',
        'im_align_diann_sage',
    ]
    im_align_imgs = [
        f'<img src="data:image/png;base64,{plots[k]}" alt="{k}">'
        for k in im_align_keys if k in plots
    ]
    im_align_section = ""
    if im_align_imgs:
        im_align_section = f"""
    <h2>Ion Mobility Alignment (Cross-Engine)</h2>
    <div class="plot-row">
        {"".join(im_align_imgs)}
    </div>"""

    # --- Assemble HTML ---
    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>QC Report: {accession}</title>
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
        .plot-row {{
            display: flex;
            flex-wrap: wrap;
            gap: 15px;
            justify-content: center;
            margin: 20px 0;
        }}
        .plot-row img {{
            max-width: 380px;
            border: 1px solid var(--border-color);
            border-radius: 6px;
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
        <div><div class="label">Generated</div><div>San Jos&eacute; Runner</div></div>
        <div><div class="label">Match Type</div><div>sequence+charge</div></div>
        <div><div class="label">PEP Filter</div><div>None (all results)</div></div>
        <div><div class="label">I/L Normalized</div><div>Yes</div></div>
    </div>

    <h2>Search Engine Results</h2>
    <div class="summary-grid">
        <div class="stat-card fragpipe">
            <div class="number">{stats['n_fragpipe']:,}</div>
            <div class="label">FragPipe</div>
            <div class="pct">{fp_val:.1%} validated</div>
        </div>
        <div class="stat-card diann">
            <div class="number">{stats['n_diann']:,}</div>
            <div class="label">DIA-NN</div>
            <div class="pct">{dn_val:.1%} validated</div>
        </div>
        <div class="stat-card sage">
            <div class="number">{stats['n_sage']:,}</div>
            <div class="label">Sage</div>
            <div class="pct">{sg_val:.1%} validated</div>
        </div>
    </div>

    <h2>Consensus Summary</h2>
    <div class="summary-grid">
        <div class="stat-card all-three">
            <div class="number">{stats['n_all_three']:,}</div>
            <div class="label">All 3 Engines</div>
            <div class="pct">{stats['three_way_rate']:.1%} of union</div>
        </div>
        <div class="stat-card two-plus">
            <div class="number">{stats['n_at_least_two']:,}</div>
            <div class="label">At Least 2 Engines</div>
            <div class="pct">{stats['at_least_two_rate']:.1%} of union</div>
        </div>
        <div class="stat-card">
            <div class="number">{stats['n_union']:,}</div>
            <div class="label">Union (Any Engine)</div>
            <div class="pct">100%</div>
        </div>
    </div>

    <h2>Overlap Breakdown</h2>
    <div class="bar-stacked">
        <div class="bar-segment bar-all" style="width: {all_three_pct}%" title="All 3">{stats['n_all_three']:,}</div>
        <div class="bar-segment bar-fp-dn" style="width: {fp_dn_only_pct}%" title="FP+DN only">{stats['n_fp_dn_only']:,}</div>
        <div class="bar-segment bar-fp-sg" style="width: {fp_sg_only_pct}%" title="FP+Sage only">{stats['n_fp_sg_only']:,}</div>
        <div class="bar-segment bar-dn-sg" style="width: {dn_sg_only_pct}%" title="DN+Sage only">{stats['n_dn_sg_only']:,}</div>
        <div class="bar-segment bar-fp" style="width: {fp_only_pct}%" title="FragPipe only">{stats['n_fragpipe_only']:,}</div>
        <div class="bar-segment bar-dn" style="width: {dn_only_pct}%" title="DIA-NN only">{stats['n_diann_only']:,}</div>
        <div class="bar-segment bar-sg" style="width: {sg_only_pct}%" title="Sage only">{stats['n_sage_only']:,}</div>
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
        <tr><td>FragPipe &cap; DIA-NN</td><td>{n_fp_dn:,}</td><td>{stats['n_fp_dn_only']:,}</td></tr>
        <tr><td>FragPipe &cap; Sage</td><td>{n_fp_sg:,}</td><td>{stats['n_fp_sg_only']:,}</td></tr>
        <tr><td>DIA-NN &cap; Sage</td><td>{n_dn_sg:,}</td><td>{stats['n_dn_sg_only']:,}</td></tr>
    </table>

    <h2>Validation Quality</h2>
    <table class="validation-table">
        <tr><th>Engine</th><th>Unique (not in others)</th><th>Validation Rate</th></tr>
        <tr><td>FragPipe</td><td>{stats['n_fragpipe_only']:,} ({fp_unique_pct:.1f}%)</td><td>{fp_val:.1%}</td></tr>
        <tr><td>DIA-NN</td><td>{stats['n_diann_only']:,} ({dn_unique_pct:.1f}%)</td><td>{dn_val:.1%}</td></tr>
        <tr class="highlight"><td>Sage</td><td>{stats['n_sage_only']:,} ({sg_unique_pct:.1f}%)</td><td>{sg_val:.1%}</td></tr>
    </table>

    <h2>Charge Distribution</h2>
    <table>
        <tr><th>Charge</th><th>FragPipe</th><th>DIA-NN</th><th>Sage</th></tr>
{charge_rows}
    </table>

    <h2>Sequence Length</h2>
    <table>
        <tr><th>Statistic</th><th>FragPipe</th><th>DIA-NN</th><th>Sage</th></tr>
        <tr><td>Min</td><td>{_sl('fragpipe','min')}</td><td>{_sl('diann','min')}</td><td>{_sl('sage','min')}</td></tr>
        <tr><td>Median</td><td>{_sl('fragpipe','median')}</td><td>{_sl('diann','median')}</td><td>{_sl('sage','median')}</td></tr>
        <tr><td>Mean</td><td>{_sl('fragpipe','mean')}</td><td>{_sl('diann','mean')}</td><td>{_sl('sage','mean')}</td></tr>
        <tr><td>Max</td><td>{_sl('fragpipe','max')}</td><td>{_sl('diann','max')}</td><td>{_sl('sage','max')}</td></tr>
    </table>

    <h2>Ion Mobility (1/K0)</h2>
    <table>
        <tr><th>Statistic</th><th>FragPipe</th><th>DIA-NN</th><th>Sage</th></tr>
        <tr><td>Min</td><td>{_mob('fragpipe','min')}</td><td>{_mob('diann','min')}</td><td>{_mob('sage','min')}</td></tr>
        <tr><td>Median</td><td>{_mob('fragpipe','median')}</td><td>{_mob('diann','median')}</td><td>{_mob('sage','median')}</td></tr>
        <tr><td>Mean</td><td>{_mob('fragpipe','mean')}</td><td>{_mob('diann','mean')}</td><td>{_mob('sage','mean')}</td></tr>
        <tr><td>Max</td><td>{_mob('fragpipe','max')}</td><td>{_mob('diann','max')}</td><td>{_mob('sage','max')}</td></tr>
    </table>

    <h2>Score Quality (PEP)</h2>
    <table>
        <tr><th>Statistic</th><th>FragPipe</th><th>DIA-NN</th><th>Sage</th></tr>
        <tr><td>Median</td><td>{_pep('fragpipe','median')}</td><td>{_pep('diann','median')}</td><td>{_pep('sage','median')}</td></tr>
        <tr><td>95th Percentile</td><td>{_pep('fragpipe','p95')}</td><td>{_pep('diann','p95')}</td><td>{_pep('sage','p95')}</td></tr>
        <tr><td>Min</td><td>{_pep('fragpipe','min')}</td><td>{_pep('diann','min')}</td><td>{_pep('sage','min')}</td></tr>
        <tr><td>Max</td><td>{_pep('fragpipe','max')}</td><td>{_pep('diann','max')}</td><td>{_pep('sage','max')}</td></tr>
        <tr><td>% PEP &lt; 0.01</td><td>{_pep('fragpipe','pct_below_001')}</td><td>{_pep('diann','pct_below_001')}</td><td>{_pep('sage','pct_below_001')}</td></tr>
    </table>

    <h2>Per-Run PSM Counts</h2>
    <table>
        <tr><th>Raw File</th><th>FragPipe</th><th>DIA-NN</th><th>Sage</th></tr>
{per_file_rows}
    </table>

    <h2>Top Modifications</h2>
    <table>
        <tr><th>Engine</th><th>Most Frequent Modifications</th></tr>
        <tr><td>FragPipe</td><td>{_mod_table('fragpipe')}</td></tr>
        <tr><td>DIA-NN</td><td>{_mod_table('diann')}</td></tr>
        <tr><td>Sage</td><td>{_mod_table('sage')}</td></tr>
    </table>
{rt_im_section}
{im_align_section}

    <div class="footer">
        <p>Generated by <strong>San Jos&eacute; Pipeline</strong> | Triple Orthogonal Validation</p>
        <p>Using I/L normalization and UNIMOD standardization</p>
    </div>
</body>
</html>"""

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
    if "n_total_precursors" in summary.data:
        print(f"  Total precursors: {summary.data['n_total_precursors']}")
    if "groups" in summary.data:
        for gid, gdata in summary.data["groups"].items():
            status = gdata.get("status", "unknown")
            n = gdata.get("n_total_precursors", "N/A")
            print(f"  Group {gid}: {status}, {n} precursors")
    if "stratified_counts" in summary.data:
        print(f"  Stratified counts: {summary.data['stratified_counts']}")
