#!/usr/bin/env python3
"""
Step 3: Stratify and Merge Third-Party Output

UNIMOD-standardizes sequences, computes consensus, stratifies by engine agreement.

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


def run_step3_stratify(
    accession: str,
    config: Dict[str, Any],
    output_base_dir: Path,
    generate_html: bool = True,
) -> StepSummary:
    """
    Execute Step 3: Stratify and merge search engine results.

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
    consensus_dir = processed_dir / "consensus"
    stratified_dir = consensus_dir / "stratified"

    try:
        # Load search engine results
        print("  Loading search engine results...")
        fp_df = _load_fragpipe(processed_dir)
        dn_df = _load_diann(processed_dir)
        sg_df = _load_sage(processed_dir)

        print(f"    FragPipe: {len(fp_df) if fp_df is not None else 0} PSMs")
        print(f"    DIA-NN: {len(dn_df) if dn_df is not None else 0} precursors")
        print(f"    Sage: {len(sg_df) if sg_df is not None else 0} PSMs")

        # Standardize sequences to UNIMOD format
        print("  Standardizing sequences to UNIMOD format...")
        if fp_df is not None:
            fp_df = _standardize_fragpipe(fp_df)
        if dn_df is not None:
            dn_df = _standardize_diann(dn_df)
        if sg_df is not None:
            sg_df = _standardize_sage(sg_df)

        # Compute overlap statistics
        print("  Computing overlap statistics...")
        overlap_stats = _compute_overlap_stats(fp_df, dn_df, sg_df)

        # Build unified precursor index
        print("  Building unified precursor index...")
        precursor_index = _build_precursor_index(fp_df, dn_df, sg_df, processed_dir)

        # Stratify by engine agreement
        print("  Stratifying by engine agreement...")
        stratified_counts = _stratify_precursors(precursor_index, stratified_dir)

        # Compute match tier statistics
        match_tiers = _compute_match_tiers(precursor_index)

        # Save precursor index
        index_path = processed_dir / "precursor_index.parquet"
        precursor_index.to_parquet(index_path, index=False)
        print(f"  Saved precursor index: {index_path}")

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


def _load_fragpipe(processed_dir: Path) -> Optional[pd.DataFrame]:
    """Load FragPipe results."""
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
    """Load DIA-NN results."""
    report_path = processed_dir / "diann" / "report.parquet"
    if report_path.exists():
        return pd.read_parquet(report_path)

    report_tsv = processed_dir / "diann" / "report.tsv"
    if report_tsv.exists():
        return pd.read_csv(report_tsv, sep="\t")

    return None


def _load_sage(processed_dir: Path) -> Optional[pd.DataFrame]:
    """Load Sage results."""
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

    if mod_col:
        df["modified_unimod"] = df[mod_col].apply(
            standardize_fragpipe_modified_peptide
        )
    elif "Peptide Sequence" in df.columns:
        df["modified_unimod"] = df["Peptide Sequence"]
    elif "Peptide" in df.columns:
        df["modified_unimod"] = df["Peptide"]
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


def _compute_overlap_stats(
    fp_df: Optional[pd.DataFrame],
    dn_df: Optional[pd.DataFrame],
    sg_df: Optional[pd.DataFrame],
) -> Dict[str, Any]:
    """Compute 3-way overlap statistics."""

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


def _build_precursor_index(
    fp_df: Optional[pd.DataFrame],
    dn_df: Optional[pd.DataFrame],
    sg_df: Optional[pd.DataFrame],
    processed_dir: Path,
) -> pd.DataFrame:
    """Build unified precursor index with all engine IDs using efficient merges.

    Schema (37 columns total):
    - Join keys: sequence_normalized, charge
    - FragPipe (10): peptide, modified, protein, probability, pep, hyperscore, qvalue, rt, mz, mobility
    - DIA-NN (13): peptide, modified, protein, qvalue, pep, global_qvalue, pg_qvalue, rt, mz, mobility, ccs
    - Sage (12): peptide, modified, protein, qvalue, pep, hyperscore, peptide_qvalue, protein_qvalue, rt, mz, mobility
    - Derived: n_engines
    """

    # Prepare each engine's data with consistent column names
    dfs_to_merge = []

    if fp_df is not None and not fp_df.empty:
        fp_summary = fp_df.groupby(["sequence_normalized", "Charge"]).first().reset_index()
        fp_summary = fp_summary.rename(columns={"Charge": "charge"})

        # Core identification columns
        fp_summary["fragpipe_modified"] = fp_summary.get("modified_unimod", "")
        fp_summary["fragpipe_peptide"] = fp_summary.get("Peptide Sequence", fp_summary.get("Peptide", ""))
        fp_summary["fragpipe_protein"] = fp_summary.get("Protein", "")

        # Quality/confidence scores
        fp_summary["fragpipe_probability"] = fp_summary.get("Probability")
        fp_summary["fragpipe_pep"] = 1.0 - fp_summary.get("Probability", 0.0)  # PEP = 1 - Probability
        fp_summary["fragpipe_hyperscore"] = fp_summary.get("Hyperscore")
        fp_summary["fragpipe_qvalue"] = fp_summary.get("Qvalue")

        # Coordinates (RT already in seconds)
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

        # Core identification columns
        dn_summary["diann_modified"] = dn_summary.get("modified_unimod", "")
        dn_summary["diann_peptide"] = dn_summary.get("Stripped.Sequence", "")
        dn_summary["diann_protein"] = dn_summary.get("Protein.Ids", "")

        # Quality/confidence scores
        dn_summary["diann_qvalue"] = dn_summary.get("Q.Value")
        dn_summary["diann_pep"] = dn_summary.get("PEP")
        dn_summary["diann_global_qvalue"] = dn_summary.get("Global.Q.Value")
        dn_summary["diann_pg_qvalue"] = dn_summary.get("PG.Q.Value")

        # Coordinates (convert RT from minutes to seconds)
        rt_minutes = dn_summary.get("RT")
        if rt_minutes is not None:
            dn_summary["diann_rt"] = rt_minutes * 60.0
        else:
            dn_summary["diann_rt"] = None
        dn_summary["diann_mz"] = dn_summary.get("Precursor.Mz")
        dn_summary["diann_mobility"] = dn_summary.get("IM")

        # DIA-NN specific: predicted CCS
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

        # Core identification columns
        sg_summary["sage_modified"] = sg_summary.get("modified_unimod", "")
        sg_summary["sage_peptide"] = sg_summary.get("stripped_peptide", "")
        sg_summary["sage_protein"] = sg_summary.get("proteins", "")

        # Quality/confidence scores
        sg_summary["sage_qvalue"] = sg_summary.get("spectrum_q")
        sg_summary["sage_pep"] = sg_summary.get("posterior_error").apply(
            lambda x: np.exp(x) if pd.notna(x) else np.nan
        ) if "posterior_error" in sg_summary.columns else np.nan
        sg_summary["sage_hyperscore"] = sg_summary.get("hyperscore")
        sg_summary["sage_peptide_qvalue"] = sg_summary.get("peptide_q")
        sg_summary["sage_protein_qvalue"] = sg_summary.get("protein_q")

        # Coordinates (RT already in seconds)
        sg_summary["sage_rt"] = sg_summary.get("rt")
        # Calculate m/z from experimental mass: mz = (M + z*H+) / z
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

    # Start with first dataframe and outer merge the rest
    _, merged = dfs_to_merge[0]
    for name, df in dfs_to_merge[1:]:
        merged = merged.merge(df, on=["sequence_normalized", "charge"], how="outer")

    # Count engines for each row
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
