#!/usr/bin/env python3
"""
San Jose Smoke Test - Pipeline execution + ground truth validation.

Runs the full 6-step pipeline on 3 simulated .d files, then validates
pipeline outputs against simulation ground truth. Generates an HTML report.

Prerequisites:
    .venv/bin/python scripts/smoke_test_setup.py   # One-time: generate simulated datasets

Usage:
    .venv/bin/python scripts/smoke_test.py                # Full run
    .venv/bin/python scripts/smoke_test.py --validate-only # Skip pipeline, validate existing outputs
    .venv/bin/python scripts/smoke_test.py --skip-package  # Skip step 6
    .venv/bin/python scripts/smoke_test.py --clean         # Clean outputs, keep sims + DIA-NN lib
    .venv/bin/python scripts/smoke_test.py --clean-all     # Clean all outputs (keep sims)
"""

import argparse
import base64
import json
import shutil
import sys
import time
from datetime import datetime
from io import BytesIO
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import yaml

# Add parent to path for runner imports
sys.path.insert(0, str(Path(__file__).parent.parent))

ACCESSION = "SIM_SMOKE"
SMOKE_BASE = Path("data/smoke_test")
SIM_DIR = SMOKE_BASE / "simulations"
VALIDATION_DIR = SMOKE_BASE / "validation"
REPLICATE_NAMES = ["SMOKE_REP_01", "SMOKE_REP_02", "SMOKE_REP_03"]
ENGINES = ["fragpipe", "diann", "sage"]


# ---------------------------------------------------------------------------
# Phase A: Pipeline execution
# ---------------------------------------------------------------------------

def check_simulations_exist() -> bool:
    """Verify that all 3 simulated .d files are present."""
    for name in REPLICATE_NAMES:
        d_path = SIM_DIR / name / name / f"{name}.d"
        db_path = SIM_DIR / name / name / "synthetic_data.db"
        if not d_path.exists() or not db_path.exists():
            return False
    return True


def create_source_symlinks() -> Path:
    """Create a directory with symlinks to the 3 simulated .d files."""
    source_dir = SMOKE_BASE / "raw_source"
    source_dir.mkdir(parents=True, exist_ok=True)

    for name in REPLICATE_NAMES:
        src = SIM_DIR / name / name / f"{name}.d"
        dst = source_dir / f"{name}.d"
        if dst.exists() or dst.is_symlink():
            dst.unlink() if dst.is_symlink() else shutil.rmtree(dst)
        dst.symlink_to(src.resolve())

    return source_dir


def build_config() -> Dict[str, Any]:
    """Load base config and apply smoke test overrides."""
    config_path = Path("config/config.yaml")
    with open(config_path) as f:
        config = yaml.safe_load(f)

    config["_config_path"] = str(config_path)

    source_dir = create_source_symlinks()

    # Override local_data to point to our symlinked simulations
    config.setdefault("local_data", {})[ACCESSION] = str(source_dir)

    # Set organism FASTA for the simulated data (HeLa)
    config.setdefault("organisms", {}).setdefault("human", {})
    config["organisms"]["human"]["local_fasta"] = str(
        Path.home() / "validate-sim/sim/hela.fasta"
    )
    config["organisms"]["human"]["includes_contaminants"] = True

    # Disable paper extraction
    config.setdefault("paper_extraction", {})["enabled"] = False

    # Dataset metadata
    config.setdefault("dataset_metadata", {})[ACCESSION] = {
        "description": "Smoke test - 3 simulated HeLa DDA replicates",
        "organism": "human",
        "instrument": "timsTOF Pro (simulated)",
        "acquisition_mode": "PASEF",
    }

    return config


def run_pipeline(config: Dict[str, Any], skip_package: bool = False) -> Dict[str, float]:
    """Run the San Jose pipeline on simulated data. Returns per-step timings."""
    from runner.run_dataset import run_dataset

    timings = {}
    t0 = time.time()

    success = run_dataset(
        accession=ACCESSION,
        config=config,
        output_base_dir=SMOKE_BASE,
        test_mode=False,
        max_files=0,
        resume=False,
        steps=None,
        num_threads=16,
        local_data_path=None,  # Handled via config local_data
        package=not skip_package,
        package_version="0.0",
    )

    timings["total"] = time.time() - t0

    if not success:
        print("\nPipeline FAILED. Cannot proceed to validation.")
        sys.exit(1)

    # Read step timings from checkpoint
    checkpoint_path = SMOKE_BASE / "checkpoints" / ACCESSION / "state.json"
    if checkpoint_path.exists():
        with open(checkpoint_path) as f:
            state = json.load(f)
        for step_name, step_data in state.get("steps", {}).items():
            if step_data.get("duration_seconds"):
                timings[step_name] = step_data["duration_seconds"]

    return timings


# ---------------------------------------------------------------------------
# Phase B: Ground truth validation
# ---------------------------------------------------------------------------

def load_ground_truth_per_file() -> Dict[str, pd.DataFrame]:
    """Load ground truth from each simulation's synthetic_data.db."""
    from imspy_simulation.timsim.validate.comparison import load_ground_truth

    gt_per_file = {}
    for name in REPLICATE_NAMES:
        db_path = str(SIM_DIR / name / name / "synthetic_data.db")
        gt = load_ground_truth(db_path, experiment_name=name)
        gt_per_file[name] = gt
    return gt_per_file


def load_pipeline_results() -> Optional[pd.DataFrame]:
    """Load precursor_index.parquet from pipeline output."""
    # Per-group mode: look in human_trypsin subdirectory
    candidates = [
        SMOKE_BASE / "processed" / ACCESSION / "human_trypsin" / "precursor_index.parquet",
        SMOKE_BASE / "processed" / ACCESSION / "precursor_index.parquet",
    ]
    for path in candidates:
        if path.exists():
            return pd.read_parquet(path)

    # Try glob fallback
    found = list((SMOKE_BASE / "processed" / ACCESSION).rglob("precursor_index.parquet"))
    if found:
        return pd.read_parquet(found[0])

    return None


def normalize_for_matching(seq: str) -> str:
    """Normalize a sequence for matching: strip UNIMOD, I->L."""
    from imspy_core.utility import remove_unimod_annotation
    from imspy_simulation.timsim.validate.parsing import replace_I_with_L

    plain = remove_unimod_annotation(seq) if "[" in seq else seq
    return replace_I_with_L(plain)


def match_engine_vs_ground_truth(
    index: pd.DataFrame,
    gt: pd.DataFrame,
    engine: str,
) -> Dict[str, Any]:
    """Match a single engine's identifications against ground truth.

    Returns metrics dict with id_rate, precision, rt_r, im_r, counts.
    """
    from scipy import stats

    peptide_col = f"{engine}_peptide"
    modified_col = f"{engine}_modified"

    # Filter to rows where this engine identified something
    if peptide_col not in index.columns:
        return {"identified": 0, "matched": 0, "id_rate": 0.0, "precision": 0.0,
                "ground_truth": 0, "true_positives": 0, "false_positives": 0}
    engine_ids = index[index[peptide_col].notna()].copy()
    if engine_ids.empty:
        return {"identified": 0, "matched": 0, "id_rate": 0.0, "precision": 0.0,
                "ground_truth": 0, "true_positives": 0, "false_positives": 0}

    # Build match keys: (normalized_plain_sequence, charge)
    engine_ids["match_key"] = engine_ids[modified_col].apply(normalize_for_matching) + "_" + engine_ids["charge"].astype(str)
    gt_copy = gt.copy()
    gt_copy["match_key"] = gt_copy["sequence_normalized"] + "_" + gt_copy["charge"].astype(str)

    engine_keys = set(engine_ids["match_key"])
    gt_keys = set(gt_copy["match_key"])
    tp_keys = engine_keys & gt_keys
    fp_keys = engine_keys - gt_keys

    n_identified = len(engine_keys)
    n_gt = len(gt_keys)
    n_tp = len(tp_keys)
    n_fp = len(fp_keys)

    id_rate = n_tp / n_gt if n_gt > 0 else 0.0
    precision = n_tp / n_identified if n_identified > 0 else 0.0

    result = {
        "ground_truth": n_gt,
        "identified": n_identified,
        "true_positives": n_tp,
        "false_positives": n_fp,
        "id_rate": id_rate,
        "precision": precision,
    }

    # RT correlation: pipeline rt_seconds vs ground truth rt (both in seconds)
    matched_engine = engine_ids[engine_ids["match_key"].isin(tp_keys)].drop_duplicates("match_key")
    matched_gt = gt_copy[gt_copy["match_key"].isin(tp_keys)].drop_duplicates("match_key")

    if len(matched_engine) >= 3 and len(matched_gt) >= 3:
        merged = matched_engine[["match_key", "rt_seconds", "mobility"]].merge(
            matched_gt[["match_key", "rt", "inverse_mobility"]],
            on="match_key",
        )

        if len(merged) >= 3:
            # RT correlation (both in seconds, convert to minutes for display)
            rt_pipe = merged["rt_seconds"].values / 60.0
            rt_gt = merged["rt"].values / 60.0
            valid_rt = ~(np.isnan(rt_pipe) | np.isnan(rt_gt))
            if valid_rt.sum() >= 3:
                r, _ = stats.pearsonr(rt_gt[valid_rt], rt_pipe[valid_rt])
                mae = np.mean(np.abs(rt_gt[valid_rt] - rt_pipe[valid_rt]))
                result["rt_pearson_r"] = float(r)
                result["rt_mae_min"] = float(mae)

            # IM correlation (both in 1/K0)
            im_pipe = merged["mobility"].values
            im_gt = merged["inverse_mobility"].values
            valid_im = ~(np.isnan(im_pipe) | np.isnan(im_gt))
            if valid_im.sum() >= 3:
                r, _ = stats.pearsonr(im_gt[valid_im], im_pipe[valid_im])
                mae = np.mean(np.abs(im_gt[valid_im] - im_pipe[valid_im]))
                result["im_pearson_r"] = float(r)
                result["im_mae"] = float(mae)

    return result


def validate_consensus(index: pd.DataFrame, gt_all: pd.DataFrame) -> Dict[str, Any]:
    """Validate consensus accuracy: what fraction of multi-engine hits are in ground truth?"""
    # Build aggregated ground truth match keys
    gt_all["match_key"] = gt_all["sequence_normalized"] + "_" + gt_all["charge"].astype(str)
    gt_keys = set(gt_all["match_key"])

    results = {}
    for tier_name, min_engines in [("1+", 1), ("2+", 2), ("3/3", 3)]:
        tier_rows = index[index["n_engines"] >= min_engines].copy()
        if tier_rows.empty:
            results[tier_name] = {"count": 0, "in_gt": 0, "accuracy": 0.0}
            continue

        tier_rows["match_key"] = tier_rows["sequence_normalized"].apply(normalize_for_matching) + "_" + tier_rows["charge"].astype(str)
        tier_keys = set(tier_rows["match_key"])
        in_gt = tier_keys & gt_keys

        results[tier_name] = {
            "count": len(tier_keys),
            "in_gt": len(in_gt),
            "accuracy": len(in_gt) / len(tier_keys) if tier_keys else 0.0,
        }

    return results


def run_validation(timings: Optional[Dict[str, float]] = None) -> Dict[str, Any]:
    """Run Phase B validation and return metrics dict."""
    print("\n" + "=" * 70)
    print("  Phase B: Ground Truth Validation")
    print("=" * 70)

    # Load data
    print("\nLoading ground truth...")
    gt_per_file = load_ground_truth_per_file()

    print("Loading pipeline results...")
    index = load_pipeline_results()
    if index is None:
        print("Error: No precursor_index.parquet found. Run pipeline first.")
        sys.exit(1)

    print(f"  Precursor index: {len(index)} rows")
    print(f"  Raw files in index: {sorted(index['raw_file'].unique())}")

    # Aggregate ground truth across all files for consensus validation
    gt_all = pd.concat(gt_per_file.values(), ignore_index=True).drop_duplicates(
        subset=["sequence_normalized", "charge"]
    )
    print(f"  Ground truth (unique precursors): {len(gt_all)}")

    all_metrics: Dict[str, Any] = {
        "timestamp": datetime.now().isoformat(),
        "accession": ACCESSION,
        "n_replicates": len(REPLICATE_NAMES),
        "timings": timings or {},
    }

    # Per-file, per-engine validation
    print("\nPer-file, per-engine validation:")
    per_file_results = {}
    for raw_file in sorted(index["raw_file"].unique()):
        # Find the matching ground truth (strip .d if present)
        gt_name = raw_file.replace(".d", "")
        if gt_name not in gt_per_file:
            print(f"  Warning: No ground truth for {raw_file}")
            continue

        gt = gt_per_file[gt_name]
        file_index = index[index["raw_file"] == raw_file]
        per_file_results[raw_file] = {}

        for engine in ENGINES:
            metrics = match_engine_vs_ground_truth(file_index, gt, engine)
            per_file_results[raw_file][engine] = metrics

            rt_r = metrics.get("rt_pearson_r", "N/A")
            im_r = metrics.get("im_pearson_r", "N/A")
            rt_str = f"RT R={rt_r:.3f}" if isinstance(rt_r, float) else "RT R=N/A"
            im_str = f"IM R={im_r:.3f}" if isinstance(im_r, float) else "IM R=N/A"
            print(f"  {raw_file} | {engine:10s} | ID rate={metrics['id_rate']:.1%} | "
                  f"Precision={metrics['precision']:.1%} | {rt_str} | {im_str}")

    all_metrics["per_file"] = per_file_results

    # Aggregate per-engine (across all files)
    print("\nAggregate per-engine:")
    per_engine_agg = {}
    for engine in ENGINES:
        agg = match_engine_vs_ground_truth(index, gt_all, engine)
        per_engine_agg[engine] = agg

        rt_r = agg.get("rt_pearson_r", "N/A")
        im_r = agg.get("im_pearson_r", "N/A")
        rt_str = f"RT R={rt_r:.3f}" if isinstance(rt_r, float) else "RT R=N/A"
        im_str = f"IM R={im_r:.3f}" if isinstance(im_r, float) else "IM R=N/A"
        print(f"  {engine:10s} | ID rate={agg['id_rate']:.1%} | "
              f"Precision={agg['precision']:.1%} | {rt_str} | {im_str}")

    all_metrics["per_engine"] = per_engine_agg

    # Consensus validation
    print("\nConsensus validation:")
    consensus = validate_consensus(index, gt_all)
    for tier, data in consensus.items():
        print(f"  {tier:5s} | {data['count']:6d} precursors | "
              f"{data['in_gt']:6d} in ground truth | accuracy={data['accuracy']:.1%}")
    all_metrics["consensus"] = consensus

    return all_metrics


# ---------------------------------------------------------------------------
# HTML Report generation
# ---------------------------------------------------------------------------

def _make_scatter_plot(
    x: np.ndarray,
    y: np.ndarray,
    xlabel: str,
    ylabel: str,
    title: str,
) -> str:
    """Create a scatter/hexbin plot and return base64-encoded PNG."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(5, 4.5))

    valid = ~(np.isnan(x) | np.isnan(y))
    x, y = x[valid], y[valid]

    if len(x) < 2:
        ax.text(0.5, 0.5, "Insufficient data", ha="center", va="center", transform=ax.transAxes)
    else:
        from scipy import stats
        r, _ = stats.pearsonr(x, y)
        mae = np.mean(np.abs(x - y))

        ax.hexbin(x, y, gridsize=40, cmap="viridis", mincnt=1, linewidths=0)
        lims = [min(x.min(), y.min()), max(x.max(), y.max())]
        ax.plot(lims, lims, "r--", alpha=0.7, linewidth=1)
        ax.text(0.05, 0.95, f"R = {r:.4f}\nMAE = {mae:.3f}\nN = {len(x)}",
                transform=ax.transAxes, ha="left", va="top", fontsize=9,
                bbox=dict(boxstyle="round", facecolor="white", alpha=0.8))

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    fig.tight_layout()

    buf = BytesIO()
    fig.savefig(buf, format="png", dpi=120)
    plt.close(fig)
    return base64.b64encode(buf.getvalue()).decode("utf-8")


def _make_bar_chart(
    labels: List[str],
    values: List[float],
    ylabel: str,
    title: str,
    colors: Optional[List[str]] = None,
    fmt: str = ".1%",
) -> str:
    """Create a bar chart and return base64-encoded PNG."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(5, 4))

    if colors is None:
        colors = ["#3498db", "#e74c3c", "#2ecc71"][:len(labels)]

    bars = ax.bar(labels, values, color=colors, edgecolor="black", linewidth=0.5)
    for bar, val in zip(bars, values):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.01,
                f"{val:{fmt}}", ha="center", va="bottom", fontsize=9)

    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.set_ylim(0, max(values) * 1.25 if values and max(values) > 0 else 1.0)
    fig.tight_layout()

    buf = BytesIO()
    fig.savefig(buf, format="png", dpi=120)
    plt.close(fig)
    return base64.b64encode(buf.getvalue()).decode("utf-8")


def generate_scatter_plots(index: pd.DataFrame, gt_all: pd.DataFrame) -> Dict[str, str]:
    """Generate RT and IM scatter plots per engine. Returns dict of base64 PNGs."""
    plots = {}

    for engine in ENGINES:
        peptide_col = f"{engine}_peptide"
        modified_col = f"{engine}_modified"

        if peptide_col not in index.columns:
            continue
        engine_ids = index[index[peptide_col].notna()].copy()
        if engine_ids.empty:
            continue

        engine_ids["match_key"] = engine_ids[modified_col].apply(normalize_for_matching) + "_" + engine_ids["charge"].astype(str)
        gt_copy = gt_all.copy()
        gt_copy["match_key"] = gt_copy["sequence_normalized"] + "_" + gt_copy["charge"].astype(str)

        tp_keys = set(engine_ids["match_key"]) & set(gt_copy["match_key"])
        matched_e = engine_ids[engine_ids["match_key"].isin(tp_keys)].drop_duplicates("match_key")
        matched_g = gt_copy[gt_copy["match_key"].isin(tp_keys)].drop_duplicates("match_key")

        merged = matched_e[["match_key", "rt_seconds", "mobility"]].merge(
            matched_g[["match_key", "rt", "inverse_mobility"]], on="match_key"
        )

        if len(merged) < 3:
            continue

        # RT scatter (in minutes)
        plots[f"{engine}_rt"] = _make_scatter_plot(
            merged["rt"].values / 60.0,
            merged["rt_seconds"].values / 60.0,
            "Ground Truth RT (min)",
            "Pipeline RT (min)",
            f"{engine.capitalize()} - RT Correlation",
        )

        # IM scatter
        plots[f"{engine}_im"] = _make_scatter_plot(
            merged["inverse_mobility"].values,
            merged["mobility"].values,
            "Ground Truth 1/K0",
            "Pipeline 1/K0",
            f"{engine.capitalize()} - IM Correlation",
        )

    return plots


def generate_html_report(metrics: Dict[str, Any], plots: Dict[str, str]) -> str:
    """Generate self-contained HTML validation report."""
    timestamp = metrics.get("timestamp", datetime.now().isoformat())
    timings = metrics.get("timings", {})
    per_engine = metrics.get("per_engine", {})
    consensus = metrics.get("consensus", {})
    per_file = metrics.get("per_file", {})

    # Bar chart: per-engine ID rates
    id_rate_chart = _make_bar_chart(
        [e.capitalize() for e in ENGINES],
        [per_engine.get(e, {}).get("id_rate", 0) for e in ENGINES],
        "Identification Rate",
        "Per-Engine Identification Rate",
    )

    # Bar chart: per-engine precision
    precision_chart = _make_bar_chart(
        [e.capitalize() for e in ENGINES],
        [per_engine.get(e, {}).get("precision", 0) for e in ENGINES],
        "Precision",
        "Per-Engine Precision",
    )

    # Bar chart: consensus accuracy
    consensus_labels = list(consensus.keys())
    consensus_acc = [consensus[k].get("accuracy", 0) for k in consensus_labels]
    consensus_chart = _make_bar_chart(
        consensus_labels, consensus_acc, "Accuracy", "Consensus Accuracy vs Ground Truth"
    )

    def _fmt(val, fmt=".3f"):
        if val is None or val == "N/A":
            return "N/A"
        if isinstance(val, float) and np.isnan(val):
            return "N/A"
        return f"{val:{fmt}}"

    def _pct(val):
        if val is None or val == "N/A":
            return "N/A"
        return f"{val:.1%}"

    # Build per-engine summary table
    engine_rows = ""
    for engine in ENGINES:
        e = per_engine.get(engine, {})
        engine_rows += f"""<tr>
            <td><strong>{engine.capitalize()}</strong></td>
            <td>{e.get('ground_truth', 'N/A')}</td>
            <td>{e.get('identified', 0)}</td>
            <td>{e.get('true_positives', 0)}</td>
            <td>{_pct(e.get('id_rate'))}</td>
            <td>{_pct(e.get('precision'))}</td>
            <td>{_fmt(e.get('rt_pearson_r'), '.4f')}</td>
            <td>{_fmt(e.get('rt_mae_min'))}</td>
            <td>{_fmt(e.get('im_pearson_r'), '.4f')}</td>
            <td>{_fmt(e.get('im_mae'), '.4f')}</td>
        </tr>"""

    # Build per-file table
    file_rows = ""
    for raw_file, engines in per_file.items():
        for engine in ENGINES:
            e = engines.get(engine, {})
            file_rows += f"""<tr>
                <td>{raw_file}</td>
                <td>{engine.capitalize()}</td>
                <td>{_pct(e.get('id_rate'))}</td>
                <td>{_pct(e.get('precision'))}</td>
                <td>{_fmt(e.get('rt_pearson_r'), '.4f')}</td>
                <td>{_fmt(e.get('im_pearson_r'), '.4f')}</td>
            </tr>"""

    # Build consensus table
    consensus_rows = ""
    for tier, data in consensus.items():
        consensus_rows += f"""<tr>
            <td>{tier}</td>
            <td>{data['count']}</td>
            <td>{data['in_gt']}</td>
            <td>{_pct(data['accuracy'])}</td>
        </tr>"""

    # Build scatter plot images
    plot_html = ""
    for key, b64 in sorted(plots.items()):
        plot_html += f"""
        <div class="plot-card">
            <img src="data:image/png;base64,{b64}" alt="{key}" />
        </div>"""

    # Timing info
    timing_html = ""
    if timings:
        timing_html = "<h3>Pipeline Timings</h3><table class='metrics-table'><tr><th>Step</th><th>Duration</th></tr>"
        for step, secs in sorted(timings.items()):
            m, s = divmod(secs, 60)
            timing_html += f"<tr><td>{step}</td><td>{int(m)}m {s:.1f}s</td></tr>"
        timing_html += "</table>"

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>San Jose Smoke Test Report</title>
<style>
body {{
    font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
    max-width: 1400px;
    margin: 0 auto;
    padding: 20px;
    background: #f5f5f5;
    color: #333;
    line-height: 1.5;
}}
h1 {{ border-bottom: 3px solid #2c3e50; padding-bottom: 10px; }}
h2 {{ color: #2c3e50; margin-top: 30px; border-bottom: 1px solid #ddd; padding-bottom: 5px; }}
.summary-box {{
    background: white;
    border-radius: 8px;
    padding: 20px;
    box-shadow: 0 2px 4px rgba(0,0,0,0.1);
    margin: 15px 0;
}}
.metrics-table {{
    width: 100%;
    border-collapse: collapse;
    margin: 10px 0;
    background: white;
    border-radius: 4px;
    overflow: hidden;
}}
.metrics-table th {{
    background: #2c3e50;
    color: white;
    padding: 10px 12px;
    text-align: left;
    font-size: 0.9em;
}}
.metrics-table td {{
    padding: 8px 12px;
    border-bottom: 1px solid #eee;
    font-size: 0.9em;
}}
.metrics-table tr:hover {{ background: #f0f7ff; }}
.plot-grid {{
    display: grid;
    grid-template-columns: repeat(auto-fill, minmax(420px, 1fr));
    gap: 15px;
    margin: 15px 0;
}}
.plot-card {{
    background: white;
    border-radius: 8px;
    padding: 10px;
    box-shadow: 0 2px 4px rgba(0,0,0,0.1);
    text-align: center;
}}
.plot-card img {{
    max-width: 100%;
    height: auto;
    border-radius: 4px;
}}
.chart-row {{
    display: grid;
    grid-template-columns: repeat(auto-fill, minmax(350px, 1fr));
    gap: 15px;
    margin: 15px 0;
}}
.chart-card {{
    background: white;
    border-radius: 8px;
    padding: 10px;
    box-shadow: 0 2px 4px rgba(0,0,0,0.1);
    text-align: center;
}}
.chart-card img {{
    max-width: 100%;
    height: auto;
}}
footer {{
    margin-top: 40px;
    padding-top: 15px;
    border-top: 1px solid #ddd;
    color: #888;
    font-size: 0.85em;
}}
</style>
</head>
<body>

<h1>San Jose Smoke Test Report</h1>

<div class="summary-box">
    <strong>Accession:</strong> {ACCESSION} |
    <strong>Replicates:</strong> {len(REPLICATE_NAMES)} |
    <strong>Engines:</strong> FragPipe, DIA-NN, Sage |
    <strong>Generated:</strong> {timestamp}
</div>

{timing_html}

<h2>Per-Engine Summary (Aggregate)</h2>
<table class="metrics-table">
<tr>
    <th>Engine</th><th>Ground Truth</th><th>Identified</th><th>True Pos.</th>
    <th>ID Rate</th><th>Precision</th>
    <th>RT R</th><th>RT MAE (min)</th>
    <th>IM R</th><th>IM MAE</th>
</tr>
{engine_rows}
</table>

<h2>Identification Rates &amp; Precision</h2>
<div class="chart-row">
    <div class="chart-card">
        <img src="data:image/png;base64,{id_rate_chart}" alt="ID Rate" />
    </div>
    <div class="chart-card">
        <img src="data:image/png;base64,{precision_chart}" alt="Precision" />
    </div>
    <div class="chart-card">
        <img src="data:image/png;base64,{consensus_chart}" alt="Consensus" />
    </div>
</div>

<h2>RT and IM Correlation Plots</h2>
<div class="plot-grid">
{plot_html}
</div>

<h2>Consensus Validation</h2>
<table class="metrics-table">
<tr><th>Tier</th><th>Precursors</th><th>In Ground Truth</th><th>Accuracy</th></tr>
{consensus_rows}
</table>

<h2>Per-File Breakdown</h2>
<table class="metrics-table">
<tr><th>Raw File</th><th>Engine</th><th>ID Rate</th><th>Precision</th><th>RT R</th><th>IM R</th></tr>
{file_rows}
</table>

<footer>
    Generated by San Jose Smoke Test | {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}
</footer>

</body>
</html>"""

    return html


# ---------------------------------------------------------------------------
# Cleanup
# ---------------------------------------------------------------------------

def clean_outputs(keep_diann_lib: bool = True) -> None:
    """Remove pipeline outputs, keep simulations."""
    dirs_to_clean = [
        SMOKE_BASE / "raw" / ACCESSION,
        SMOKE_BASE / "raw_source",
        SMOKE_BASE / "extracted" / ACCESSION,
        SMOKE_BASE / "merged" / ACCESSION,
        SMOKE_BASE / "checkpoints" / ACCESSION,
        SMOKE_BASE / "metadata" / ACCESSION,
        SMOKE_BASE / "validation",
    ]

    processed_dir = SMOKE_BASE / "processed" / ACCESSION

    if keep_diann_lib:
        # Preserve DIA-NN library files, remove everything else in processed/
        if processed_dir.exists():
            for item in processed_dir.iterdir():
                if item.is_dir():
                    for sub in item.iterdir():
                        if sub.is_dir() and sub.name == "diann":
                            # Keep .speclib files
                            for f in sub.iterdir():
                                if not f.name.endswith(".speclib"):
                                    if f.is_dir():
                                        shutil.rmtree(f)
                                    else:
                                        f.unlink()
                        else:
                            if sub.is_dir():
                                shutil.rmtree(sub)
                            else:
                                sub.unlink()
                else:
                    item.unlink()
    else:
        dirs_to_clean.append(processed_dir)

    for d in dirs_to_clean:
        if d.exists():
            shutil.rmtree(d)
            print(f"  Removed: {d}")

    # Clean packages
    packages_dir = SMOKE_BASE / "packages"
    if packages_dir.exists():
        for f in packages_dir.glob(f"{ACCESSION}*"):
            f.unlink()
            print(f"  Removed: {f}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="San Jose Smoke Test")
    parser.add_argument("--validate-only", action="store_true",
                        help="Skip pipeline, validate existing outputs")
    parser.add_argument("--skip-package", action="store_true",
                        help="Skip step 6 (packaging)")
    parser.add_argument("--clean", action="store_true",
                        help="Clean outputs, keep simulations + DIA-NN library")
    parser.add_argument("--clean-all", action="store_true",
                        help="Clean all outputs (keep simulations)")
    args = parser.parse_args()

    if args.clean or args.clean_all:
        print("Cleaning smoke test outputs...")
        clean_outputs(keep_diann_lib=not args.clean_all)
        print("Done.")
        return

    # Check prerequisites
    if not check_simulations_exist():
        print("Error: Simulated datasets not found.")
        print("Run:  .venv/bin/python scripts/smoke_test_setup.py")
        sys.exit(1)

    timings = None

    # Phase A: Pipeline execution
    if not args.validate_only:
        print("=" * 70)
        print("  Phase A: Pipeline Execution")
        print("=" * 70)

        config = build_config()
        timings = run_pipeline(config, skip_package=args.skip_package)

    # Phase B: Validation
    all_metrics = run_validation(timings)

    # Generate scatter plots
    print("\nGenerating scatter plots...")
    gt_per_file = load_ground_truth_per_file()
    gt_all = pd.concat(gt_per_file.values(), ignore_index=True).drop_duplicates(
        subset=["sequence_normalized", "charge"]
    )
    index = load_pipeline_results()
    scatter_plots = generate_scatter_plots(index, gt_all) if index is not None else {}

    # Write outputs
    VALIDATION_DIR.mkdir(parents=True, exist_ok=True)

    # JSON metrics
    json_path = VALIDATION_DIR / "smoke_test_report.json"
    with open(json_path, "w") as f:
        json.dump(all_metrics, f, indent=2, default=str)
    print(f"\nMetrics saved to: {json_path}")

    # HTML report
    html_content = generate_html_report(all_metrics, scatter_plots)
    html_path = VALIDATION_DIR / "smoke_test_report.html"
    with open(html_path, "w") as f:
        f.write(html_content)
    print(f"HTML report saved to: {html_path}")

    # Save plots as separate files too
    plots_dir = VALIDATION_DIR / "plots"
    plots_dir.mkdir(exist_ok=True)
    for key, b64 in scatter_plots.items():
        plot_path = plots_dir / f"{key}.png"
        with open(plot_path, "wb") as f:
            f.write(base64.b64decode(b64))

    print(f"\nSmoke test complete. Open {html_path} to view the report.")


if __name__ == "__main__":
    main()
