#!/usr/bin/env python3
"""
CCS Prediction Accuracy vs Cross-Sample Reproducibility

Compares predicted 1/K0 (from local DeepPeptideIonMobilityApex model) against
observed per-peptide mean 1/K0. The key question: how does prediction error
compare to within-dataset cross-sample reproducibility?

Hypothesis: The predictor was trained on old timsTOF data (Meier 2021), so
prediction error may be ~1 order of magnitude worse than cross-sample CV%.

Output:
    - ccs_prediction_accuracy.parquet: per-peptide predicted vs observed
    - ccs_prediction_accuracy_report.pdf: diagnostic plots
"""

import argparse
import re
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)

RUSTIMS_PREDICTORS = Path("/home/administrator/Documents/promotion/rust/rustims/packages/imspy-predictors/src")


def clean_sequence(seq: str) -> str | None:
    """Clean sequence for the CCS predictor (all UNIMOD mods supported)."""
    if not seq or not isinstance(seq, str):
        return None
    cleaned = seq.replace("-", "").replace("[]", "")
    stripped = re.sub(r"\[UNIMOD:\d+\]", "", cleaned)
    if len(stripped) < 7 or len(stripped) > 30:
        return None
    return cleaned


# CCS model charge embedding only covers 2-4 (trained on Meier 2021 data)
_VALID_CHARGES = {2, 3, 4}


def load_peptide_mobilities(
    store_path: str,
    repro_path: str | None = None,
    min_files: int = 3,
) -> pd.DataFrame:
    """Load per-peptide mean mobility and mz from precursor_store.

    If repro_path is given, merges aligned mobility and CV% from the
    CCS reproducibility analysis.
    """
    cols = ["sequence_normalized", "charge", "mobility", "mono_mz", "raw_file",
            "n_engines"]
    pf = pq.ParquetFile(store_path)
    df = pf.read(columns=cols).to_pandas()

    # Filter to identified precursors with valid mobility
    mask = (df["n_engines"] >= 1) & df["mobility"].notna() & (df["mobility"] > 0)
    df = df.loc[mask].copy()

    # Per-peptide aggregation
    agg = df.groupby(["sequence_normalized", "charge"]).agg(
        mean_mobility=("mobility", "mean"),
        mean_mz=("mono_mz", "mean"),
        n_files=("raw_file", "nunique"),
        n_obs=("mobility", "count"),
    ).reset_index()

    agg = agg[agg["n_files"] >= min_files]

    # Merge with CCS reproducibility if available
    if repro_path and Path(repro_path).exists():
        repro = pd.read_parquet(repro_path)
        merge_cols = ["sequence_normalized", "charge"]
        extra = ["cv_pct", "cv_pct_aligned", "mean_mobility_aligned"]
        extra = [c for c in extra if c in repro.columns]
        agg = agg.merge(repro[merge_cols + extra], on=merge_cols, how="left")
        print(f"  Merged CCS reproducibility for {agg['cv_pct'].notna().sum():,} peptides",
              flush=True)

    return agg


def predict_mobilities(
    df: pd.DataFrame,
    batch_size: int = 2048,
) -> pd.Series:
    """Predict 1/K0 for each peptide using DeepPeptideIonMobilityApex."""
    sys.path.insert(0, str(RUSTIMS_PREDICTORS))
    from imspy_predictors.ccs.predictors import DeepPeptideIonMobilityApex

    predictor = DeepPeptideIonMobilityApex()
    print(f"  CCS model loaded", flush=True)

    # Clean sequences and filter to supported charges (2-4)
    seqs_clean = df["sequence_normalized"].apply(clean_sequence)
    charge_ok = df["charge"].isin(_VALID_CHARGES)
    valid_mask = seqs_clean.notna() & charge_ok
    valid_df = df.loc[valid_mask]
    valid_seqs = seqs_clean.loc[valid_mask].tolist()

    print(f"  {len(valid_seqs):,} / {len(df):,} peptides with valid sequences",
          flush=True)

    predicted = pd.Series(np.nan, index=df.index)

    if not valid_seqs:
        return predicted

    charges = valid_df["charge"].values.astype(np.int64)
    mzs = valid_df["mean_mz"].values.astype(np.float64)

    inv_mob = predictor.simulate_ion_mobilities(
        sequences=valid_seqs,
        charges=charges,
        mz=mzs,
        batch_size=batch_size,
    )

    predicted.loc[valid_mask] = inv_mob
    return predicted


def generate_report(df: pd.DataFrame, output_path: Path):
    """Generate CCS prediction accuracy PDF report."""
    plt.rcParams.update({
        "figure.dpi": 150, "font.size": 9,
        "axes.titlesize": 11, "axes.labelsize": 10,
        "figure.figsize": (8, 6),
    })

    has_cv = "cv_pct_aligned" in df.columns and df["cv_pct_aligned"].notna().any()
    has_cv_raw = "cv_pct" in df.columns and df["cv_pct"].notna().any()

    with PdfPages(str(output_path)) as pdf:
        # --- Page 1: Predicted (aligned) vs Observed scatter ---
        fig, ax = plt.subplots()
        valid = df.dropna(subset=["predicted_aligned"])
        hb = ax.hexbin(valid["mean_mobility_obs"], valid["predicted_aligned"],
                       gridsize=80, cmap="YlGnBu", mincnt=1)
        fig.colorbar(hb, ax=ax, label="Count")
        lims = [valid["mean_mobility_obs"].min() * 0.95,
                valid["mean_mobility_obs"].max() * 1.05]
        ax.plot(lims, lims, "r--", alpha=0.5, linewidth=1, label="y = x")
        ax.set_xlabel("Observed 1/K0 (Vs/cm²)")
        ax.set_ylabel("Predicted 1/K0 (Vs/cm², aligned)")
        r2 = np.corrcoef(valid["mean_mobility_obs"], valid["predicted_aligned"])[0, 1] ** 2
        mae = valid["abs_error"].median()
        ax.set_title(f"Predicted vs Observed 1/K0 (aligned) — R² = {r2:.4f}, median |error| = {mae:.4f}")
        ax.legend()
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)

        # --- Page 2: Relative error distribution vs CV% ---
        fig, ax = plt.subplots()
        rel_err = valid["rel_error_pct"].dropna()
        ax.hist(rel_err, bins=200, alpha=0.6, density=True,
                label=f"Prediction error (med={rel_err.median():.2f}%)", color="#dc2626")
        if has_cv:
            cv_vals = df["cv_pct_aligned"].dropna()
            ax.hist(cv_vals, bins=200, alpha=0.6, density=True,
                    label=f"Cross-sample CV% aligned (med={cv_vals.median():.2f}%)",
                    color="#2563eb")
        elif has_cv_raw:
            cv_vals = df["cv_pct"].dropna()
            ax.hist(cv_vals, bins=200, alpha=0.6, density=True,
                    label=f"Cross-sample CV% (med={cv_vals.median():.2f}%)",
                    color="#2563eb")
        ax.set_xlabel("Error / CV (%)")
        ax.set_ylabel("Density")
        ax.set_xlim(0, min(rel_err.quantile(0.99), 20))
        ax.set_title("CCS Prediction Error vs Cross-Sample Reproducibility")
        ax.legend(fontsize=8)
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)

        # --- Page 3: Relative error by charge state ---
        fig, ax = plt.subplots()
        charge_vals = sorted(valid["charge"].unique())
        charge_vals = [c for c in charge_vals if 2 <= c <= 5]
        data_by_charge = [
            valid[valid["charge"] == c]["rel_error_pct"].values
            for c in charge_vals
        ]
        if data_by_charge:
            parts = ax.violinplot(data_by_charge, positions=range(len(charge_vals)),
                                  showmedians=True, showextrema=False)
            for pc in parts["bodies"]:
                pc.set_facecolor("#6366f1")
                pc.set_alpha(0.6)
            ax.set_xticks(range(len(charge_vals)))
            ax.set_xticklabels([f"{c}+" for c in charge_vals])
            for i, c in enumerate(charge_vals):
                med = np.median(data_by_charge[i])
                n = len(data_by_charge[i])
                ax.text(i, med + 0.2, f"{med:.2f}%\nn={n:,}",
                        ha="center", fontsize=7)
        ax.set_xlabel("Charge State")
        ax.set_ylabel("Relative Error (%)")
        ax.set_title("CCS Prediction Error by Charge State")
        ax.set_ylim(0, min(valid["rel_error_pct"].quantile(0.98), 20))
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)

        # --- Page 4: Residual vs observed mobility (bias plot) ---
        fig, ax = plt.subplots()
        hb = ax.hexbin(valid["mean_mobility_obs"], valid["residual"],
                       gridsize=60, cmap="RdBu_r", mincnt=1,
                       reduce_C_function=np.median)
        fig.colorbar(hb, ax=ax, label="Median residual")
        ax.axhline(0, color="black", linestyle="--", alpha=0.3)
        ax.set_xlabel("Observed 1/K0 (Vs/cm²)")
        ax.set_ylabel("Residual (predicted - observed)")
        ax.set_title("Prediction Bias vs Mobility Range")
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)

        # --- Page 5: Residual vs observed mobility (relative, %) ---
        fig, ax = plt.subplots()
        hb = ax.hexbin(valid["mean_mobility_obs"],
                       valid["residual"] / valid["mean_mobility_obs"] * 100,
                       gridsize=60, cmap="RdBu_r", mincnt=1,
                       reduce_C_function=np.median)
        fig.colorbar(hb, ax=ax, label="Median relative residual (%)")
        ax.axhline(0, color="black", linestyle="--", alpha=0.3)
        ax.set_xlabel("Observed 1/K0 (Vs/cm²)")
        ax.set_ylabel("Relative residual (%)")
        ax.set_title("Relative Prediction Bias vs Mobility Range")
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)

        # --- Page 6: Summary table ---
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.axis("off")

        rows = []
        # Prediction accuracy
        rel = valid["rel_error_pct"]
        rows.append(("Prediction Error", "ALL",
                      f"{len(valid):,}",
                      f"{rel.median():.3f}%",
                      f"{rel.mean():.3f}%",
                      f"{rel.quantile(0.25):.3f}%",
                      f"{rel.quantile(0.75):.3f}%"))

        for c in [2, 3, 4]:
            s = valid[valid["charge"] == c]["rel_error_pct"]
            if len(s) > 0:
                rows.append(("", f"  Charge {c}+",
                              f"{len(s):,}",
                              f"{s.median():.3f}%",
                              f"{s.mean():.3f}%",
                              f"{s.quantile(0.25):.3f}%",
                              f"{s.quantile(0.75):.3f}%"))

        # Cross-sample reproducibility
        if has_cv:
            cv = df["cv_pct_aligned"].dropna()
            rows.append(("Cross-Sample CV%", "ALL (aligned)",
                          f"{len(cv):,}",
                          f"{cv.median():.3f}%",
                          f"{cv.mean():.3f}%",
                          f"{cv.quantile(0.25):.3f}%",
                          f"{cv.quantile(0.75):.3f}%"))
        if has_cv_raw:
            cv_r = df["cv_pct"].dropna()
            rows.append(("Cross-Sample CV%", "ALL (raw)",
                          f"{len(cv_r):,}",
                          f"{cv_r.median():.3f}%",
                          f"{cv_r.mean():.3f}%",
                          f"{cv_r.quantile(0.25):.3f}%",
                          f"{cv_r.quantile(0.75):.3f}%"))

        # Ratio
        if has_cv:
            ratio = rel.median() / cv.median()
            rows.append(("", "", "", "", "", "", ""))
            rows.append(("Ratio", "pred_error / CV_aligned",
                          "", f"{ratio:.1f}x", "", "", ""))

        col_labels = ["Metric", "Stratum", "N", "Median", "Mean", "Q25", "Q75"]
        table = ax.table(cellText=rows, colLabels=col_labels, loc="center",
                         cellLoc="center")
        table.auto_set_font_size(False)
        table.set_fontsize(8)
        table.scale(1.0, 1.5)

        for j in range(len(col_labels)):
            table[0, j].set_facecolor("#1e293b")
            table[0, j].set_text_props(color="white", fontweight="bold")

        ax.set_title("CCS Prediction Accuracy vs Cross-Sample Reproducibility",
                     fontsize=13, fontweight="bold", pad=20)
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)


def main():
    parser = argparse.ArgumentParser(
        description="CCS prediction accuracy vs cross-sample reproducibility"
    )
    parser.add_argument("--store", required=True,
                        help="Path to precursor_store.parquet")
    parser.add_argument("--repro",
                        default="notebook/analysis/ccs_reproducibility/ccs_reproducibility.parquet",
                        help="Path to ccs_reproducibility.parquet (for CV%% data)")
    parser.add_argument("--output",
                        default="notebook/analysis/ccs_prediction_accuracy/",
                        help="Output directory")
    parser.add_argument("--min-files", type=int, default=3,
                        help="Minimum distinct raw files per peptide (default: 3)")
    args = parser.parse_args()

    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load per-peptide mobilities
    print("Loading per-peptide mobilities ...", flush=True)
    df = load_peptide_mobilities(args.store, args.repro, min_files=args.min_files)
    print(f"  {len(df):,} peptides with >= {args.min_files} files", flush=True)

    # Predict 1/K0
    print("\nPredicting 1/K0 with DeepPeptideIonMobilityApex ...", flush=True)
    df["predicted_mobility"] = predict_mobilities(df)

    # Use aligned mobility if available, else raw
    if "mean_mobility_aligned" in df.columns and df["mean_mobility_aligned"].notna().any():
        df["mean_mobility_obs"] = df["mean_mobility_aligned"].fillna(df["mean_mobility"])
        print("  Using aligned mobility as reference", flush=True)
    else:
        df["mean_mobility_obs"] = df["mean_mobility"]
        print("  Using raw mean mobility as reference", flush=True)

    # Compute raw errors (before alignment)
    valid_mask = df["predicted_mobility"].notna() & df["mean_mobility_obs"].notna()
    df["residual_raw"] = df["predicted_mobility"] - df["mean_mobility_obs"]

    valid_raw = df.loc[valid_mask]
    median_offset = valid_raw["residual_raw"].median()
    print(f"\n  Raw prediction offset (median residual): {median_offset:+.5f} Vs/cm²",
          flush=True)

    # Align predictions: subtract global median offset (analogous to per-file CCS alignment)
    df["predicted_aligned"] = df["predicted_mobility"] - median_offset
    df["residual"] = df["predicted_aligned"] - df["mean_mobility_obs"]
    df["abs_error"] = df["residual"].abs()
    df["rel_error_pct"] = (df["abs_error"] / df["mean_mobility_obs"]) * 100

    valid = df.loc[valid_mask]
    med_err = valid["rel_error_pct"].median()
    print(f"\n  After alignment (median offset removed):", flush=True)
    print(f"  Median relative error: {med_err:.3f}%", flush=True)
    print(f"  Mean relative error:  {valid['rel_error_pct'].mean():.3f}%", flush=True)
    print(f"  Q25/Q75:              {valid['rel_error_pct'].quantile(0.25):.3f}% / "
          f"{valid['rel_error_pct'].quantile(0.75):.3f}%", flush=True)
    print(f"  Median abs error:     {valid['abs_error'].median():.5f} Vs/cm²", flush=True)

    r2 = np.corrcoef(valid["mean_mobility_obs"], valid["predicted_aligned"])[0, 1] ** 2
    print(f"  R²:                   {r2:.6f}", flush=True)

    if "cv_pct_aligned" in df.columns:
        cv_med = df["cv_pct_aligned"].dropna().median()
        ratio = med_err / cv_med
        print(f"\n  Cross-sample CV% (aligned): {cv_med:.3f}%", flush=True)
        print(f"  Ratio (pred_error / CV):    {ratio:.1f}x", flush=True)

    # Save
    out_cols = ["sequence_normalized", "charge", "n_files", "n_obs",
                "mean_mobility", "mean_mz", "mean_mobility_obs",
                "predicted_mobility", "predicted_aligned",
                "residual", "abs_error", "rel_error_pct"]
    if "cv_pct" in df.columns:
        out_cols.append("cv_pct")
    if "cv_pct_aligned" in df.columns:
        out_cols.append("cv_pct_aligned")
    out_cols = [c for c in out_cols if c in df.columns]

    parquet_path = output_dir / "ccs_prediction_accuracy.parquet"
    df[out_cols].to_parquet(parquet_path, index=False)
    print(f"\nSaved {parquet_path} ({len(df):,} rows)", flush=True)

    # Report
    pdf_path = output_dir / "ccs_prediction_accuracy_report.pdf"
    print(f"Generating report: {pdf_path} ...", flush=True)
    generate_report(df, pdf_path)
    print(f"Done. Report saved to {pdf_path}", flush=True)


if __name__ == "__main__":
    main()
