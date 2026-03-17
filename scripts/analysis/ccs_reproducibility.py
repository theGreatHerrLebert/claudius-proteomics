#!/usr/bin/env python3
"""
CCS (1/K0) Reproducibility Analysis

Within-dataset 1/K0 reproducibility for timsTOF data. For each (sequence, charge),
computes CV% of 1/K0 across raw files. Stratifies by quality tier and engine consensus.

Includes per-file CCS alignment: uses reference peptides (seen in many files) to
compute a per-file 1/K0 offset, then reports before/after CV% improvement.

Output:
    - ccs_reproducibility.parquet: one row per (sequence, charge) peptide precursor
    - ccs_reproducibility_report.pdf: 8-page diagnostic report (with alignment pages)
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


# ---------------------------------------------------------------------------
# Data loading & aggregation
# ---------------------------------------------------------------------------

def load_precursor_data(store_path: str, columns: list[str] | None = None) -> pd.DataFrame:
    """Load required columns from precursor_store.parquet."""
    needed = columns or [
        "sequence_normalized", "charge", "mobility",
        "raw_file", "n_engines", "is_high_quality",
        "ms1_im_r2", "isotope_cosim",
    ]
    pf = pq.ParquetFile(store_path)
    # Read only needed columns across all row groups
    table = pf.read(columns=needed)
    df = table.to_pandas()
    return df


def compute_reproducibility(df: pd.DataFrame, min_files: int = 3) -> pd.DataFrame:
    """Group by (sequence_normalized, charge), compute 1/K0 CV% stats.

    Returns one row per peptide precursor.
    """
    # Filter: identified (n_engines >= 1) and non-null mobility
    mask = (df["n_engines"] >= 1) & df["mobility"].notna() & (df["mobility"] > 0)
    df = df.loc[mask].copy()

    # Unique raw files per peptide
    grouped = df.groupby(["sequence_normalized", "charge"])

    agg = grouped.agg(
        n_obs=("mobility", "size"),
        n_files=("raw_file", "nunique"),
        mean_mobility=("mobility", "mean"),
        std_mobility=("mobility", "std"),
        median_im_r2=("ms1_im_r2", "median"),
        median_isotope_cosim=("isotope_cosim", "median"),
        # consensus info — take max n_engines seen for this peptide
        max_n_engines=("n_engines", "max"),
        # fraction of observations that are high quality
        frac_high_quality=("is_high_quality", "mean"),
    ).reset_index()

    # Require peptide seen in >= min_files distinct raw files
    agg = agg[agg["n_files"] >= min_files].copy()

    # CV%
    agg["cv_pct"] = (agg["std_mobility"] / agg["mean_mobility"]) * 100.0
    # Replace NaN CV (from single observation std=NaN) with 0
    agg["cv_pct"] = agg["cv_pct"].fillna(0.0)

    # Quality tier (majority vote)
    agg["quality_tier"] = np.where(agg["frac_high_quality"] >= 0.5, "High Quality", "Standard")

    # Consensus tier
    agg["consensus_tier"] = pd.cut(
        agg["max_n_engines"],
        bins=[0, 1, 2, 3],
        labels=["1/3 engines", "2/3 engines", "3/3 engines"],
        include_lowest=True,
    )

    return agg


# ---------------------------------------------------------------------------
# Per-file CCS alignment
# ---------------------------------------------------------------------------

def compute_file_offsets(
    df: pd.DataFrame,
    min_ref_files: int = 30,
) -> tuple[pd.DataFrame, dict[str, float]]:
    """Compute per-file 1/K0 alignment offsets from reference peptides.

    Reference peptides are (sequence_normalized, charge) pairs seen in >= min_ref_files
    distinct raw files. For each, compute global mean 1/K0. Then for each raw file,
    compute median residual (observation - global mean) → per-file offset.

    Args:
        df: Raw precursor observations (filtered to n_engines >= 1, mobility > 0)
        min_ref_files: Minimum number of files a peptide must appear in to be
                       used as a reference (default: 30)

    Returns:
        (df_corrected, offsets_dict)
        df_corrected: copy of df with 'mobility_aligned' column added
        offsets_dict: {raw_file: offset_value} mapping
    """
    mask = (df["n_engines"] >= 1) & df["mobility"].notna() & (df["mobility"] > 0)
    df_id = df.loc[mask].copy()

    # Find reference peptides (seen in many files)
    peptide_files = df_id.groupby(["sequence_normalized", "charge"])["raw_file"].nunique()
    ref_peptides = set(peptide_files[peptide_files >= min_ref_files].index)
    print(f"  {len(ref_peptides):,} reference peptides (seen in >= {min_ref_files} files)",
          flush=True)

    if len(ref_peptides) == 0:
        print(f"  WARNING: No reference peptides found, skipping alignment", flush=True)
        df["mobility_aligned"] = df["mobility"]
        return df, {}

    # Compute global mean 1/K0 per reference peptide
    ref_mask = df_id.set_index(["sequence_normalized", "charge"]).index.isin(ref_peptides)
    df_ref = df_id.loc[ref_mask].copy()

    global_means = df_ref.groupby(["sequence_normalized", "charge"])["mobility"].mean()
    global_means.name = "global_mean_mobility"

    # Merge global mean back to reference observations
    df_ref = df_ref.join(global_means, on=["sequence_normalized", "charge"])
    df_ref["residual"] = df_ref["mobility"] - df_ref["global_mean_mobility"]

    # Per-file median residual = alignment offset
    file_offsets = df_ref.groupby("raw_file")["residual"].median()
    offsets_dict = file_offsets.to_dict()

    print(f"  Per-file offsets: min={file_offsets.min():.6f}, "
          f"max={file_offsets.max():.6f}, "
          f"range={file_offsets.max() - file_offsets.min():.6f}", flush=True)

    # Apply correction to ALL observations
    df = df.copy()
    df["file_offset"] = df["raw_file"].map(offsets_dict).fillna(0.0)
    df["mobility_aligned"] = df["mobility"] - df["file_offset"]
    df = df.drop(columns=["file_offset"])

    return df, offsets_dict


def compute_reproducibility_aligned(
    df: pd.DataFrame, min_files: int = 3,
) -> pd.DataFrame:
    """Same as compute_reproducibility but on mobility_aligned column."""
    mask = (df["n_engines"] >= 1) & df["mobility_aligned"].notna() & (df["mobility_aligned"] > 0)
    df = df.loc[mask].copy()

    grouped = df.groupby(["sequence_normalized", "charge"])

    agg = grouped.agg(
        n_obs=("mobility_aligned", "size"),
        n_files=("raw_file", "nunique"),
        mean_mobility_aligned=("mobility_aligned", "mean"),
        std_mobility_aligned=("mobility_aligned", "std"),
        max_n_engines=("n_engines", "max"),
        frac_high_quality=("is_high_quality", "mean"),
    ).reset_index()

    agg = agg[agg["n_files"] >= min_files].copy()
    agg["cv_pct_aligned"] = (agg["std_mobility_aligned"] / agg["mean_mobility_aligned"]) * 100.0
    agg["cv_pct_aligned"] = agg["cv_pct_aligned"].fillna(0.0)

    agg["quality_tier"] = np.where(agg["frac_high_quality"] >= 0.5, "High Quality", "Standard")
    agg["consensus_tier"] = pd.cut(
        agg["max_n_engines"],
        bins=[0, 1, 2, 3],
        labels=["1/3 engines", "2/3 engines", "3/3 engines"],
        include_lowest=True,
    )

    return agg


# ---------------------------------------------------------------------------
# Report generation (8-page PDF)
# ---------------------------------------------------------------------------

def _set_style():
    plt.rcParams.update({
        "figure.dpi": 150,
        "font.size": 9,
        "axes.titlesize": 11,
        "axes.labelsize": 10,
        "figure.figsize": (8, 6),
    })


def generate_report(
    df: pd.DataFrame,
    output_path: Path,
    df_aligned: pd.DataFrame | None = None,
    file_offsets: dict | None = None,
):
    """Generate CCS reproducibility PDF report (6 base pages + 2 alignment pages)."""
    _set_style()

    with PdfPages(str(output_path)) as pdf:
        # --- Page 1: CV% histogram by quality tier ---
        fig, ax = plt.subplots()
        for tier, color in [("High Quality", "#2563eb"), ("Standard", "#dc2626")]:
            subset = df[df["quality_tier"] == tier]["cv_pct"]
            if len(subset) == 0:
                continue
            # Clip for display
            clipped = subset.clip(upper=5.0)
            ax.hist(clipped, bins=100, alpha=0.6, label=f"{tier} (n={len(subset):,})",
                    color=color, density=True)
        ax.set_xlabel("CV% of 1/K0")
        ax.set_ylabel("Density")
        ax.set_title("1/K0 Reproducibility: CV% Distribution by Quality Tier")
        ax.legend()
        ax.set_xlim(0, 5)
        med_all = df["cv_pct"].median()
        ax.axvline(med_all, color="black", linestyle="--", alpha=0.5,
                   label=f"Overall median = {med_all:.3f}%")
        ax.legend()
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)

        # --- Page 2: CV% by charge state (violin/box) ---
        fig, ax = plt.subplots()
        charge_vals = sorted(df["charge"].unique())
        # Limit to charges 2-5 for clarity
        charge_vals = [c for c in charge_vals if 2 <= c <= 5]
        data_by_charge = [df[df["charge"] == c]["cv_pct"].clip(upper=5.0).values for c in charge_vals]
        if data_by_charge:
            parts = ax.violinplot(data_by_charge, positions=range(len(charge_vals)),
                                  showmedians=True, showextrema=False)
            for pc in parts["bodies"]:
                pc.set_facecolor("#6366f1")
                pc.set_alpha(0.6)
            ax.set_xticks(range(len(charge_vals)))
            ax.set_xticklabels([f"{c}+" for c in charge_vals])
            # Annotate medians
            for i, c in enumerate(charge_vals):
                med = np.median(data_by_charge[i])
                n = len(data_by_charge[i])
                ax.text(i, med + 0.1, f"{med:.3f}%\nn={n:,}", ha="center", fontsize=7)
        ax.set_xlabel("Charge State")
        ax.set_ylabel("CV% of 1/K0")
        ax.set_title("1/K0 CV% by Precursor Charge State")
        ax.set_ylim(0, 5)
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)

        # --- Page 3: CV% vs n_observations ---
        fig, ax = plt.subplots()
        # Bin n_obs for clarity
        df_plot = df[["n_obs", "cv_pct"]].copy()
        df_plot["cv_clipped"] = df_plot["cv_pct"].clip(upper=5.0)
        max_obs = df_plot["n_obs"].quantile(0.99)
        df_plot = df_plot[df_plot["n_obs"] <= max_obs]
        # Hexbin for density
        hb = ax.hexbin(df_plot["n_obs"], df_plot["cv_clipped"],
                       gridsize=40, cmap="YlOrRd", mincnt=1,
                       reduce_C_function=np.median)
        cb = fig.colorbar(hb, ax=ax)
        cb.set_label("Median CV%")
        ax.set_xlabel("Number of Observations")
        ax.set_ylabel("CV% of 1/K0")
        ax.set_title("Reproducibility vs Observation Count")
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)

        # --- Page 4: CV% vs mean_mobility ---
        fig, ax = plt.subplots()
        hb = ax.hexbin(df["mean_mobility"], df["cv_pct"].clip(upper=3.0),
                       gridsize=50, cmap="viridis", mincnt=1,
                       reduce_C_function=np.median)
        cb = fig.colorbar(hb, ax=ax)
        cb.set_label("Median CV%")
        ax.set_xlabel("Mean 1/K0 (V s cm$^{-2}$)")
        ax.set_ylabel("CV% of 1/K0")
        ax.set_title("Mobility-Dependent Reproducibility")
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)

        # --- Page 5: CV% by engine consensus ---
        fig, ax = plt.subplots()
        consensus_tiers = ["1/3 engines", "2/3 engines", "3/3 engines"]
        data_by_consensus = []
        labels = []
        for tier in consensus_tiers:
            subset = df[df["consensus_tier"] == tier]["cv_pct"].clip(upper=5.0).values
            if len(subset) > 0:
                data_by_consensus.append(subset)
                labels.append(tier)
        if data_by_consensus:
            parts = ax.violinplot(data_by_consensus, positions=range(len(labels)),
                                  showmedians=True, showextrema=False)
            for pc in parts["bodies"]:
                pc.set_facecolor("#10b981")
                pc.set_alpha(0.6)
            ax.set_xticks(range(len(labels)))
            ax.set_xticklabels(labels)
            for i, lbl in enumerate(labels):
                med = np.median(data_by_consensus[i])
                n = len(data_by_consensus[i])
                ax.text(i, med + 0.1, f"{med:.3f}%\nn={n:,}", ha="center", fontsize=7)
        ax.set_xlabel("Engine Consensus")
        ax.set_ylabel("CV% of 1/K0")
        ax.set_title("1/K0 CV% by Search Engine Agreement")
        ax.set_ylim(0, 5)
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)

        # --- Page 6: Per-file alignment offsets ---
        if file_offsets and df_aligned is not None:
            fig, axes = plt.subplots(1, 2, figsize=(12, 5))

            # Left: bar chart of per-file offsets
            ax = axes[0]
            files_sorted = sorted(file_offsets.items(), key=lambda x: x[1])
            file_names = [f.replace(".d", "").split("/")[-1][-12:] for f, _ in files_sorted]
            offsets_vals = [v for _, v in files_sorted]
            colors = ["#dc2626" if v > 0 else "#2563eb" for v in offsets_vals]
            ax.barh(range(len(file_names)), offsets_vals, color=colors, alpha=0.7, height=0.8)
            ax.set_yticks(range(len(file_names)))
            ax.set_yticklabels(file_names, fontsize=5)
            ax.set_xlabel("1/K0 Offset (V s cm$^{-2}$)")
            ax.set_title("Per-File 1/K0 Alignment Offsets")
            ax.axvline(0, color="black", linestyle="-", alpha=0.3)

            # Right: before/after CV% histograms
            ax = axes[1]
            cv_before = df["cv_pct"].clip(upper=5.0)
            cv_after = df_aligned["cv_pct_aligned"].clip(upper=5.0)
            ax.hist(cv_before, bins=100, alpha=0.5, density=True,
                    label=f"Before (median={df['cv_pct'].median():.4f}%)", color="#dc2626")
            ax.hist(cv_after, bins=100, alpha=0.5, density=True,
                    label=f"After (median={df_aligned['cv_pct_aligned'].median():.4f}%)",
                    color="#2563eb")
            ax.set_xlabel("CV% of 1/K0")
            ax.set_ylabel("Density")
            ax.set_title("CV% Distribution: Before vs After Alignment")
            ax.set_xlim(0, 5)
            ax.legend(fontsize=8)

            fig.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)

            # --- Page 7: Aligned CV% by quality tier ---
            fig, axes = plt.subplots(1, 2, figsize=(12, 5))

            for ax_idx, (tier, color) in enumerate(
                [("High Quality", "#2563eb"), ("Standard", "#dc2626")]
            ):
                ax = axes[ax_idx]
                before = df[df["quality_tier"] == tier]["cv_pct"].clip(upper=5.0)
                after = df_aligned[df_aligned["quality_tier"] == tier]["cv_pct_aligned"].clip(upper=5.0)
                if len(before) > 0:
                    ax.hist(before, bins=80, alpha=0.5, density=True,
                            label=f"Before (med={before.median():.4f}%)", color="#dc2626")
                    ax.hist(after, bins=80, alpha=0.5, density=True,
                            label=f"After (med={after.median():.4f}%)", color="#2563eb")
                    improvement = (1 - after.median() / before.median()) * 100
                    ax.set_title(f"{tier}: {improvement:.1f}% improvement")
                    ax.set_xlabel("CV% of 1/K0")
                    ax.set_ylabel("Density")
                    ax.set_xlim(0, 5)
                    ax.legend(fontsize=8)

            fig.suptitle("Per-File Alignment Effect by Quality Tier", fontsize=12, fontweight="bold")
            fig.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)

        # --- Summary table page ---
        fig, ax = plt.subplots(figsize=(8, 7))
        ax.axis("off")

        rows = []
        has_aligned = df_aligned is not None and "cv_pct_aligned" in df_aligned.columns

        # Overall
        row = ["ALL", f"{len(df):,}", f"{df['cv_pct'].median():.4f}",
               f"{df['cv_pct'].mean():.4f}",
               f"{(df['cv_pct'] < 1.0).sum() / len(df) * 100:.1f}%"]
        if has_aligned:
            row.append(f"{df_aligned['cv_pct_aligned'].median():.4f}")
        rows.append(tuple(row))

        # By quality
        for tier in ["High Quality", "Standard"]:
            s = df[df["quality_tier"] == tier]
            if len(s) > 0:
                row = [f"  {tier}", f"{len(s):,}", f"{s['cv_pct'].median():.4f}",
                       f"{s['cv_pct'].mean():.4f}",
                       f"{(s['cv_pct'] < 1.0).sum() / len(s) * 100:.1f}%"]
                if has_aligned:
                    sa = df_aligned[df_aligned["quality_tier"] == tier]
                    row.append(f"{sa['cv_pct_aligned'].median():.4f}" if len(sa) > 0 else "—")
                rows.append(tuple(row))

        # By consensus
        consensus_tiers = ["1/3 engines", "2/3 engines", "3/3 engines"]
        for tier in consensus_tiers:
            s = df[df["consensus_tier"] == tier]
            if len(s) > 0:
                row = [f"  {tier}", f"{len(s):,}", f"{s['cv_pct'].median():.4f}",
                       f"{s['cv_pct'].mean():.4f}",
                       f"{(s['cv_pct'] < 1.0).sum() / len(s) * 100:.1f}%"]
                if has_aligned:
                    sa = df_aligned[df_aligned["consensus_tier"] == tier]
                    row.append(f"{sa['cv_pct_aligned'].median():.4f}" if len(sa) > 0 else "—")
                rows.append(tuple(row))

        # By charge
        for c in sorted(df["charge"].unique()):
            if c < 2 or c > 5:
                continue
            s = df[df["charge"] == c]
            if len(s) > 0:
                row = [f"  Charge {c}+", f"{len(s):,}", f"{s['cv_pct'].median():.4f}",
                       f"{s['cv_pct'].mean():.4f}",
                       f"{(s['cv_pct'] < 1.0).sum() / len(s) * 100:.1f}%"]
                if has_aligned:
                    sa = df_aligned[df_aligned["charge"] == c]
                    row.append(f"{sa['cv_pct_aligned'].median():.4f}" if len(sa) > 0 else "—")
                rows.append(tuple(row))

        col_labels = ["Stratum", "N Peptides", "Median CV%", "Mean CV%", "% < 1% CV"]
        if has_aligned:
            col_labels.append("Aligned CV%")
        table = ax.table(cellText=rows, colLabels=col_labels, loc="center",
                         cellLoc="center")
        table.auto_set_font_size(False)
        table.set_fontsize(9)
        table.scale(1.0, 1.4)

        # Style header
        for j in range(len(col_labels)):
            table[0, j].set_facecolor("#1e293b")
            table[0, j].set_text_props(color="white", fontweight="bold")

        ax.set_title("CCS Reproducibility Summary", fontsize=13, fontweight="bold", pad=20)
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="CCS (1/K0) reproducibility analysis across raw files"
    )
    parser.add_argument("--store", required=True,
                        help="Path to precursor_store.parquet")
    parser.add_argument("--output", default="notebook/analysis/ccs_reproducibility/",
                        help="Output directory")
    parser.add_argument("--min-files", type=int, default=3,
                        help="Minimum distinct raw files per peptide (default: 3)")
    args = parser.parse_args()

    store_path = Path(args.store)
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    if not store_path.exists():
        print(f"ERROR: Store not found: {store_path}")
        sys.exit(1)

    # Load data
    print(f"Loading precursor data from {store_path} ...")
    df = load_precursor_data(str(store_path))
    print(f"  Loaded {len(df):,} precursor observations")
    print(f"  {df['raw_file'].nunique()} raw files")
    print(f"  {df['n_engines'].ge(1).sum():,} identified (n_engines >= 1)")

    # Compute reproducibility (raw)
    print(f"\nComputing reproducibility (min_files={args.min_files}) ...")
    result = compute_reproducibility(df, min_files=args.min_files)
    print(f"  {len(result):,} peptide precursors with >= {args.min_files} raw files")
    print(f"  Median CV%: {result['cv_pct'].median():.4f}")
    print(f"  Mean CV%:   {result['cv_pct'].mean():.4f}")
    print(f"  < 1% CV:    {(result['cv_pct'] < 1.0).sum():,} "
          f"({(result['cv_pct'] < 1.0).sum() / len(result) * 100:.1f}%)")

    # Per-file CCS alignment
    print(f"\nComputing per-file 1/K0 alignment ...")
    min_ref = max(int(df["raw_file"].nunique() * 0.4), 3)
    df_corrected, file_offsets = compute_file_offsets(df, min_ref_files=min_ref)

    result_aligned = None
    if file_offsets:
        result_aligned = compute_reproducibility_aligned(
            df_corrected, min_files=args.min_files
        )
        print(f"  After alignment:")
        print(f"    Median CV%: {result_aligned['cv_pct_aligned'].median():.4f}")
        print(f"    Mean CV%:   {result_aligned['cv_pct_aligned'].mean():.4f}")
        improvement = (1 - result_aligned['cv_pct_aligned'].median() /
                       result['cv_pct'].median()) * 100
        print(f"    Improvement: {improvement:.1f}%")

    # Save parquet (include aligned CV% if available)
    if result_aligned is not None:
        result = result.merge(
            result_aligned[["sequence_normalized", "charge", "cv_pct_aligned",
                            "mean_mobility_aligned"]],
            on=["sequence_normalized", "charge"],
            how="left",
        )

    parquet_path = output_dir / "ccs_reproducibility.parquet"
    result.to_parquet(parquet_path, index=False)
    print(f"\nSaved {parquet_path} ({len(result):,} rows)")

    # Generate report
    pdf_path = output_dir / "ccs_reproducibility_report.pdf"
    print(f"Generating report: {pdf_path} ...")
    generate_report(result, pdf_path, df_aligned=result_aligned,
                    file_offsets=file_offsets)
    print(f"Done. Report saved to {pdf_path}")


if __name__ == "__main__":
    main()
