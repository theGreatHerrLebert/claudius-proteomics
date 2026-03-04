#!/usr/bin/env python3
"""
Spectral Reproducibility Analysis

Within-dataset MS2 spectral reproducibility for timsTOF data. For each
(sequence, charge) peptide seen in multiple raw files, computes pairwise
spectral angles between observed fragment spectra. Stratifies by quality
tier, engine consensus, and charge state — mirroring the CCS reproducibility
analysis but for MS2 spectral agreement.

Spectral angle (Prosit convention): SA = 1 - (2/π) * arccos(cosine_similarity).
Range [0, 1], higher = more reproducible.

Output:
    - spectral_reproducibility.parquet: one row per (sequence, charge)
    - spectral_reproducibility_report.pdf: 8-page diagnostic report
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

# For optional predicted spectrum comparison
RUSTIMS_PREDICTORS = Path("/home/administrator/Documents/promotion/rust/rustims/packages/imspy-predictors/src")


# ---------------------------------------------------------------------------
# Spectral angle computation
# ---------------------------------------------------------------------------

def build_spectrum_dict(
    ion_types: list[str],
    ion_numbers: list[int],
    ion_charges: list[int],
    intensities: list[float],
) -> dict[tuple, float]:
    """Build {(ion_type, ion_number, charge): sqrt(intensity)} dict.

    Applies sqrt transform (Prosit convention) to compress dynamic range
    before spectral angle computation.
    """
    spec = {}
    for itype, inum, ichg, iint in zip(ion_types, ion_numbers, ion_charges, intensities):
        if iint > 0:
            spec[(itype, inum, ichg)] = np.sqrt(float(iint))
    return spec


def spectral_angle(spec_a: dict, spec_b: dict) -> float:
    """Compute normalized spectral angle (Prosit convention) between two spectra.

    SA = 1 - (2/π) * arccos(cosine_similarity)
    Range [0, 1], where 1 = identical, 0 = orthogonal.

    Input spectra should already be sqrt-transformed. Vectors are built
    over the union of fragment positions (0 for missing).

    Returns SA in [0, 1], or NaN if either spectrum is empty.
    """
    keys = set(spec_a.keys()) | set(spec_b.keys())
    if not keys:
        return np.nan

    vec_a = np.array([spec_a.get(k, 0.0) for k in keys])
    vec_b = np.array([spec_b.get(k, 0.0) for k in keys])

    norm_a = np.linalg.norm(vec_a)
    norm_b = np.linalg.norm(vec_b)

    if norm_a == 0 or norm_b == 0:
        return 0.0

    cosim = np.dot(vec_a, vec_b) / (norm_a * norm_b)
    cosim = np.clip(cosim, -1.0, 1.0)
    return float(1.0 - (2.0 / np.pi) * np.arccos(cosim))


# ---------------------------------------------------------------------------
# Predicted spectrum comparison (optional, requires imspy-predictors)
# ---------------------------------------------------------------------------

_SUPPORTED_MODS = {"4", "35"}


def _clean_seq_for_predictor(seq: str) -> str | None:
    """Clean sequence for the local intensity predictor."""
    if not seq:
        return None
    cleaned = seq.replace("-", "").replace("[]", "")
    mods = re.findall(r"\[UNIMOD:(\d+)\]", cleaned)
    if any(m not in _SUPPORTED_MODS for m in mods):
        return None
    stripped = re.sub(r"\[UNIMOD:\d+\]", "", cleaned)
    if len(stripped) < 7 or len(stripped) > 30:
        return None
    return cleaned


def predict_and_compare(
    result_df: pd.DataFrame,
    consensus_spectra: dict,
    ce_offset: float = 0.0,
    batch_size: int = 512,
) -> pd.Series:
    """Predict spectra for consensus peptides and compute SA vs consensus.

    Args:
        result_df: DataFrame with sequence, charge, mean_collision_energy columns
        consensus_spectra: {(sequence, charge): {(ion_type, ion_num, ion_chg): sqrt_intensity}}
        ce_offset: CE calibration offset to apply
        batch_size: prediction batch size

    Returns:
        Series of SA values aligned with result_df index.
    """
    sys.path.insert(0, str(RUSTIMS_PREDICTORS))
    from imspy_predictors.intensity.predictors import DeepPeptideIntensityPredictor

    predictor = DeepPeptideIntensityPredictor(verbose=False)

    # Prepare inputs — only peptides with supported sequences
    sa_predicted = pd.Series(np.nan, index=result_df.index)

    seqs_clean = result_df["sequence"].apply(_clean_seq_for_predictor)
    valid_mask = seqs_clean.notna()
    valid_df = result_df.loc[valid_mask].copy()
    valid_seqs = seqs_clean.loc[valid_mask].tolist()

    if not valid_seqs:
        return sa_predicted

    charges = np.array(valid_df["charge"].tolist(), dtype=np.int64)
    ces = np.array(valid_df["mean_collision_energy"].tolist(), dtype=np.float64) + ce_offset

    print(f"  Predicting {len(valid_seqs):,} spectra with local model ...", flush=True)
    raw_preds = predictor.predict_intensities(
        sequences=valid_seqs,
        charges=charges,
        collision_energies=ces,
        batch_size=batch_size,
        flatten=False,
    )

    # Convert predictions to spectrum dicts and compute SA vs consensus
    for idx_pos, (df_idx, row) in enumerate(valid_df.iterrows()):
        arr = raw_preds[idx_pos]
        pred_spec = {}
        for pos in range(arr.shape[0]):
            for ion_type_idx, ion_type in enumerate(["y", "b"]):
                for chg_idx in range(arr.shape[2]):
                    intensity = float(arr[pos, ion_type_idx, chg_idx])
                    if intensity < 0:
                        continue
                    # sqrt-transform predicted to match observed consensus
                    pred_spec[(ion_type, pos + 1, chg_idx + 1)] = np.sqrt(max(intensity, 0.0))

        key = (row["sequence"], row["charge"])
        if key in consensus_spectra:
            sa_val = spectral_angle(pred_spec, consensus_spectra[key])
            sa_predicted.at[df_idx] = sa_val

    return sa_predicted


# ---------------------------------------------------------------------------
# Data loading and pairwise comparison
# ---------------------------------------------------------------------------

def load_fragment_data(path: str) -> pd.DataFrame:
    """Load fragment_peaks.parquet with needed columns."""
    needed = [
        "sequence", "charge", "raw_file", "n_engines", "collision_energy",
        "ion_type", "ion_number", "ion_charge", "intensity_mono",
        "n_matched", "coverage_b", "coverage_y", "intensity_explained",
    ]
    pf = pq.ParquetFile(path)
    return pf.read(columns=needed).to_pandas()


def compute_spectral_reproducibility(
    df: pd.DataFrame,
    min_files: int = 3,
    max_pairs_per_peptide: int = 20,
    rng_seed: int = 42,
) -> pd.DataFrame:
    """For each peptide precursor seen in multiple files, compute pairwise SA.

    Strategy per peptide group:
      - Pick one PSM per raw file (highest n_matched)
      - If > max_pairs_per_peptide files, subsample
      - Compute all pairwise spectral angles
      - Report median, mean, min, max SA

    Returns one row per (sequence, charge) with SA statistics.
    """
    rng = np.random.RandomState(rng_seed)

    # Filter to identified PSMs with decent fragment coverage
    mask = (df["n_engines"] >= 1) & (df["n_matched"] >= 3)
    df = df.loc[mask].copy()

    # Pre-build spectrum dicts for all PSMs (vectorized prep)
    print("  Building spectrum dicts ...", flush=True)
    spectra = []
    for _, row in df.iterrows():
        spectra.append(build_spectrum_dict(
            row["ion_type"], row["ion_number"],
            row["ion_charge"], row["intensity_mono"],
        ))
    df["_spec"] = spectra

    # Group by peptide identity
    grouped = df.groupby(["sequence", "charge"])

    results = []
    consensus_spectra_dict = {}
    n_skipped = 0

    for (seq, chg), group in grouped:
        # One PSM per raw file (best coverage)
        best_per_file = group.sort_values("n_matched", ascending=False).drop_duplicates(
            subset=["raw_file"], keep="first"
        )
        n_files = len(best_per_file)
        if n_files < min_files:
            n_skipped += 1
            continue

        # Subsample if too many files
        if n_files > max_pairs_per_peptide:
            best_per_file = best_per_file.sample(
                n=max_pairs_per_peptide, random_state=rng
            )
            n_files = max_pairs_per_peptide

        # Compute all pairwise spectral angles
        specs = best_per_file["_spec"].tolist()
        angles = []
        for i in range(n_files):
            for j in range(i + 1, n_files):
                sa = spectral_angle(specs[i], specs[j])
                if not np.isnan(sa):
                    angles.append(sa)

        if not angles:
            continue

        # Build consensus spectrum (mean intensity per fragment across runs)
        all_keys = set()
        for s in specs:
            all_keys.update(s.keys())
        consensus = {}
        for k in all_keys:
            vals = [s[k] for s in specs if k in s]
            consensus[k] = float(np.mean(vals))

        # Store for optional predicted comparison
        consensus_spectra_dict[(seq, chg)] = consensus

        # Compute individual-vs-consensus SA
        consensus_angles = []
        for s in specs:
            sa_c = spectral_angle(s, consensus)
            if not np.isnan(sa_c):
                consensus_angles.append(sa_c)

        # Aggregate per-peptide stats
        max_engines = int(group["n_engines"].max())
        mean_ce = float(group["collision_energy"].mean())
        mean_coverage = float(
            (group["coverage_b"] + group["coverage_y"]).mean() / 2
        )

        results.append({
            "sequence": seq,
            "charge": chg,
            "n_files": n_files,
            "n_pairs": len(angles),
            "median_sa": float(np.median(angles)),
            "mean_sa": float(np.mean(angles)),
            "std_sa": float(np.std(angles)),
            "min_sa": float(np.min(angles)),
            "max_sa": float(np.max(angles)),
            "median_sa_consensus": float(np.median(consensus_angles)) if consensus_angles else np.nan,
            "mean_sa_consensus": float(np.mean(consensus_angles)) if consensus_angles else np.nan,
            "max_n_engines": max_engines,
            "mean_collision_energy": mean_ce,
            "mean_coverage": mean_coverage,
            "n_obs": len(group),
        })

    print(f"  {len(results):,} peptides with >= {min_files} files "
          f"(skipped {n_skipped:,})", flush=True)

    result_df = pd.DataFrame(results)

    # Quality and consensus tiers
    result_df["consensus_tier"] = pd.cut(
        result_df["max_n_engines"],
        bins=[0, 1, 2, 3],
        labels=["1/3 engines", "2/3 engines", "3/3 engines"],
        include_lowest=True,
    )

    return result_df, consensus_spectra_dict


# ---------------------------------------------------------------------------
# Report generation
# ---------------------------------------------------------------------------

def _set_style():
    plt.rcParams.update({
        "figure.dpi": 150,
        "font.size": 9,
        "axes.titlesize": 11,
        "axes.labelsize": 10,
        "figure.figsize": (8, 6),
    })


def generate_report(df: pd.DataFrame, output_path: Path):
    """Generate spectral reproducibility PDF report."""
    _set_style()

    with PdfPages(str(output_path)) as pdf:
        # --- Page 1: SA histogram by consensus tier ---
        fig, ax = plt.subplots()
        colors = {"3/3 engines": "#2563eb", "2/3 engines": "#f59e0b", "1/3 engines": "#dc2626"}
        for tier in ["3/3 engines", "2/3 engines", "1/3 engines"]:
            subset = df[df["consensus_tier"] == tier]["median_sa"]
            if len(subset) == 0:
                continue
            ax.hist(subset, bins=100, alpha=0.5, density=True,
                    label=f"{tier} (n={len(subset):,})", color=colors[tier])
        ax.set_xlabel("Median Spectral Angle (SA)")
        ax.set_ylabel("Density")
        ax.set_title("MS2 Spectral Reproducibility: SA Distribution by Engine Consensus")
        med_all = df["median_sa"].median()
        ax.axvline(med_all, color="black", linestyle="--", alpha=0.5,
                   label=f"Overall median = {med_all:.3f}")
        ax.legend(fontsize=8)
        ax.set_xlim(0, 1)
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)

        # --- Page 2: SA by charge state (violin) ---
        fig, ax = plt.subplots()
        charge_vals = sorted(df["charge"].unique())
        charge_vals = [c for c in charge_vals if 2 <= c <= 5]
        data_by_charge = [
            df[df["charge"] == c]["median_sa"].values
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
                ax.text(i, med + 0.02, f"{med:.3f}\nn={n:,}",
                        ha="center", fontsize=7)
        ax.set_xlabel("Charge State")
        ax.set_ylabel("Median Spectral Angle (SA)")
        ax.set_title("Spectral Angle by Precursor Charge State")
        ax.set_ylim(0, 1)
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)

        # --- Page 3: SA vs number of raw files ---
        fig, ax = plt.subplots()
        df_plot = df[["n_files", "median_sa"]].copy()
        max_files = df_plot["n_files"].quantile(0.99)
        df_plot = df_plot[df_plot["n_files"] <= max_files]
        hb = ax.hexbin(df_plot["n_files"], df_plot["median_sa"],
                       gridsize=40, cmap="YlGnBu", mincnt=1,
                       reduce_C_function=np.median)
        cb = fig.colorbar(hb, ax=ax)
        cb.set_label("Median SA")
        ax.set_xlabel("Number of Raw Files")
        ax.set_ylabel("Median Spectral Angle (SA)")
        ax.set_title("Spectral Reproducibility vs File Count")
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)

        # --- Page 4: SA vs mean collision energy ---
        fig, ax = plt.subplots()
        hb = ax.hexbin(df["mean_collision_energy"],
                       df["median_sa"],
                       gridsize=50, cmap="viridis", mincnt=1,
                       reduce_C_function=np.median)
        cb = fig.colorbar(hb, ax=ax)
        cb.set_label("Median SA")
        ax.set_xlabel("Mean Collision Energy (eV)")
        ax.set_ylabel("Median Spectral Angle (SA)")
        ax.set_title("CE-Dependent Spectral Reproducibility")
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)

        # --- Page 5: SA by engine consensus (violin) ---
        fig, ax = plt.subplots()
        consensus_tiers = ["1/3 engines", "2/3 engines", "3/3 engines"]
        data_by_consensus = []
        labels = []
        for tier in consensus_tiers:
            subset = df[df["consensus_tier"] == tier]["median_sa"].values
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
                ax.text(i, med + 0.02, f"{med:.3f}\nn={n:,}",
                        ha="center", fontsize=7)
        ax.set_xlabel("Engine Consensus")
        ax.set_ylabel("Median Spectral Angle (SA)")
        ax.set_title("Spectral Angle by Search Engine Agreement")
        ax.set_ylim(0, 1)
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)

        # --- Page 6: SA vs fragment coverage ---
        fig, ax = plt.subplots()
        hb = ax.hexbin(df["mean_coverage"],
                       df["median_sa"],
                       gridsize=50, cmap="magma", mincnt=1,
                       reduce_C_function=np.median)
        cb = fig.colorbar(hb, ax=ax)
        cb.set_label("Median SA")
        ax.set_xlabel("Mean Fragment Coverage (avg of b + y)")
        ax.set_ylabel("Median Spectral Angle (SA)")
        ax.set_title("Coverage-Dependent Spectral Reproducibility")
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)

        # --- Page 7: SA within-peptide variability (std_sa vs median_sa) ---
        fig, ax = plt.subplots()
        mask = df["n_pairs"] >= 3  # Need enough pairs for meaningful std
        df_var = df.loc[mask]
        if len(df_var) > 0:
            hb = ax.hexbin(df_var["median_sa"],
                           df_var["std_sa"],
                           gridsize=50, cmap="plasma", mincnt=1)
            cb = fig.colorbar(hb, ax=ax)
            cb.set_label("Count")
            ax.set_xlabel("Median Spectral Angle (SA)")
            ax.set_ylabel("Std of Spectral Angle")
            ax.set_title("Within-Peptide SA Variability")
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)

        # --- Page 8: Consensus comparison ---
        has_consensus = "median_sa_consensus" in df.columns and df["median_sa_consensus"].notna().any()
        if has_consensus:
            fig, axes = plt.subplots(1, 2, figsize=(12, 5))

            sa_pair = df["median_sa"].dropna()
            sa_cons = df["median_sa_consensus"].dropna()

            # Left: overlay histograms
            ax = axes[0]
            ax.hist(sa_pair, bins=80, alpha=0.5, density=True,
                    label=f"Pairwise (med={sa_pair.median():.3f})", color="#dc2626")
            ax.hist(sa_cons, bins=80, alpha=0.5, density=True,
                    label=f"vs Consensus (med={sa_cons.median():.3f})", color="#2563eb")
            ax.set_xlabel("Spectral Angle (SA)")
            ax.set_ylabel("Density")
            ax.set_title("Pairwise vs Consensus SA")
            ax.set_xlim(0, 1)
            ax.legend(fontsize=8)

            # Right: scatter pairwise vs consensus
            ax = axes[1]
            both = df[["median_sa", "median_sa_consensus"]].dropna()
            hb = ax.hexbin(both["median_sa"], both["median_sa_consensus"],
                           gridsize=50, cmap="YlGnBu", mincnt=1)
            fig.colorbar(hb, ax=ax, label="Count")
            ax.plot([0, 1], [0, 1], "k--", alpha=0.3)
            ax.set_xlabel("Pairwise SA (median)")
            ax.set_ylabel("Consensus SA (median)")
            ax.set_title("Pairwise vs Consensus Agreement")
            ax.set_xlim(0, 1)
            ax.set_ylim(0, 1)

            fig.suptitle("Consensus Spectrum Comparison", fontsize=12, fontweight="bold")
            fig.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)

        # --- Page: Predicted vs Consensus (if available) ---
        has_predicted = "sa_predicted" in df.columns and df["sa_predicted"].notna().any()
        if has_predicted:
            fig, axes = plt.subplots(1, 2, figsize=(12, 5))

            sa_pred = df["sa_predicted"].dropna()

            ax = axes[0]
            if has_consensus:
                sa_cons = df["median_sa_consensus"].dropna()
                ax.hist(sa_cons, bins=80, alpha=0.5, density=True,
                        label=f"Obs vs Consensus (med={sa_cons.median():.3f})",
                        color="#2563eb")
            ax.hist(sa_pred, bins=80, alpha=0.5, density=True,
                    label=f"Pred vs Consensus (med={sa_pred.median():.3f})",
                    color="#10b981")
            ax.set_xlabel("Spectral Angle (SA)")
            ax.set_ylabel("Density")
            ax.set_title("Predicted vs Consensus Spectra")
            ax.set_xlim(0, 1)
            ax.legend(fontsize=8)

            # By charge
            ax = axes[1]
            charge_vals = sorted(df["charge"].unique())
            charge_vals = [c for c in charge_vals if 2 <= c <= 4]
            data_pred = [df[df["charge"] == c]["sa_predicted"].dropna().values
                         for c in charge_vals]
            if any(len(d) > 0 for d in data_pred):
                parts = ax.violinplot(
                    [d for d in data_pred if len(d) > 0],
                    positions=range(sum(1 for d in data_pred if len(d) > 0)),
                    showmedians=True, showextrema=False,
                )
                for pc in parts["bodies"]:
                    pc.set_facecolor("#10b981")
                    pc.set_alpha(0.6)
                valid_labels = [f"{c}+" for c, d in zip(charge_vals, data_pred) if len(d) > 0]
                ax.set_xticks(range(len(valid_labels)))
                ax.set_xticklabels(valid_labels)
                for i, (c, d) in enumerate(
                    [(c, d) for c, d in zip(charge_vals, data_pred) if len(d) > 0]
                ):
                    ax.text(i, np.median(d) + 0.02, f"{np.median(d):.3f}\nn={len(d):,}",
                            ha="center", fontsize=7)
            ax.set_xlabel("Charge State")
            ax.set_ylabel("SA (predicted vs consensus)")
            ax.set_title("Prediction Quality by Charge")
            ax.set_ylim(0, 1)

            fig.suptitle("Model Prediction vs Observed Consensus",
                         fontsize=12, fontweight="bold")
            fig.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)

        # --- Summary table ---
        fig, ax = plt.subplots(figsize=(12, 7))
        ax.axis("off")

        rows = []
        # Overall
        row = ["ALL", f"{len(df):,}",
               f"{df['median_sa'].median():.3f}"]
        if has_consensus:
            row.append(f"{df['median_sa_consensus'].median():.3f}")
        if has_predicted:
            sa_p = df["sa_predicted"].dropna()
            row.append(f"{sa_p.median():.3f}" if len(sa_p) > 0 else "—")
        rows.append(tuple(row))

        # By consensus tier
        consensus_tiers = ["1/3 engines", "2/3 engines", "3/3 engines"]
        for tier in consensus_tiers:
            s = df[df["consensus_tier"] == tier]
            if len(s) > 0:
                row = [f"  {tier}", f"{len(s):,}",
                       f"{s['median_sa'].median():.3f}"]
                if has_consensus:
                    row.append(f"{s['median_sa_consensus'].median():.3f}")
                if has_predicted:
                    sp = s["sa_predicted"].dropna()
                    row.append(f"{sp.median():.3f}" if len(sp) > 0 else "—")
                rows.append(tuple(row))

        # By charge
        for c in sorted(df["charge"].unique()):
            if c < 2 or c > 5:
                continue
            s = df[df["charge"] == c]
            if len(s) > 0:
                row = [f"  Charge {c}+", f"{len(s):,}",
                       f"{s['median_sa'].median():.3f}"]
                if has_consensus:
                    row.append(f"{s['median_sa_consensus'].median():.3f}")
                if has_predicted:
                    sp = s["sa_predicted"].dropna()
                    row.append(f"{sp.median():.3f}" if len(sp) > 0 else "—")
                rows.append(tuple(row))

        col_labels = ["Stratum", "N Peptides", "Pairwise SA"]
        if has_consensus:
            col_labels.append("vs Consensus")
        if has_predicted:
            col_labels.append("vs Predicted")
        table = ax.table(cellText=rows, colLabels=col_labels, loc="center",
                         cellLoc="center")
        table.auto_set_font_size(False)
        table.set_fontsize(9)
        table.scale(1.0, 1.4)

        for j in range(len(col_labels)):
            table[0, j].set_facecolor("#1e293b")
            table[0, j].set_text_props(color="white", fontweight="bold")

        ax.set_title("Spectral Reproducibility Summary", fontsize=13,
                     fontweight="bold", pad=20)
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="MS2 spectral reproducibility (spectral angle) across raw files"
    )
    parser.add_argument("--fragments", required=True,
                        help="Path to fragment_peaks.parquet")
    parser.add_argument("--output",
                        default="notebook/analysis/spectral_reproducibility/",
                        help="Output directory")
    parser.add_argument("--min-files", type=int, default=3,
                        help="Minimum distinct raw files per peptide (default: 3)")
    parser.add_argument("--max-pairs", type=int, default=20,
                        help="Max files per peptide for pairwise comparison (default: 20)")
    parser.add_argument("--predict", action="store_true",
                        help="Also compare consensus spectra against predicted (requires imspy-predictors)")
    parser.add_argument("--ce-offset", type=float, default=0.0,
                        help="CE calibration offset in eV (default: 0)")
    args = parser.parse_args()

    fragments_path = Path(args.fragments)
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    if not fragments_path.exists():
        print(f"ERROR: Fragment peaks not found: {fragments_path}")
        sys.exit(1)

    # Load data
    print(f"Loading fragment peaks from {fragments_path} ...", flush=True)
    df = load_fragment_data(str(fragments_path))
    print(f"  Loaded {len(df):,} PSMs", flush=True)
    print(f"  {df['raw_file'].nunique()} raw files", flush=True)

    # Compute spectral reproducibility
    print(f"\nComputing spectral reproducibility (min_files={args.min_files}) ...",
          flush=True)
    result, consensus_spectra = compute_spectral_reproducibility(
        df, min_files=args.min_files, max_pairs_per_peptide=args.max_pairs
    )
    med = result['median_sa'].median()
    print(f"  Median pairwise SA:  {med:.3f}", flush=True)
    print(f"  Mean pairwise SA:    {result['median_sa'].mean():.3f}", flush=True)
    print(f"  Q25/Q75:             {result['median_sa'].quantile(0.25):.3f} / "
          f"{result['median_sa'].quantile(0.75):.3f}", flush=True)

    med_cons = result["median_sa_consensus"].median()
    print(f"  Median vs consensus: {med_cons:.3f}", flush=True)

    # Predicted comparison (optional)
    if args.predict:
        print(f"\nComputing predicted vs consensus (CE offset: {args.ce_offset:+.0f} eV) ...",
              flush=True)
        sa_predicted = predict_and_compare(
            result, consensus_spectra, ce_offset=args.ce_offset
        )
        result["sa_predicted"] = sa_predicted
        valid = sa_predicted.dropna()
        print(f"  {len(valid):,} / {len(result):,} peptides with predictions", flush=True)
        print(f"  Median SA (pred vs consensus): {valid.median():.3f}", flush=True)

    # Save parquet
    parquet_path = output_dir / "spectral_reproducibility.parquet"
    result.to_parquet(parquet_path, index=False)
    print(f"\nSaved {parquet_path} ({len(result):,} rows)", flush=True)

    # Generate report
    pdf_path = output_dir / "spectral_reproducibility_report.pdf"
    print(f"Generating report: {pdf_path} ...", flush=True)
    generate_report(result, pdf_path)
    print(f"Done. Report saved to {pdf_path}", flush=True)


if __name__ == "__main__":
    main()
