#!/usr/bin/env python3
"""
Collision Energy Calibration via Spectral Prediction

For each dataset, predicts MS2 spectra with a local PyTorch model
(DeepPeptideIntensityPredictor / PROSPECT fine-tuned transformer) and finds
the CE offset that maximizes cosine similarity between predicted and observed
spectra. This calibrates all datasets to the predictor's CE scale as a
universal reference.

Approach:
    1. Sample high-confidence PSMs with observed fragment spectra
    2. Grid-search CE offsets (default -10 to +10 eV)
    3. At each offset: predict spectra at CE + offset, compute cosine similarity
    4. Optimal offset = argmax(mean cosine similarity)

Output:
    - ce_calibration.json: per-dataset CE offset and similarity curve
    - ce_calibration_report.pdf: diagnostic plots with before/after cossim
"""

import os
import argparse
import json
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

# Rustims / imspy paths
# Local checkout of https://github.com/theGreatHerrLebert/rustims.
# Override with RUSTIMS_ROOT; defaults to a sibling of this repository.
RUSTIMS_ROOT = Path(
    os.environ.get("RUSTIMS_ROOT", Path(__file__).resolve().parents[2] / ".." / "rustims")
).expanduser().resolve()
RUSTIMS_PREDICTORS = RUSTIMS_ROOT / "packages" / "imspy-predictors" / "src"
sys.path.insert(0, str(RUSTIMS_PREDICTORS))

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)


# ---------------------------------------------------------------------------
# Local model prediction
# ---------------------------------------------------------------------------

def get_local_predictor():
    """Initialize local PyTorch intensity predictor (PROSPECT fine-tuned)."""
    from imspy_predictors.intensity.predictors import DeepPeptideIntensityPredictor
    return DeepPeptideIntensityPredictor(verbose=False)


def predict_spectra_local(
    predictor,
    sequences: list[str],
    charges: list[int],
    collision_energies: list[float],
    batch_size: int = 512,
) -> list[dict]:
    """Predict fragment spectra locally, returning per-PSM fragment dicts.

    Returns list of dicts, one per input PSM:
        {(ion_type, ion_number, charge): predicted_intensity, ...}

    The local model returns (29, 2, 3) arrays per PSM:
        arr[pos, 0, chg_idx] = y(pos+1) at charge (chg_idx+1)
        arr[pos, 1, chg_idx] = b(pos+1) at charge (chg_idx+1)

    We convert these to the same dict format as our fragment_peaks.parquet
    uses: (ion_type, ion_number, charge) keys.
    """
    charges_arr = np.array(charges, dtype=np.int64)
    ces_arr = np.array(collision_energies, dtype=np.float64)

    # Local model prediction — no network calls
    raw_predictions = predictor.predict_intensities(
        sequences=sequences,
        charges=charges_arr,
        collision_energies=ces_arr,
        batch_size=batch_size,
        flatten=False,  # Returns list of (29, 2, 3) arrays
    )

    # Convert (29, 2, 3) arrays to {(ion_type, ion_number, charge): intensity} dicts
    spectra = []
    for arr in raw_predictions:
        spec = {}
        for pos in range(arr.shape[0]):
            for ion_type_idx, ion_type in enumerate(["y", "b"]):
                for chg_idx in range(arr.shape[2]):
                    intensity = float(arr[pos, ion_type_idx, chg_idx])
                    # Skip masked/invalid positions (model uses -1.0)
                    if intensity < 0:
                        continue
                    ion_number = pos + 1   # 1-indexed
                    ion_charge = chg_idx + 1  # 1-indexed
                    spec[(ion_type, ion_number, ion_charge)] = intensity
        spectra.append(spec)

    return spectra


# ---------------------------------------------------------------------------
# Sequence cleaning for the local model
# ---------------------------------------------------------------------------

def clean_sequence_for_predictor(seq: str) -> str | None:
    """Clean and validate a sequence for the local intensity predictor.

    The local model (ProformaTokenizer) supports all UNIMOD modifications.
    Only constraint is sequence length 7-30 AA.

    Returns cleaned sequence or None if unsupported.
    """
    if not seq or not isinstance(seq, str):
        return None

    # Remove hyphens and empty terminal brackets from canonical format
    # []-PEPTIDEK-[] → PEPTIDEK
    # [UNIMOD:1]-PEPTIDEK-[] → [UNIMOD:1]PEPTIDEK
    cleaned = seq.replace("-", "").replace("[]", "")

    # Check sequence length (strip mods to count AAs)
    stripped = re.sub(r"\[UNIMOD:\d+\]", "", cleaned)
    if len(stripped) < 7 or len(stripped) > 30:
        return None

    return cleaned


# ---------------------------------------------------------------------------
# Cosine similarity between predicted and observed
# ---------------------------------------------------------------------------

def spectral_angle_fragments(
    predicted: dict,
    observed_ions: list[str],
    observed_numbers: list[int],
    observed_charges: list[int],
    observed_intensities: list[float],
) -> float:
    """Compute normalized spectral angle between predicted and observed spectra.

    SA = 1 - (2/pi) * arccos(cosine_similarity)  [Prosit convention]
    Range [0, 1], where 1 = identical, 0 = orthogonal.

    Vectors are constructed over the predicted fragment ion positions
    (model output), matching against observed intensities.

    Args:
        predicted: {(ion_type, ion_number, charge): intensity} from model
        observed_ions: list of ion types ('b' or 'y')
        observed_numbers: list of ion numbers
        observed_charges: list of ion charges
        observed_intensities: list of observed monoisotopic intensities

    Returns:
        Normalized spectral angle [0, 1], or NaN if not computable.
    """
    if not predicted:
        return np.nan

    # Build observed dict from list columns
    observed = {}
    for itype, inum, ichg, iint in zip(
        observed_ions, observed_numbers, observed_charges, observed_intensities
    ):
        if iint > 0:
            observed[(itype, inum, ichg)] = iint

    # Compute similarity over predicted fragment positions
    keys = list(predicted.keys())
    if not keys:
        return np.nan

    pred_vec = np.array([predicted[k] for k in keys])
    obs_vec = np.array([observed.get(k, 0.0) for k in keys])

    pred_norm = np.linalg.norm(pred_vec)
    obs_norm = np.linalg.norm(obs_vec)

    if pred_norm == 0 or obs_norm == 0:
        return 0.0

    sa = np.dot(pred_vec, obs_vec) / (pred_norm * obs_norm)
    sa = np.clip(sa, -1.0, 1.0)
    return float(1.0 - (2.0 / np.pi) * np.arccos(sa))


# ---------------------------------------------------------------------------
# CE calibration grid search
# ---------------------------------------------------------------------------

def calibrate_ce(
    predictor,
    sample_df: pd.DataFrame,
    ce_range: tuple[int, int] = (-10, 10),
    batch_size: int = 512,
    verbose: bool = True,
) -> tuple[float, list[tuple[int, float]], list[float]]:
    """Grid-search CE offset to maximize predicted-vs-observed spectral similarity.

    Args:
        predictor: DeepPeptideIntensityPredictor instance
        sample_df: DataFrame with columns: sequence, charge, collision_energy,
                   ion_type, ion_number, ion_charge, intensity_mono
        ce_range: (lower, upper) CE offset range to search
        batch_size: prediction batch size
        verbose: print progress

    Returns:
        (best_offset, [(offset, mean_similarity), ...], per_psm_sims_at_zero)
    """
    sequences = sample_df["sequence"].tolist()
    charges = sample_df["charge"].tolist()
    base_ces = sample_df["collision_energy"].tolist()

    # Pre-extract observed spectra (list columns)
    observed_data = []
    for _, row in sample_df.iterrows():
        observed_data.append((
            row["ion_type"],
            row["ion_number"],
            row["ion_charge"],
            row["intensity_mono"],
        ))

    results = []
    per_psm_sims_at_zero = None

    for offset in range(ce_range[0], ce_range[1] + 1):
        ces_shifted = [ce + offset for ce in base_ces]

        # Predict at shifted CE (local, fast)
        predicted_spectra = predict_spectra_local(
            predictor, sequences, charges, ces_shifted, batch_size
        )

        # Compute cosine similarity per PSM
        sims = []
        for i, pred in enumerate(predicted_spectra):
            obs_ions, obs_nums, obs_chgs, obs_ints = observed_data[i]
            sim = spectral_angle_fragments(
                pred, obs_ions, obs_nums, obs_chgs, obs_ints
            )
            sims.append(sim if not np.isnan(sim) else 0.0)

        # Save per-PSM similarities at offset=0 for before/after comparison
        if offset == 0:
            per_psm_sims_at_zero = list(sims)

        valid_sims = [s for s in sims if s > 0]
        mean_sim = float(np.mean(valid_sims)) if valid_sims else 0.0
        results.append((offset, mean_sim))

        if verbose:
            print(f"  CE offset {offset:+3d}: mean SA = {mean_sim:.4f} "
                  f"({len(valid_sims)} valid)", flush=True)

    best_idx = int(np.argmax([r[1] for r in results]))
    best_offset = results[best_idx][0]
    best_sim = results[best_idx][1]

    if verbose:
        print(f"\n  Best offset: {best_offset:+d} eV "
              f"(sa = {best_sim:.4f})", flush=True)

    return best_offset, results, per_psm_sims_at_zero


def compute_per_psm_sims(
    predictor,
    sample_df: pd.DataFrame,
    ce_offset: float,
    batch_size: int = 512,
) -> list[float]:
    """Compute per-PSM cosine similarities at a specific CE offset."""
    sequences = sample_df["sequence"].tolist()
    charges = sample_df["charge"].tolist()
    ces = [ce + ce_offset for ce in sample_df["collision_energy"].tolist()]

    predicted_spectra = predict_spectra_local(
        predictor, sequences, charges, ces, batch_size
    )

    sims = []
    for i, pred in enumerate(predicted_spectra):
        row = sample_df.iloc[i]
        sim = spectral_angle_fragments(
            pred,
            row["ion_type"], row["ion_number"],
            row["ion_charge"], row["intensity_mono"],
        )
        sims.append(sim if not np.isnan(sim) else 0.0)

    return sims


# ---------------------------------------------------------------------------
# Report generation
# ---------------------------------------------------------------------------

def generate_report(
    results: dict,
    output_path: Path,
):
    """Generate CE calibration report PDF with before/after cossim distributions."""
    plt.rcParams.update({
        "figure.dpi": 150,
        "font.size": 9,
        "axes.titlesize": 11,
        "axes.labelsize": 10,
    })

    with PdfPages(str(output_path)) as pdf:
        # --- Page 1: Similarity vs CE offset curve ---
        fig, ax = plt.subplots(figsize=(8, 5))

        for dataset_name, data in results.items():
            offsets = [r[0] for r in data["curve"]]
            sims = [r[1] for r in data["curve"]]
            best = data["best_offset"]
            best_sim = sims[offsets.index(best)]

            ax.plot(offsets, sims, label=f'{dataset_name} (best: {best:+d} eV)',
                    linewidth=1.5)
            ax.axvline(best, linestyle="--", alpha=0.4)
            ax.plot(best, best_sim, "o", markersize=6)

        ax.set_xlabel("CE Offset (eV)")
        ax.set_ylabel("Mean Spectral Angle (predicted vs observed)")
        ax.set_title("Collision Energy Calibration: Model vs Observed Spectra")
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)

        # --- Page 2: Before/after cossim distributions ---
        fig, axes = plt.subplots(1, 2, figsize=(12, 5))

        for dataset_name, data in results.items():
            sims_before = np.array(data.get("sims_before", []))
            sims_after = np.array(data.get("sims_after", []))

            if len(sims_before) > 0:
                sims_before_valid = sims_before[sims_before > 0]
                sims_after_valid = sims_after[sims_after > 0]

                # Histogram overlay
                ax = axes[0]
                ax.hist(sims_before_valid, bins=80, alpha=0.5, density=True,
                        label=f"Before (offset=0)", color="#dc2626")
                ax.hist(sims_after_valid, bins=80, alpha=0.5, density=True,
                        label=f"After (offset={data['best_offset']:+d})", color="#2563eb")
                ax.set_xlabel("Spectral Angle")
                ax.set_ylabel("Density")
                ax.set_title(f"Spectral Angle Distribution: {dataset_name}")
                ax.axvline(np.median(sims_before_valid), color="#dc2626",
                           linestyle="--", alpha=0.7,
                           label=f"Median before: {np.median(sims_before_valid):.4f}")
                ax.axvline(np.median(sims_after_valid), color="#2563eb",
                           linestyle="--", alpha=0.7,
                           label=f"Median after: {np.median(sims_after_valid):.4f}")
                ax.legend(fontsize=7)

                # Delta distribution
                ax = axes[1]
                delta = sims_after - sims_before
                delta_nonzero = delta[sims_before > 0]
                ax.hist(delta_nonzero, bins=80, alpha=0.7, color="#6366f1", density=True)
                ax.axvline(0, color="black", linestyle="-", alpha=0.3)
                ax.axvline(np.median(delta_nonzero), color="#6366f1", linestyle="--",
                           label=f"Median delta: {np.median(delta_nonzero):+.4f}")
                improved = (delta_nonzero > 0).sum()
                ax.set_xlabel("Delta Spectral Angle (after - before)")
                ax.set_ylabel("Density")
                ax.set_title(f"Per-PSM Improvement ({improved}/{len(delta_nonzero)} improved)")
                ax.legend(fontsize=8)

        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)

        # --- Page 3: Summary table ---
        fig, ax = plt.subplots(figsize=(8, 4))
        ax.axis("off")

        rows = []
        for dataset_name, data in results.items():
            offsets = [r[0] for r in data["curve"]]
            sims_curve = [r[1] for r in data["curve"]]
            best = data["best_offset"]
            best_sim = sims_curve[offsets.index(best)]
            zero_sim = sims_curve[offsets.index(0)] if 0 in offsets else 0

            sims_before = np.array(data.get("sims_before", []))
            sims_after = np.array(data.get("sims_after", []))
            med_before = np.median(sims_before[sims_before > 0]) if len(sims_before) > 0 else 0
            med_after = np.median(sims_after[sims_after > 0]) if len(sims_after) > 0 else 0

            rows.append((
                dataset_name,
                f"{data['n_psms']:,}",
                f"{best:+d}",
                f"{med_before:.4f}",
                f"{med_after:.4f}",
                f"{med_after - med_before:+.4f}",
            ))

        col_labels = ["Dataset", "N PSMs", "Best Offset",
                       "Median SA (before)", "Median SA (after)", "Improvement"]
        table = ax.table(cellText=rows, colLabels=col_labels, loc="center",
                         cellLoc="center")
        table.auto_set_font_size(False)
        table.set_fontsize(9)
        table.scale(1.0, 1.6)
        for j in range(len(col_labels)):
            table[0, j].set_facecolor("#1e293b")
            table[0, j].set_text_props(color="white", fontweight="bold")

        ax.set_title("CE Calibration Summary", fontsize=13,
                     fontweight="bold", pad=20)
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Calibrate collision energy via local spectral prediction"
    )
    parser.add_argument("--fragments", required=True,
                        help="Path to fragment_peaks.parquet")
    parser.add_argument("--output", default="notebook/analysis/ce_calibration/",
                        help="Output directory")
    parser.add_argument("--n-sample", type=int, default=2000,
                        help="Number of PSMs to sample for calibration (default: 2000)")
    parser.add_argument("--ce-lower", type=int, default=-10,
                        help="Lower CE offset bound (default: -10)")
    parser.add_argument("--ce-upper", type=int, default=10,
                        help="Upper CE offset bound (default: 10)")
    parser.add_argument("--batch-size", type=int, default=512,
                        help="Prediction batch size (default: 512)")
    parser.add_argument("--dataset-name", default=None,
                        help="Dataset name label (default: inferred from path)")
    args = parser.parse_args()

    fragments_path = Path(args.fragments)
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    if not fragments_path.exists():
        print(f"ERROR: Fragment peaks not found: {fragments_path}")
        sys.exit(1)

    # Infer dataset name
    dataset_name = args.dataset_name or fragments_path.parent.name

    # Load fragment peaks
    print(f"Loading fragment peaks from {fragments_path} ...", flush=True)
    needed_cols = [
        "sequence", "charge", "collision_energy", "n_engines",
        "ion_type", "ion_number", "ion_charge", "intensity_mono",
        "n_matched", "coverage_b", "coverage_y",
    ]
    pf = pq.ParquetFile(str(fragments_path))
    df = pf.read(columns=needed_cols).to_pandas()
    print(f"  Loaded {len(df):,} PSMs", flush=True)

    # Filter: need reasonable coverage and valid CE
    mask = (
        (df["n_matched"] >= 5)
        & (df["collision_energy"].notna())
        & (df["collision_energy"] > 0)
        & (df["n_engines"] >= 2)
    )
    df_good = df[mask].copy()
    print(f"  {len(df_good):,} PSMs with >= 5 matched ions and valid CE", flush=True)

    # Clean and validate sequences for the local predictor
    df_good["sequence_clean"] = df_good["sequence"].apply(clean_sequence_for_predictor)
    df_good = df_good[df_good["sequence_clean"].notna()].copy()
    df_good["sequence"] = df_good["sequence_clean"]
    df_good = df_good.drop(columns=["sequence_clean"])
    print(f"  {len(df_good):,} PSMs with predictor-compatible sequences", flush=True)

    # Sample for calibration (stratified by charge, prefer charges 2-4)
    df_good = df_good[df_good["charge"].between(2, 4)]
    n_sample = min(args.n_sample, len(df_good))
    sample_parts = []
    for charge, group in df_good.groupby("charge"):
        n_take = max(1, int(n_sample * len(group) / len(df_good)))
        n_take = min(n_take, len(group))
        sample_parts.append(group.sample(n=n_take, random_state=42))
    sample = pd.concat(sample_parts, ignore_index=False)
    if len(sample) > n_sample:
        sample = sample.sample(n=n_sample, random_state=42)
    sample = sample.reset_index(drop=True)
    print(f"  Sampled {len(sample):,} PSMs for calibration", flush=True)
    print(f"  Charge distribution: "
          f"{sample['charge'].value_counts().sort_index().to_dict()}", flush=True)
    print(f"  CE range: {sample['collision_energy'].min():.1f} - "
          f"{sample['collision_energy'].max():.1f} eV", flush=True)

    # Initialize local predictor
    print(f"\nInitializing local intensity predictor ...", flush=True)
    predictor = get_local_predictor()
    print(f"  Model loaded (local PyTorch, no API calls)", flush=True)

    # Run calibration grid search
    print(f"\nSearching CE offset [{args.ce_lower}, {args.ce_upper}] ...", flush=True)
    best_offset, curve, sims_at_zero = calibrate_ce(
        predictor,
        sample,
        ce_range=(args.ce_lower, args.ce_upper),
        batch_size=args.batch_size,
        verbose=True,
    )

    # Compute per-PSM similarities at best offset (for before/after comparison)
    print(f"\nComputing per-PSM similarities at best offset ({best_offset:+d}) ...",
          flush=True)
    sims_at_best = compute_per_psm_sims(
        predictor, sample, best_offset, args.batch_size
    )

    # Save results
    results = {
        dataset_name: {
            "best_offset": int(best_offset),
            "n_psms": len(sample),
            "curve": [(int(o), float(s)) for o, s in curve],
            "sims_before": sims_at_zero,
            "sims_after": sims_at_best,
        }
    }

    json_path = output_dir / "ce_calibration.json"
    with open(json_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nSaved {json_path}", flush=True)

    # Generate report
    pdf_path = output_dir / "ce_calibration_report.pdf"
    print(f"Generating report: {pdf_path} ...", flush=True)
    generate_report(results, pdf_path)
    print(f"Done. Report saved to {pdf_path}", flush=True)


if __name__ == "__main__":
    main()
