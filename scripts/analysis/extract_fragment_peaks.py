#!/usr/bin/env python3
"""
Search-Engine-Independent Fragment Peak Extractor

For high-confidence PSMs: generates theoretical b/y ions via imspy, visits the
raw PASEF spectrum (from blob), and extracts monoisotopic + M+1/M+2 intensities.
Annotates with collision energy for downstream CE-dependent analysis.

Output:
    - fragment_peaks.parquet: one row per PSM with list columns for fragment arrays
    - fragment_peaks_report.pdf: 6-page diagnostic report
"""

import argparse
import gzip
import io
import json
import re
import sys
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

# Add project paths
PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))

from fragment_matching import FragmentMatcher, MatchConfig, sage_to_imspy_sequence

try:
    import zstandard as zstd
    HAS_ZSTD = True
except ImportError:
    HAS_ZSTD = False


# Neutron mass for isotope spacing (C13-C12)
NEUTRON_MASS = 1.003355

# Valid amino acid characters for imspy
_VALID_AA = set("ARNDCQEGHILKMFPSTWYVU")


def clean_sequence_for_imspy(seq: str) -> Optional[str]:
    """Clean a modified peptide sequence so imspy can parse it.

    Handles:
    - Removes hyphens from N/C-term mods (e.g. [UNIMOD:1]-M -> [UNIMOD:1]M)
    - Validates only valid amino acids + [UNIMOD:X] remain

    Returns None if the sequence cannot be cleaned.
    """
    if not seq or not isinstance(seq, str):
        return None

    # Remove hyphens (from N-term/C-term separator notation)
    cleaned = seq.replace("-", "")

    # Validate: after stripping [UNIMOD:X] tags, only valid AAs should remain
    stripped = re.sub(r"\[UNIMOD:\d+\]", "", cleaned)
    if not stripped:  # empty after stripping mods
        return None
    if not all(c in _VALID_AA for c in stripped):
        return None

    return cleaned


# ---------------------------------------------------------------------------
# Blob reading (adapted from dashboard/backend/main.py)
# ---------------------------------------------------------------------------

def read_blob(blob_path: Path, offset: int, size: int) -> Optional[dict]:
    """Read and decompress a precursor blob from blobs.bin.

    Returns dict with arrays: frag_mz, frag_intensity (numpy arrays).
    """
    if not blob_path.exists():
        return None

    try:
        with open(blob_path, "rb") as f:
            file_size = f.seek(0, 2)
            if offset + size > file_size:
                return None
            f.seek(offset)
            compressed = f.read(size)

        if len(compressed) == 0:
            return None

        # Detect compression
        if compressed[:4] == b"\x28\xb5\x2f\xfd":
            if not HAS_ZSTD:
                return None
            dctx = zstd.ZstdDecompressor()
            combined = dctx.decompress(compressed)
        elif compressed[:2] == b"\x1f\x8b":
            combined = gzip.decompress(compressed)
        else:
            return None

        # Parse: [4 bytes metadata_len][metadata JSON][npz arrays]
        metadata_len = int.from_bytes(combined[:4], "little")
        npz_bytes = combined[4 + metadata_len:]
        arrays = np.load(io.BytesIO(npz_bytes))

        return {
            "frag_mz": arrays["frag_mz"] if "frag_mz" in arrays else np.array([], dtype=np.float32),
            "frag_intensity": arrays["frag_intensity"] if "frag_intensity" in arrays else np.array([], dtype=np.float32),
        }
    except Exception:
        return None


def resolve_blob_path(blob_dir: Path, raw_file: str) -> Path:
    """Resolve blobs.bin path for a raw file."""
    raw_clean = raw_file.replace(".d", "")
    return blob_dir / f"{raw_clean}.d" / "blobs.bin"


# ---------------------------------------------------------------------------
# Fragment peak extraction
# ---------------------------------------------------------------------------

def extract_peaks_for_psm(
    frag_mz: np.ndarray,
    frag_intensity: np.ndarray,
    theoretical_fragments: dict,
    tolerance_ppm: float = 20.0,
) -> dict:
    """Extract monoisotopic + M+1/M+2 intensities for each theoretical ion.

    Args:
        frag_mz: Raw observed m/z array (sorted)
        frag_intensity: Raw observed intensity array
        theoretical_fragments: Dict from FragmentMatcher.generate_theoretical_fragments()
        tolerance_ppm: Mass tolerance in ppm

    Returns:
        Dict with per-ion arrays and summary metrics.
    """
    ion_types = []
    ion_numbers = []
    ion_charges = []
    mz_theoretical = []
    mz_matched = []
    intensity_mono = []
    intensity_m1 = []
    intensity_m2 = []
    error_ppm = []

    n_matched = 0
    n_theoretical = 0
    n_b_theo = 0
    n_y_theo = 0
    n_b_matched = 0
    n_y_matched = 0
    matched_intensity_sum = 0.0
    total_intensity = float(np.sum(frag_intensity)) if len(frag_intensity) > 0 else 0.0

    # Ensure mz is sorted for binary search
    if len(frag_mz) > 0:
        sort_idx = np.argsort(frag_mz)
        frag_mz = frag_mz[sort_idx]
        frag_intensity = frag_intensity[sort_idx]

    for ion_type, ions in theoretical_fragments.items():
        for ion_num, ion_charge, theo_mz, _seq in ions:
            n_theoretical += 1
            if ion_type == "b":
                n_b_theo += 1
            else:
                n_y_theo += 1

            ion_types.append(ion_type)
            ion_numbers.append(ion_num)
            ion_charges.append(ion_charge)
            mz_theoretical.append(theo_mz)

            # Search for monoisotopic peak
            mono_int, mono_mz, mono_err = _find_peak(
                frag_mz, frag_intensity, theo_mz, tolerance_ppm
            )

            # M+1 isotope
            m1_target = theo_mz + NEUTRON_MASS / ion_charge
            m1_int, _, _ = _find_peak(frag_mz, frag_intensity, m1_target, tolerance_ppm)

            # M+2 isotope
            m2_target = theo_mz + 2 * NEUTRON_MASS / ion_charge
            m2_int, _, _ = _find_peak(frag_mz, frag_intensity, m2_target, tolerance_ppm)

            mz_matched.append(mono_mz)
            intensity_mono.append(mono_int)
            intensity_m1.append(m1_int)
            intensity_m2.append(m2_int)
            error_ppm.append(mono_err)

            if mono_int > 0:
                n_matched += 1
                matched_intensity_sum += mono_int
                if ion_type == "b":
                    n_b_matched += 1
                else:
                    n_y_matched += 1

    coverage_b = n_b_matched / n_b_theo if n_b_theo > 0 else 0.0
    coverage_y = n_y_matched / n_y_theo if n_y_theo > 0 else 0.0
    intensity_explained = matched_intensity_sum / total_intensity if total_intensity > 0 else 0.0

    # Median absolute ppm error for matched ions
    matched_errors = [abs(e) for e, m in zip(error_ppm, intensity_mono) if m > 0 and not np.isnan(e)]
    median_ppm = float(np.median(matched_errors)) if matched_errors else np.nan

    return {
        "ion_type": ion_types,
        "ion_number": ion_numbers,
        "ion_charge": ion_charges,
        "mz_theoretical": mz_theoretical,
        "mz_matched": mz_matched,
        "intensity_mono": intensity_mono,
        "intensity_m1": intensity_m1,
        "intensity_m2": intensity_m2,
        "error_ppm": error_ppm,
        "n_matched": n_matched,
        "n_theoretical": n_theoretical,
        "coverage_b": coverage_b,
        "coverage_y": coverage_y,
        "intensity_explained": intensity_explained,
        "median_ppm_error": median_ppm,
    }


def _find_peak(
    mz_array: np.ndarray,
    int_array: np.ndarray,
    target_mz: float,
    tolerance_ppm: float,
) -> tuple[float, float, float]:
    """Find closest peak within ppm tolerance.

    Returns (intensity, observed_mz, error_ppm). All 0.0 if not found.
    """
    if len(mz_array) == 0:
        return 0.0, 0.0, np.nan

    tol_da = target_mz * tolerance_ppm * 1e-6

    # Binary search for window
    lo = np.searchsorted(mz_array, target_mz - tol_da)
    hi = np.searchsorted(mz_array, target_mz + tol_da, side="right")

    if lo >= hi:
        return 0.0, 0.0, np.nan

    # Find closest peak in window
    window_mz = mz_array[lo:hi]
    window_int = int_array[lo:hi]
    errors = np.abs(window_mz - target_mz)
    best_idx = np.argmin(errors)

    obs_mz = float(window_mz[best_idx])
    obs_int = float(window_int[best_idx])
    err_ppm = (obs_mz - target_mz) / target_mz * 1e6

    return obs_int, obs_mz, err_ppm


# ---------------------------------------------------------------------------
# Batch processing
# ---------------------------------------------------------------------------

def get_sequence_for_psm(row: pd.Series) -> Optional[str]:
    """Pick the best modified sequence for theoretical fragment generation.

    Prefers sage_modified > consensus_peptide.
    Cleans sequences to be imspy-compatible (removes hyphens, validates AAs).
    Returns None if no usable sequence.
    """
    # sage_modified already has [UNIMOD:X] format in our store
    sage_mod = row.get("sage_modified")
    if pd.notna(sage_mod) and isinstance(sage_mod, str) and len(sage_mod) > 0:
        cleaned = clean_sequence_for_imspy(sage_mod)
        if cleaned is not None:
            return cleaned

    # consensus_peptide may also have [UNIMOD:X] mods
    consensus = row.get("consensus_peptide")
    if pd.notna(consensus) and isinstance(consensus, str) and len(consensus) > 0:
        cleaned = clean_sequence_for_imspy(consensus)
        if cleaned is not None:
            return cleaned

    return None


def process_batch(
    batch_df: pd.DataFrame,
    blob_dir: Path,
    matcher: FragmentMatcher,
    tolerance_ppm: float = 20.0,
) -> list[dict]:
    """Process a batch of PSMs, extracting fragment peaks from blobs.

    Returns list of result dicts (one per successfully processed PSM).
    """
    results = []

    # Group by raw_file to minimize blob file opens
    for raw_file, group in batch_df.groupby("raw_file"):
        blob_path = resolve_blob_path(blob_dir, raw_file)
        if not blob_path.exists():
            continue

        for _, row in group.iterrows():
            offset = int(row["blob_offset"])
            size = int(row["blob_size"])

            if size <= 0:
                continue

            # Read blob
            blob = read_blob(blob_path, offset, size)
            if blob is None:
                continue

            frag_mz = blob["frag_mz"]
            frag_intensity = blob["frag_intensity"]
            if len(frag_mz) == 0:
                continue

            # Get modified sequence
            sequence = get_sequence_for_psm(row)
            if sequence is None:
                continue

            charge = int(row["charge"])

            # Generate theoretical fragments
            try:
                theoretical = matcher.generate_theoretical_fragments(
                    sequence, charge
                )
            except Exception:
                continue

            # Extract peaks
            peaks = extract_peaks_for_psm(
                frag_mz, frag_intensity, theoretical, tolerance_ppm
            )

            result = {
                "precursor_id": int(row["precursor_id"]),
                "raw_file": raw_file,
                "sequence": sequence,
                "charge": charge,
                "collision_energy": float(row["collision_energy"]) if pd.notna(row.get("collision_energy")) else np.nan,
                "mobility": float(row["mobility"]) if pd.notna(row.get("mobility")) else np.nan,
                "n_engines": int(row["n_engines"]),
                **peaks,
            }
            results.append(result)

    return results


def build_pyarrow_table(results: list[dict]) -> pa.Table:
    """Convert results list to a PyArrow table with list columns."""
    # Separate scalar and list columns
    scalar_cols = [
        "precursor_id", "raw_file", "sequence", "charge",
        "collision_energy", "mobility", "n_engines",
        "n_matched", "n_theoretical", "coverage_b", "coverage_y",
        "intensity_explained", "median_ppm_error",
    ]
    list_cols = [
        "ion_type", "ion_number", "ion_charge",
        "mz_theoretical", "mz_matched",
        "intensity_mono", "intensity_m1", "intensity_m2",
        "error_ppm",
    ]

    arrays = {}
    for col in scalar_cols:
        arrays[col] = [r[col] for r in results]
    for col in list_cols:
        arrays[col] = [r[col] for r in results]

    # Define schema
    fields = [
        pa.field("precursor_id", pa.int64()),
        pa.field("raw_file", pa.string()),
        pa.field("sequence", pa.string()),
        pa.field("charge", pa.int32()),
        pa.field("collision_energy", pa.float32()),
        pa.field("mobility", pa.float32()),
        pa.field("n_engines", pa.int32()),
        pa.field("ion_type", pa.list_(pa.string())),
        pa.field("ion_number", pa.list_(pa.int32())),
        pa.field("ion_charge", pa.list_(pa.int32())),
        pa.field("mz_theoretical", pa.list_(pa.float64())),
        pa.field("mz_matched", pa.list_(pa.float64())),
        pa.field("intensity_mono", pa.list_(pa.float32())),
        pa.field("intensity_m1", pa.list_(pa.float32())),
        pa.field("intensity_m2", pa.list_(pa.float32())),
        pa.field("error_ppm", pa.list_(pa.float32())),
        pa.field("n_matched", pa.int32()),
        pa.field("n_theoretical", pa.int32()),
        pa.field("coverage_b", pa.float32()),
        pa.field("coverage_y", pa.float32()),
        pa.field("intensity_explained", pa.float32()),
        pa.field("median_ppm_error", pa.float32()),
    ]
    schema = pa.schema(fields)

    pa_arrays = []
    for f in fields:
        col = f.name
        if pa.types.is_list(f.type):
            inner_type = f.type.value_type
            pa_arrays.append(pa.array(arrays[col], type=f.type))
        else:
            pa_arrays.append(pa.array(arrays[col], type=f.type))

    return pa.table(pa_arrays, schema=schema)


# ---------------------------------------------------------------------------
# Report generation (6-page PDF)
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
    """Generate 6-page fragment peak extraction report."""
    _set_style()

    with PdfPages(str(output_path)) as pdf:
        # --- Page 1: Collection summary ---
        fig, axes = plt.subplots(2, 2, figsize=(10, 7))

        # Top-left: summary text
        ax = axes[0, 0]
        ax.axis("off")
        n_spectra = len(df)
        n_peptides = df["sequence"].nunique()
        ce_min = df["collision_energy"].min()
        ce_max = df["collision_energy"].max()
        med_coverage_b = df["coverage_b"].median()
        med_coverage_y = df["coverage_y"].median()
        med_matched = df["n_matched"].median()
        med_ppm = df["median_ppm_error"].median()
        text = (
            f"Fragment Peak Collection Summary\n"
            f"{'='*40}\n\n"
            f"Total spectra:        {n_spectra:,}\n"
            f"Unique peptides:      {n_peptides:,}\n"
            f"CE range:             {ce_min:.1f} - {ce_max:.1f} eV\n"
            f"Raw files:            {df['raw_file'].nunique()}\n\n"
            f"Median b coverage:    {med_coverage_b:.1%}\n"
            f"Median y coverage:    {med_coverage_y:.1%}\n"
            f"Median ions matched:  {med_matched:.0f}\n"
            f"Median ppm error:     {med_ppm:.1f}\n"
        )
        ax.text(0.05, 0.95, text, transform=ax.transAxes, fontsize=9,
                verticalalignment="top", fontfamily="monospace")

        # Top-right: charge distribution
        ax = axes[0, 1]
        charge_counts = df["charge"].value_counts().sort_index()
        ax.bar(charge_counts.index.astype(str) + "+", charge_counts.values,
               color="#6366f1")
        ax.set_xlabel("Charge")
        ax.set_ylabel("Count")
        ax.set_title("Charge Distribution")

        # Bottom-left: n_engines distribution
        ax = axes[1, 0]
        eng_counts = df["n_engines"].value_counts().sort_index()
        ax.bar(eng_counts.index.astype(str), eng_counts.values, color="#10b981")
        ax.set_xlabel("Engine Consensus")
        ax.set_ylabel("Count")
        ax.set_title("Engine Agreement")

        # Bottom-right: matched vs theoretical
        ax = axes[1, 1]
        ax.hist(df["n_matched"], bins=50, alpha=0.7, label="Matched", color="#2563eb")
        ax.hist(df["n_theoretical"], bins=50, alpha=0.4, label="Theoretical", color="#dc2626")
        ax.set_xlabel("Number of Ions")
        ax.set_ylabel("Count")
        ax.set_title("Matched vs Theoretical Ions")
        ax.legend()

        fig.suptitle("Fragment Peak Extraction Report", fontsize=14, fontweight="bold")
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)

        # --- Page 2: CE distribution ---
        fig, ax = plt.subplots()
        ce_valid = df["collision_energy"].dropna()
        ax.hist(ce_valid, bins=100, color="#f59e0b", edgecolor="none")
        ax.set_xlabel("Collision Energy (eV)")
        ax.set_ylabel("Count")
        ax.set_title(f"Collision Energy Distribution (n={len(ce_valid):,})")
        ax.axvline(ce_valid.median(), color="red", linestyle="--",
                   label=f"Median = {ce_valid.median():.1f} eV")
        ax.legend()
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)

        # --- Page 3: Fragment coverage vs CE ---
        fig, axes = plt.subplots(1, 2, figsize=(10, 5))

        for i, (col, label) in enumerate([("coverage_b", "b-ion Coverage"),
                                           ("coverage_y", "y-ion Coverage")]):
            ax = axes[i]
            valid = df[["collision_energy", col]].dropna()
            if len(valid) > 0:
                hb = ax.hexbin(valid["collision_energy"], valid[col],
                               gridsize=40, cmap="YlOrRd", mincnt=1)
                fig.colorbar(hb, ax=ax, label="Count")
            ax.set_xlabel("Collision Energy (eV)")
            ax.set_ylabel(label)
            ax.set_title(f"{label} vs CE")

        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)

        # --- Page 4: Isotope ratio analysis ---
        fig, axes = plt.subplots(1, 2, figsize=(10, 5))

        # Explode list columns for isotope analysis
        # Sample up to 50k rows for plotting performance
        sample_df = df.sample(n=min(50000, len(df)), random_state=42) if len(df) > 50000 else df

        # Compute per-ion M+1/mono ratios
        ratios_mz = []
        ratios_val = []
        for _, row in sample_df.iterrows():
            monos = row["intensity_mono"]
            m1s = row["intensity_m1"]
            mzs = row["mz_theoretical"]
            if not isinstance(monos, (list, np.ndarray)):
                continue
            for mono, m1, mz in zip(monos, m1s, mzs):
                if mono > 0 and m1 > 0:
                    ratios_mz.append(mz)
                    ratios_val.append(m1 / mono)

        ax = axes[0]
        if ratios_mz:
            ax.hexbin(ratios_mz, ratios_val, gridsize=50, cmap="viridis",
                      mincnt=1, extent=[0, 1700, 0, 2.0])
            ax.set_xlabel("Fragment m/z")
            ax.set_ylabel("M+1 / Monoisotopic Ratio")
            ax.set_title("Isotope Ratio vs Fragment m/z")
        else:
            ax.text(0.5, 0.5, "No isotope data", transform=ax.transAxes, ha="center")

        # M+2/mono ratio
        ratios_mz2 = []
        ratios_val2 = []
        for _, row in sample_df.iterrows():
            monos = row["intensity_mono"]
            m2s = row["intensity_m2"]
            mzs = row["mz_theoretical"]
            if not isinstance(monos, (list, np.ndarray)):
                continue
            for mono, m2, mz in zip(monos, m2s, mzs):
                if mono > 0 and m2 > 0:
                    ratios_mz2.append(mz)
                    ratios_val2.append(m2 / mono)

        ax = axes[1]
        if ratios_mz2:
            ax.hexbin(ratios_mz2, ratios_val2, gridsize=50, cmap="viridis",
                      mincnt=1, extent=[0, 1700, 0, 1.0])
            ax.set_xlabel("Fragment m/z")
            ax.set_ylabel("M+2 / Monoisotopic Ratio")
            ax.set_title("M+2 Isotope Ratio vs Fragment m/z")
        else:
            ax.text(0.5, 0.5, "No M+2 data", transform=ax.transAxes, ha="center")

        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)

        # --- Page 5: PPM error distribution ---
        fig, axes = plt.subplots(1, 2, figsize=(10, 5))

        ax = axes[0]
        ppm_valid = df["median_ppm_error"].dropna()
        if len(ppm_valid) > 0:
            ax.hist(ppm_valid.clip(-20, 20), bins=100, color="#8b5cf6", edgecolor="none")
            ax.set_xlabel("Median Absolute PPM Error")
            ax.set_ylabel("Count")
            ax.set_title(f"PPM Error Distribution (median={ppm_valid.median():.2f})")
            ax.axvline(ppm_valid.median(), color="red", linestyle="--")

        # Per-ion ppm error (sampled)
        ax = axes[1]
        all_errors = []
        for _, row in sample_df.iterrows():
            errs = row["error_ppm"]
            monos = row["intensity_mono"]
            if not isinstance(errs, (list, np.ndarray)):
                continue
            for e, m in zip(errs, monos):
                if m > 0 and not np.isnan(e):
                    all_errors.append(e)

        if all_errors:
            all_errors = np.array(all_errors)
            ax.hist(all_errors.clip(-20, 20), bins=200, color="#06b6d4", edgecolor="none")
            ax.set_xlabel("PPM Error (signed)")
            ax.set_ylabel("Count")
            ax.set_title(f"Per-Ion PPM Error (n={len(all_errors):,})")
            ax.axvline(0, color="black", linestyle="-", alpha=0.3)
            ax.axvline(np.median(all_errors), color="red", linestyle="--",
                       label=f"Median = {np.median(all_errors):.2f}")
            ax.legend()

        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)

        # --- Page 6: Intensity explained vs quality ---
        fig, axes = plt.subplots(1, 2, figsize=(10, 5))

        ax = axes[0]
        for n_eng in sorted(df["n_engines"].unique()):
            subset = df[df["n_engines"] == n_eng]["intensity_explained"]
            if len(subset) > 0:
                ax.hist(subset, bins=50, alpha=0.5,
                        label=f"{n_eng} engines (n={len(subset):,})", density=True)
        ax.set_xlabel("Intensity Explained")
        ax.set_ylabel("Density")
        ax.set_title("Intensity Explained by Engine Consensus")
        ax.legend(fontsize=7)

        # Coverage vs n_engines
        ax = axes[1]
        eng_vals = sorted(df["n_engines"].unique())
        coverage_data = [df[df["n_engines"] == n]["coverage_b"].values for n in eng_vals]
        if coverage_data:
            parts = ax.violinplot(coverage_data, positions=range(len(eng_vals)),
                                  showmedians=True, showextrema=False)
            for pc in parts["bodies"]:
                pc.set_facecolor("#10b981")
                pc.set_alpha(0.6)
            ax.set_xticks(range(len(eng_vals)))
            ax.set_xticklabels([str(n) for n in eng_vals])
        ax.set_xlabel("Engine Consensus (n_engines)")
        ax.set_ylabel("b-ion Coverage")
        ax.set_title("Fragment Coverage by Engine Agreement")

        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Search-engine-independent fragment peak extractor"
    )
    parser.add_argument("--store", required=True,
                        help="Path to precursor_store.parquet")
    parser.add_argument("--blob-dir", required=True,
                        help="Directory containing extracted/{raw_file}.d/blobs.bin")
    parser.add_argument("--output", default="notebook/analysis/fragment_collection/",
                        help="Output directory")
    parser.add_argument("--min-engines", type=int, default=2,
                        help="Minimum engine consensus for PSM selection (default: 2)")
    parser.add_argument("--batch-size", type=int, default=5000,
                        help="PSMs per processing batch (default: 5000)")
    parser.add_argument("--tolerance-ppm", type=float, default=20.0,
                        help="m/z tolerance in ppm (default: 20.0)")
    parser.add_argument("--max-psms", type=int, default=None,
                        help="Maximum PSMs to process (for testing)")
    args = parser.parse_args()

    store_path = Path(args.store)
    blob_dir = Path(args.blob_dir)
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    if not store_path.exists():
        print(f"ERROR: Store not found: {store_path}")
        sys.exit(1)
    if not blob_dir.exists():
        print(f"ERROR: Blob directory not found: {blob_dir}")
        sys.exit(1)

    # Load PSM selection columns — row-group-at-a-time to limit memory
    print(f"Loading precursor data from {store_path} ...", flush=True)
    select_cols = [
        "precursor_id", "raw_file", "charge", "mobility",
        "collision_energy", "n_engines", "is_high_quality",
        "blob_offset", "blob_size",
        "consensus_peptide", "sage_modified",
    ]
    pf = pq.ParquetFile(str(store_path))
    chunks = []
    for rg in range(pf.metadata.num_row_groups):
        rg_df = pf.read_row_group(rg, columns=select_cols).to_pandas()
        mask = (
            (rg_df["n_engines"] >= args.min_engines)
            & (rg_df["is_high_quality"] == True)
            & (rg_df["blob_size"] > 0)
        )
        chunks.append(rg_df[mask])
        print(f"  Row group {rg + 1}/{pf.metadata.num_row_groups}: "
              f"{mask.sum():,} / {len(rg_df):,} pass filter", flush=True)
    df_sel = pd.concat(chunks, ignore_index=True)
    del chunks
    print(f"  Selected {len(df_sel):,} high-confidence PSMs "
          f"(n_engines >= {args.min_engines}, is_high_quality, blob_size > 0)",
          flush=True)

    if args.max_psms and len(df_sel) > args.max_psms:
        df_sel = df_sel.sample(n=args.max_psms, random_state=42)
        print(f"  Subsampled to {len(df_sel):,} PSMs (--max-psms)", flush=True)

    if len(df_sel) == 0:
        print("ERROR: No PSMs pass selection criteria.")
        sys.exit(1)

    # Initialize fragment matcher
    config = MatchConfig(
        mz_tolerance_ppm=args.tolerance_ppm,
        ion_types=["b", "y"],
        max_fragment_charge=2,
    )
    matcher = FragmentMatcher(config)
    print(f"  Fragment matcher: {args.tolerance_ppm} ppm, b/y ions, max charge 2",
          flush=True)

    # Process in batches
    n_total = len(df_sel)
    n_batches = (n_total + args.batch_size - 1) // args.batch_size
    all_results = []
    n_processed = 0
    n_failed = 0

    print(f"\nProcessing {n_total:,} PSMs in {n_batches} batches ...", flush=True)

    for batch_idx in range(n_batches):
        start = batch_idx * args.batch_size
        end = min(start + args.batch_size, n_total)
        batch = df_sel.iloc[start:end]

        batch_results = process_batch(batch, blob_dir, matcher, args.tolerance_ppm)
        all_results.extend(batch_results)

        n_processed += len(batch)
        batch_ok = len(batch_results)
        n_failed += len(batch) - batch_ok

        if (batch_idx + 1) % 10 == 0 or batch_idx == n_batches - 1:
            print(f"  Batch {batch_idx + 1}/{n_batches}: "
                  f"{n_processed:,}/{n_total:,} processed, "
                  f"{len(all_results):,} extracted, "
                  f"{n_failed:,} failed", flush=True)

    print(f"\nExtraction complete: {len(all_results):,} PSMs with fragment peaks",
          flush=True)

    if len(all_results) == 0:
        print("ERROR: No results extracted.")
        sys.exit(1)

    # Save parquet
    print("Building PyArrow table ...", flush=True)
    table = build_pyarrow_table(all_results)

    parquet_path = output_dir / "fragment_peaks.parquet"
    pq.write_table(table, str(parquet_path), compression="zstd")
    print(f"Saved {parquet_path} ({len(all_results):,} rows)", flush=True)

    # Generate report
    # Convert to pandas for plotting (list columns preserved)
    report_df = table.to_pandas()

    pdf_path = output_dir / "fragment_peaks_report.pdf"
    print(f"Generating report: {pdf_path} ...", flush=True)
    generate_report(report_df, pdf_path)
    print(f"Done. Report saved to {pdf_path}", flush=True)


if __name__ == "__main__":
    main()
