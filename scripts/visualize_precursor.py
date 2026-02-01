#!/usr/bin/env python3
"""
Precursor Visualization Tool

Visualizes extracted precursors with their search engine identifications.
Shows: XIC (chromatogram), mobilogram, isotope envelope, and fragment spectrum.

Usage:
    # Using unified index (recommended)
    python visualize_precursor.py --index data/processed/PXD019086/precursor_index.parquet \\
        --raw /path/to/data.d --precursor-id 68381

    # Interactive mode - browse precursors
    python visualize_precursor.py --index data/processed/PXD019086/precursor_index.parquet \\
        --raw /path/to/data.d --interactive

    # Filter to show disagreements between engines
    python visualize_precursor.py --index data/processed/PXD019086/precursor_index.parquet \\
        --raw /path/to/data.d --show-disagreements
"""

import sys
import argparse
from pathlib import Path
from typing import Optional, Dict, Any, Tuple
import json
import io

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import matplotlib.patches as mpatches

# imspy imports
from imspy_core.timstof import TimsDatasetDDA
from imspy_core.timstof.frame import TimsFrame
from imspy_connector import py_dda

# Optional zstd compression
try:
    import zstandard as zstd
    HAS_ZSTD = True
except ImportError:
    import gzip
    HAS_ZSTD = False


def deserialize_blob(data: bytes) -> dict:
    """Deserialize blob to dict with arrays."""
    if HAS_ZSTD:
        dctx = zstd.ZstdDecompressor()
        combined = dctx.decompress(data)
    else:
        combined = gzip.decompress(data)

    metadata_len = int.from_bytes(combined[:4], "little")
    metadata_bytes = combined[4:4+metadata_len]
    metadata = json.loads(metadata_bytes.decode("utf-8"))

    npz_bytes = combined[4+metadata_len:]
    npz_buffer = io.BytesIO(npz_bytes)
    arrays = dict(np.load(npz_buffer))

    return {"metadata": metadata, "arrays": arrays}


def load_blob_by_precursor_id(
    blobs_path: Path,
    index_df: pd.DataFrame,
    precursor_id: int
) -> Optional[dict]:
    """Load a single precursor blob by ID."""
    row = index_df[index_df["precursor_id"] == precursor_id]
    if len(row) == 0:
        return None

    row = row.iloc[0]
    offset = int(row["blob_offset"])
    size = int(row["blob_size"])

    with open(blobs_path, "rb") as f:
        f.seek(offset)
        data = f.read(size)

    return deserialize_blob(data)


def load_psm_data(psm_path: Path) -> pd.DataFrame:
    """Load FragPipe PSM data."""
    df = pd.read_csv(psm_path, sep="\t")

    # Parse precursor ID from Spectrum column
    # Format: rawfile.scannum.scannum.charge
    def parse_precursor_id(spectrum):
        parts = spectrum.split(".")
        if len(parts) >= 3:
            return int(parts[-2])
        return None

    df["precursor_id"] = df["Spectrum"].apply(parse_precursor_id)
    return df


def load_diann_data(diann_path: Path) -> pd.DataFrame:
    """Load DIA-NN report data."""
    df = pd.read_parquet(diann_path)
    # DIA-NN uses Precursor.Id format: SEQUENCE_CHARGE
    return df


def load_sage_data(sage_path: Path) -> pd.DataFrame:
    """Load Sage results data."""
    df = pd.read_parquet(sage_path)
    return df


class PrecursorVisualizer:
    """Visualize precursor data with search engine identifications."""

    def __init__(
        self,
        raw_data_path: Path,
        psm_path: Optional[Path] = None,
        diann_path: Optional[Path] = None,
        sage_path: Optional[Path] = None,
        extracted_path: Optional[Path] = None,
        use_bruker_sdk: bool = False,
    ):
        """
        Initialize visualizer.

        Args:
            raw_data_path: Path to .d folder
            psm_path: Path to FragPipe psm.tsv
            diann_path: Path to DIA-NN report.parquet
            sage_path: Path to Sage results.sage.parquet
            extracted_path: Path to extracted data directory (index.parquet + blobs.bin)
            use_bruker_sdk: Whether to use Bruker SDK for raw data access
        """
        self.raw_data_path = Path(raw_data_path)

        # Load dataset
        print(f"Loading dataset: {self.raw_data_path.name}")
        self.dataset = TimsDatasetDDA(
            str(self.raw_data_path),
            in_memory=False,
            use_bruker_sdk=use_bruker_sdk,
        )

        # Get precursor metadata
        self.precursors_df = self.dataset.get_pasef_fragments(num_threads=1)
        self.precursor_meta = self.dataset.fragmented_precursors

        # Build precursor ID -> metadata lookup
        self.precursor_lookup = {}
        for _, row in self.precursor_meta.iterrows():
            self.precursor_lookup[row["precursor_id"]] = row.to_dict()

        # Load search engine results
        self.psm_df = None
        self.diann_df = None
        self.sage_df = None

        if psm_path and Path(psm_path).exists():
            print(f"Loading FragPipe PSMs: {psm_path}")
            self.psm_df = load_psm_data(Path(psm_path))
            print(f"  Loaded {len(self.psm_df)} PSMs")

        if diann_path and Path(diann_path).exists():
            print(f"Loading DIA-NN results: {diann_path}")
            self.diann_df = load_diann_data(Path(diann_path))
            print(f"  Loaded {len(self.diann_df)} precursors")

        if sage_path and Path(sage_path).exists():
            print(f"Loading Sage results: {sage_path}")
            self.sage_df = load_sage_data(Path(sage_path))
            print(f"  Loaded {len(self.sage_df)} PSMs")

        # Load extracted data if available
        self.extracted_index = None
        self.extracted_blobs_path = None

        if extracted_path:
            extracted_path = Path(extracted_path)
            index_path = extracted_path / "index.parquet"
            blobs_path = extracted_path / "blobs.bin"

            if index_path.exists():
                print(f"Loading extracted index: {index_path}")
                self.extracted_index = pd.read_parquet(index_path)
                self.extracted_blobs_path = blobs_path
                print(f"  Loaded {len(self.extracted_index)} precursors")

    def get_precursor_identification(self, precursor_id: int) -> Dict[str, Any]:
        """Get search engine identification for a precursor."""
        result = {
            "fragpipe": None,
            "diann": None,
            "sage": None,
        }

        if self.psm_df is not None:
            psm = self.psm_df[self.psm_df["precursor_id"] == precursor_id]
            if len(psm) > 0:
                # Take best PSM (highest probability)
                best = psm.loc[psm["Probability"].idxmax()]
                result["fragpipe"] = {
                    "peptide": best.get("Modified Peptide") or best.get("Peptide", ""),
                    "protein": best.get("Protein", ""),
                    "probability": best.get("Probability", 0),
                    "qvalue": best.get("Qvalue", 1),
                    "hyperscore": best.get("Hyperscore", 0),
                }

        # For DIA-NN, we need to match by m/z and RT window
        if self.diann_df is not None and precursor_id in self.precursor_lookup:
            meta = self.precursor_lookup[precursor_id]
            mz = meta.get("mono_mz") or meta.get("largest_peak_mz", 0)
            charge = meta.get("charge", 0)

            if mz > 0 and charge > 0:
                # Match by m/z (within 20 ppm) and charge
                ppm_tol = 20
                mz_tol = mz * ppm_tol / 1e6

                matches = self.diann_df[
                    (self.diann_df["Precursor.Charge"] == charge) &
                    (abs(self.diann_df["Precursor.Mz"] - mz) < mz_tol)
                ]

                if len(matches) > 0:
                    best = matches.iloc[0]
                    result["diann"] = {
                        "peptide": best.get("Modified.Sequence", ""),
                        "protein": best.get("Protein.Ids", ""),
                        "qvalue": best.get("Q.Value", 1) if "Q.Value" in best else None,
                    }

        # Sage matching (by scan number in scannr column)
        if self.sage_df is not None:
            # Sage uses scannr which should match precursor_id
            if "scannr" in self.sage_df.columns:
                sage_match = self.sage_df[self.sage_df["scannr"] == precursor_id]
                if len(sage_match) > 0:
                    best = sage_match.iloc[0]
                    result["sage"] = {
                        "peptide": best.get("peptide", ""),
                        "protein": best.get("proteins", ""),
                        "hyperscore": best.get("hyperscore", 0),
                        "qvalue": best.get("spectrum_q", 1),
                    }

        return result

    def extract_ms1_signal(
        self,
        precursor_id: int,
        rt_window_sec: float = 30.0,
        mz_tol_ppm: float = 20.0,
        im_window: float = 0.1,
        n_isotopes: int = 5,
    ) -> Optional[Dict[str, np.ndarray]]:
        """Extract MS1 signal for a precursor using Rust backend."""
        if precursor_id not in self.precursor_lookup:
            return None

        meta = self.precursor_lookup[precursor_id]

        # Get precursor properties
        mz = meta.get("mono_mz") or meta.get("largest_peak_mz", 0)
        charge = meta.get("charge") or 2

        # Get RT from PASEF frame
        frame_id = meta.get("frame_id", 0)
        frame_meta = self.dataset.get_frame(int(frame_id))
        rt_sec = frame_meta.retention_time * 60.0  # Convert to seconds

        # Get mobility from precursor
        mobility = meta.get("inverse_ion_mobility", 0.9)

        # Build coordinate
        coord = py_dda.PyPrecursorCoord(
            precursor_id=int(precursor_id),
            mz=float(mz),
            rt_seconds=float(rt_sec),
            mobility=float(mobility),
            charge=int(charge),
        )

        # Extract using Rust
        rust_dataset = self.dataset.get_py_ptr()
        signals = rust_dataset.extract_precursor_ms1_signals(
            [coord],
            rt_window_sec=rt_window_sec,
            mz_tol_ppm=mz_tol_ppm,
            im_window=im_window,
            n_isotopes=n_isotopes,
            num_threads=4,
        )

        if len(signals) == 0:
            return None

        sig = signals[0]
        return {
            "rt_coords": np.array(sig.rt_coords),
            "rt_intensities": np.array(sig.rt_intensities),
            "im_coords": np.array(sig.im_coords),
            "im_intensities": np.array(sig.im_intensities),
            "isotope_mz": np.array(sig.isotope_mz),
            "isotope_intensity": np.array(sig.isotope_intensity),
            "rt_apex": sig.rt_moments.apex,
            "rt_fwhm": sig.rt_moments.fwhm,
            "im_apex": sig.im_moments.apex,
            "im_fwhm": sig.im_moments.fwhm,
        }

    def get_fragment_spectrum(self, precursor_id: int) -> Optional[TimsFrame]:
        """Get merged fragment spectrum for a precursor."""
        grouped = self.precursors_df.groupby('precursor_id').agg({
            'raw_data': 'sum',
        })

        if precursor_id not in grouped.index:
            return None

        return grouped.loc[precursor_id, 'raw_data']

    def visualize(
        self,
        precursor_id: int,
        save_path: Optional[Path] = None,
        rt_window_sec: float = 30.0,
        show: bool = True,
    ) -> plt.Figure:
        """
        Create comprehensive visualization for a precursor.

        Args:
            precursor_id: Precursor ID to visualize
            save_path: Optional path to save figure
            rt_window_sec: RT window for MS1 extraction
            show: Whether to display the plot

        Returns:
            matplotlib Figure object
        """
        if precursor_id not in self.precursor_lookup:
            raise ValueError(f"Precursor {precursor_id} not found")

        meta = self.precursor_lookup[precursor_id]

        # Get identification
        idents = self.get_precursor_identification(precursor_id)

        # Extract MS1 signal
        ms1_signal = self.extract_ms1_signal(
            precursor_id,
            rt_window_sec=rt_window_sec,
        )

        # Get fragment spectrum
        fragment = self.get_fragment_spectrum(precursor_id)

        # Create figure
        fig = plt.figure(figsize=(16, 12))
        gs = GridSpec(3, 3, figure=fig, height_ratios=[1, 1, 1.5])

        # Title with identification info
        mz = meta.get("mono_mz") or meta.get("largest_peak_mz", 0)
        charge = meta.get("charge") or "?"

        title_parts = [f"Precursor {precursor_id}"]
        title_parts.append(f"m/z: {mz:.4f}")
        title_parts.append(f"z: {charge}")

        # Add identification info
        for engine, ident in idents.items():
            if ident:
                peptide = ident.get("peptide", "")
                if peptide:
                    title_parts.append(f"\n{engine}: {peptide}")

        fig.suptitle(" | ".join(title_parts[:3]) + ("".join(title_parts[3:]) if len(title_parts) > 3 else ""),
                     fontsize=12, fontweight="bold")

        # Panel 1: XIC (Chromatogram)
        ax1 = fig.add_subplot(gs[0, 0])
        if ms1_signal and len(ms1_signal["rt_coords"]) > 0:
            rt_min = ms1_signal["rt_coords"] / 60.0
            ax1.fill_between(rt_min, ms1_signal["rt_intensities"], alpha=0.3, color="blue")
            ax1.plot(rt_min, ms1_signal["rt_intensities"], "b-", linewidth=1.5)

            # Mark apex
            apex_min = ms1_signal["rt_apex"] / 60.0
            ax1.axvline(apex_min, color="red", linestyle="--", alpha=0.7,
                       label=f"Apex: {apex_min:.2f} min")

            # Mark FWHM
            fwhm_min = ms1_signal["rt_fwhm"] / 60.0
            ax1.axhline(max(ms1_signal["rt_intensities"]) / 2, color="green",
                       linestyle=":", alpha=0.5, label=f"FWHM: {fwhm_min:.2f} min")

            ax1.legend(fontsize=8)
        else:
            ax1.text(0.5, 0.5, "No MS1 signal", ha="center", va="center", transform=ax1.transAxes)

        ax1.set_xlabel("Retention Time (min)")
        ax1.set_ylabel("Intensity")
        ax1.set_title("XIC (Chromatogram)")
        ax1.ticklabel_format(style='scientific', axis='y', scilimits=(0,0))

        # Panel 2: Mobilogram
        ax2 = fig.add_subplot(gs[0, 1])
        if ms1_signal and len(ms1_signal["im_coords"]) > 0:
            ax2.fill_between(ms1_signal["im_coords"], ms1_signal["im_intensities"],
                           alpha=0.3, color="green")
            ax2.plot(ms1_signal["im_coords"], ms1_signal["im_intensities"], "g-", linewidth=1.5)

            # Mark apex
            ax2.axvline(ms1_signal["im_apex"], color="red", linestyle="--", alpha=0.7,
                       label=f"Apex: {ms1_signal['im_apex']:.3f}")
            ax2.axhline(max(ms1_signal["im_intensities"]) / 2, color="orange",
                       linestyle=":", alpha=0.5, label=f"FWHM: {ms1_signal['im_fwhm']:.3f}")
            ax2.legend(fontsize=8)
        else:
            ax2.text(0.5, 0.5, "No mobility signal", ha="center", va="center", transform=ax2.transAxes)

        ax2.set_xlabel("Ion Mobility (1/K₀)")
        ax2.set_ylabel("Intensity")
        ax2.set_title("Mobilogram")
        ax2.ticklabel_format(style='scientific', axis='y', scilimits=(0,0))

        # Panel 3: Isotope envelope
        ax3 = fig.add_subplot(gs[0, 2])
        if ms1_signal and len(ms1_signal["isotope_mz"]) > 0:
            iso_mz = ms1_signal["isotope_mz"]
            iso_int = ms1_signal["isotope_intensity"]

            # Filter out zero intensities
            mask = iso_int > 0
            if mask.any():
                ax3.bar(range(len(iso_mz[mask])), iso_int[mask], width=0.6, color="purple", alpha=0.7)
                ax3.set_xticks(range(len(iso_mz[mask])))
                ax3.set_xticklabels([f"M+{i}\n{iso_mz[mask][i]:.3f}" for i in range(len(iso_mz[mask]))],
                                   fontsize=8)
            else:
                ax3.text(0.5, 0.5, "No isotope signal", ha="center", va="center", transform=ax3.transAxes)
        else:
            ax3.text(0.5, 0.5, "No isotope data", ha="center", va="center", transform=ax3.transAxes)

        ax3.set_xlabel("Isotope")
        ax3.set_ylabel("Intensity")
        ax3.set_title("Isotope Envelope")
        ax3.ticklabel_format(style='scientific', axis='y', scilimits=(0,0))

        # Panel 4-5: Search engine comparison
        ax_engines = fig.add_subplot(gs[1, :2])
        engine_data = []
        colors = {"fragpipe": "blue", "diann": "green", "sage": "orange"}

        for engine, ident in idents.items():
            if ident:
                engine_data.append({
                    "engine": engine.upper(),
                    "peptide": ident.get("peptide", "N/A"),
                    "protein": ident.get("protein", "N/A")[:40],
                    "qvalue": ident.get("qvalue", "N/A"),
                    "score": ident.get("hyperscore") or ident.get("probability") or "N/A",
                })

        if engine_data:
            # Create table
            cell_text = [[d["engine"], d["peptide"], d["protein"],
                         f"{d['qvalue']:.4f}" if isinstance(d['qvalue'], float) else str(d['qvalue']),
                         f"{d['score']:.2f}" if isinstance(d['score'], float) else str(d['score'])]
                        for d in engine_data]

            table = ax_engines.table(
                cellText=cell_text,
                colLabels=["Engine", "Peptide", "Protein", "Q-value", "Score"],
                loc="center",
                cellLoc="left",
            )
            table.auto_set_font_size(False)
            table.set_fontsize(9)
            table.scale(1.2, 1.5)

            # Color code by engine
            for i, d in enumerate(engine_data):
                for j in range(5):
                    table[(i+1, j)].set_facecolor(colors.get(d["engine"].lower(), "gray") + "20")
        else:
            ax_engines.text(0.5, 0.5, "No search engine identifications",
                          ha="center", va="center", fontsize=12)

        ax_engines.axis("off")
        ax_engines.set_title("Search Engine Results", fontsize=11, fontweight="bold")

        # Panel 6: Precursor metadata
        ax_meta = fig.add_subplot(gs[1, 2])
        meta_text = [
            f"Precursor ID: {precursor_id}",
            f"Frame ID: {meta.get('frame_id', 'N/A')}",
            f"m/z (mono): {meta.get('mono_mz', 'N/A')}",
            f"m/z (largest): {meta.get('largest_peak_mz', 'N/A'):.4f}",
            f"Charge: {meta.get('charge', 'N/A')}",
            f"Mobility: {meta.get('inverse_ion_mobility', 'N/A'):.4f}",
            f"Intensity: {meta.get('precuror_total_intensity', 'N/A'):.0f}",
            f"Isolation m/z: {meta.get('isolation_mz', 'N/A'):.4f}",
            f"Isolation width: {meta.get('isolation_width', 'N/A'):.1f}",
        ]
        ax_meta.text(0.05, 0.95, "\n".join(meta_text), transform=ax_meta.transAxes,
                    fontsize=10, verticalalignment="top", fontfamily="monospace",
                    bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.5))
        ax_meta.axis("off")
        ax_meta.set_title("Precursor Metadata", fontsize=11, fontweight="bold")

        # Panel 7: Fragment spectrum (full width)
        ax_frag = fig.add_subplot(gs[2, :])
        if fragment is not None and len(fragment.mz) > 0:
            mz_vals = fragment.mz
            intensities = fragment.intensity

            # Normalize intensities
            max_int = max(intensities) if len(intensities) > 0 else 1
            norm_int = intensities / max_int * 100

            # Plot as stem
            markerline, stemlines, baseline = ax_frag.stem(
                mz_vals, norm_int, linefmt="b-", markerfmt="", basefmt="k-"
            )
            plt.setp(stemlines, linewidth=0.5)

            # Annotate top peaks
            top_indices = np.argsort(intensities)[-10:]
            for idx in top_indices:
                if norm_int[idx] > 10:
                    ax_frag.annotate(f"{mz_vals[idx]:.2f}",
                                    xy=(mz_vals[idx], norm_int[idx]),
                                    xytext=(0, 5), textcoords="offset points",
                                    fontsize=7, ha="center", rotation=45)

            ax_frag.set_xlim(min(mz_vals) - 50, max(mz_vals) + 50)
        else:
            ax_frag.text(0.5, 0.5, "No fragment spectrum", ha="center", va="center",
                        transform=ax_frag.transAxes, fontsize=12)

        ax_frag.set_xlabel("m/z")
        ax_frag.set_ylabel("Relative Intensity (%)")
        ax_frag.set_title(f"Fragment Spectrum (MS2) - {len(fragment.mz) if fragment else 0} peaks",
                         fontsize=11, fontweight="bold")

        plt.tight_layout()

        if save_path:
            plt.savefig(save_path, dpi=150, bbox_inches="tight")
            print(f"Saved figure to: {save_path}")

        if show:
            plt.show()

        return fig

    def list_identified_precursors(self, min_probability: float = 0.9) -> pd.DataFrame:
        """List precursors with high-confidence identifications."""
        if self.psm_df is None:
            return pd.DataFrame()

        # Filter by probability
        good_psms = self.psm_df[self.psm_df["Probability"] >= min_probability].copy()

        # Add metadata
        def get_meta(pid):
            if pid in self.precursor_lookup:
                m = self.precursor_lookup[pid]
                return pd.Series({
                    "mono_mz": m.get("mono_mz"),
                    "charge": m.get("charge"),
                    "mobility": m.get("inverse_ion_mobility"),
                })
            return pd.Series({"mono_mz": None, "charge": None, "mobility": None})

        meta_df = good_psms["precursor_id"].apply(get_meta)
        good_psms = pd.concat([good_psms, meta_df], axis=1)

        return good_psms[[
            "precursor_id", "Modified Peptide", "Peptide", "Charge", "Protein",
            "Probability", "mono_mz", "mobility"
        ]].sort_values("Probability", ascending=False)


def load_unified_index(index_path: Path, raw_file_filter: Optional[str] = None) -> pd.DataFrame:
    """Load the unified precursor index."""
    df = pd.read_parquet(index_path)
    if raw_file_filter:
        df = df[df["raw_file"].str.contains(raw_file_filter, case=False, na=False)]
    return df


def find_disagreements(index_df: pd.DataFrame) -> pd.DataFrame:
    """Find precursors where search engines disagree."""
    # Only look at precursors with 2+ engines
    multi = index_df[index_df["n_engines"] >= 2].copy()

    # Check if peptides match (after normalization)
    def check_agreement(row):
        peptides = []
        if pd.notna(row.get("fragpipe_peptide")):
            peptides.append(row["fragpipe_peptide"].replace("I", "L").upper())
        if pd.notna(row.get("sage_peptide")):
            peptides.append(row["sage_peptide"].replace("I", "L").upper())
        if pd.notna(row.get("diann_peptide")):
            peptides.append(row["diann_peptide"].replace("I", "L").upper())
        return len(set(peptides)) > 1

    multi["is_disagreement"] = multi.apply(check_agreement, axis=1)
    return multi[multi["is_disagreement"]]


def main():
    parser = argparse.ArgumentParser(
        description="Visualize precursor with search engine identifications",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Visualize using unified index (recommended)
  python visualize_precursor.py --index data/processed/PXD019086/precursor_index.parquet \\
      --raw /path/to/data.d --precursor-id 68381

  # List top identified precursors
  python visualize_precursor.py --index data/processed/PXD019086/precursor_index.parquet \\
      --raw /path/to/data.d --list-top 20

  # Interactive mode (pick from list)
  python visualize_precursor.py --index data/processed/PXD019086/precursor_index.parquet \\
      --raw /path/to/data.d --interactive

  # Show disagreements between engines
  python visualize_precursor.py --index data/processed/PXD019086/precursor_index.parquet \\
      --raw /path/to/data.d --show-disagreements
"""
    )

    parser.add_argument("raw_data", type=Path, nargs="?", help="Path to .d folder")
    parser.add_argument("psm_file", type=Path, nargs="?", help="Path to FragPipe psm.tsv")
    parser.add_argument("--index", type=Path, help="Path to unified precursor_index.parquet")
    parser.add_argument("--raw", type=Path, help="Path to .d folder (alternative to positional)")
    parser.add_argument("--precursor-id", "-p", type=int, help="Precursor ID to visualize")
    parser.add_argument("--diann", type=Path, help="Path to DIA-NN report.parquet")
    parser.add_argument("--sage", type=Path, help="Path to Sage results.sage.parquet")
    parser.add_argument("--extracted", type=Path, help="Path to extracted data directory")
    parser.add_argument("--list-top", type=int, help="List top N identified precursors")
    parser.add_argument("--show-disagreements", action="store_true", help="Show precursors where engines disagree")
    parser.add_argument("--interactive", "-i", action="store_true", help="Interactive mode")
    parser.add_argument("--save", type=Path, help="Save figure to file")
    parser.add_argument("--rt-window", type=float, default=30.0, help="RT window in seconds")
    parser.add_argument("--no-bruker-sdk", action="store_true", help="Don't use Bruker SDK")

    args = parser.parse_args()

    # Resolve raw data path
    raw_data_path = args.raw or args.raw_data

    # Load unified index if provided
    index_df = None
    if args.index and args.index.exists():
        print(f"Loading unified index: {args.index}")
        raw_file_name = raw_data_path.stem if raw_data_path else None
        index_df = load_unified_index(args.index, raw_file_name)
        print(f"  Loaded {len(index_df)} precursors")

        if args.show_disagreements:
            disagreements = find_disagreements(index_df)
            print(f"\nFound {len(disagreements)} precursors with engine disagreements:")
            print("=" * 100)
            cols = ["precursor_id", "fragpipe_peptide", "sage_peptide", "fragpipe_probability", "sage_qvalue"]
            cols = [c for c in cols if c in disagreements.columns]
            print(disagreements[cols].head(30).to_string())
            return

    # Check we have enough args
    if not raw_data_path:
        parser.print_help()
        print("\nPlease specify raw data path with --raw or as positional argument")
        return

    # Initialize visualizer
    viz = PrecursorVisualizer(
        raw_data_path=raw_data_path,
        psm_path=args.psm_file,
        diann_path=args.diann,
        sage_path=args.sage,
        extracted_path=args.extracted,
        use_bruker_sdk=not args.no_bruker_sdk,
    )

    # Inject unified index data
    if index_df is not None:
        viz.index_df = index_df

    # List mode
    if args.list_top:
        if index_df is not None:
            print(f"\nTop {args.list_top} precursors from unified index:")
            print("=" * 100)
            # Sort by fragpipe_probability or n_engines
            sorted_df = index_df.sort_values(
                ["n_engines", "fragpipe_probability"],
                ascending=[False, False]
            )
            cols = ["precursor_id", "consensus_peptide", "n_engines",
                    "fragpipe_peptide", "sage_peptide", "fragpipe_probability"]
            cols = [c for c in cols if c in sorted_df.columns]
            print(sorted_df[cols].head(args.list_top).to_string())
        else:
            print(f"\nTop {args.list_top} identified precursors:")
            print("=" * 80)
            top_df = viz.list_identified_precursors()
            if len(top_df) > 0:
                print(top_df.head(args.list_top).to_string())
            else:
                print("No identified precursors found.")
        return

    # Interactive mode
    if args.interactive:
        if index_df is not None:
            print("\nPrecursors from unified index:")
            sorted_df = index_df.sort_values(
                ["n_engines", "fragpipe_probability"],
                ascending=[False, False]
            )
            cols = ["precursor_id", "consensus_peptide", "n_engines", "fragpipe_probability"]
            cols = [c for c in cols if c in sorted_df.columns]
            print(sorted_df[cols].head(30).to_string())
        else:
            print("\nListing top identified precursors...")
            top_df = viz.list_identified_precursors()
            if len(top_df) == 0:
                print("No identified precursors found.")
                return
            print("\nTop 20 precursors:")
            print(top_df.head(20).to_string())

        while True:
            try:
                pid = input("\nEnter precursor ID to visualize (or 'q' to quit): ")
                if pid.lower() == 'q':
                    break

                pid = int(pid)
                viz.visualize(pid, rt_window_sec=args.rt_window)
            except ValueError:
                print("Invalid input. Please enter a number or 'q'.")
            except KeyboardInterrupt:
                break
        return

    # Single precursor mode
    if args.precursor_id:
        viz.visualize(
            args.precursor_id,
            save_path=args.save,
            rt_window_sec=args.rt_window,
        )
    else:
        parser.print_help()
        print("\nPlease specify --precursor-id, --list-top, --show-disagreements, or --interactive")


if __name__ == "__main__":
    main()
