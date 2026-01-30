#!/usr/bin/env python3
"""
Spectrum Annotation Wrapper

Provides batch fragment annotation using fragment_matching.py, with output
formatted as parallel arrays suitable for Parquet storage.

The annotation produces parallel arrays where index i corresponds to peak i:
- fragment_ion_type[i]: "b", "y", or "other" (unmatched)
- fragment_ion_number[i]: Ion position (1-indexed for b/y, 0 for unmatched)
- fragment_ion_charge[i]: Fragment charge state
- fragment_theoretical_mz[i]: Theoretical m/z (0 for unmatched)
- fragment_error_ppm[i]: Mass error in ppm (0 for unmatched)

Usage:
    from annotate_spectrum import SpectrumAnnotator

    annotator = SpectrumAnnotator(mz_tolerance_ppm=20.0)

    result = annotator.annotate(
        sequence="PEPTIDEK",
        charge=2,
        fragment_mz=np.array([...]),
        fragment_intensity=np.array([...]),
    )

    # Parallel arrays for Parquet
    print(result.ion_type)       # ["b", "other", "y", ...]
    print(result.ion_number)     # [3, 0, 2, ...]
    print(result.intensity_explained)  # 0.75
"""

import sys
from dataclasses import dataclass, field
from typing import List, Optional, Tuple
from pathlib import Path

import numpy as np
import pandas as pd

# Import fragment matching
from fragment_matching import (
    FragmentMatcher,
    MatchConfig,
    FragmentMatchResult,
    sage_to_imspy_sequence,
    IMSPY_AVAILABLE,
)


@dataclass
class AnnotatedSpectrum:
    """
    Annotated spectrum with parallel arrays for each peak.

    All arrays have the same length as the input spectrum.
    """
    # Parallel arrays (same length as spectrum)
    ion_type: List[str]              # "b", "y", "other" per peak
    ion_number: np.ndarray           # Position (0 = unmatched)
    ion_charge: np.ndarray           # Fragment charge (0 = unmatched)
    theoretical_mz: np.ndarray       # Theoretical m/z (0 = unmatched)
    error_ppm: np.ndarray            # Mass error (0 = unmatched)

    # Summary metrics
    is_annotated: bool               # Has sequence for annotation
    n_peaks: int                     # Total peaks in spectrum
    n_matched_peaks: int             # Peaks matched to theoretical
    intensity_explained: float       # Fraction of intensity matched
    sequence_coverage_b: float       # Fraction of b ions found
    sequence_coverage_y: float       # Fraction of y ions found

    # Optional: keep the original match result for debugging
    match_result: Optional[FragmentMatchResult] = None

    def to_dict(self) -> dict:
        """Convert to dictionary for Parquet serialization."""
        return {
            'fragment_ion_type': self.ion_type,
            'fragment_ion_number': self.ion_number.tolist(),
            'fragment_ion_charge': self.ion_charge.tolist(),
            'fragment_theoretical_mz': self.theoretical_mz.tolist(),
            'fragment_error_ppm': self.error_ppm.tolist(),
            'is_annotated': self.is_annotated,
            'n_matched_peaks': self.n_matched_peaks,
            'intensity_explained': self.intensity_explained,
            'sequence_coverage_b': self.sequence_coverage_b,
            'sequence_coverage_y': self.sequence_coverage_y,
        }


class SpectrumAnnotator:
    """
    Annotates spectra with b/y ion identifications.

    Uses FragmentMatcher from fragment_matching.py for theoretical fragment
    generation and matching, then reformats output as parallel arrays.
    """

    def __init__(
        self,
        mz_tolerance_ppm: float = 20.0,
        max_fragment_charge: int = 2,
        min_intensity: float = 0.0,
    ):
        """
        Initialize annotator.

        Args:
            mz_tolerance_ppm: Mass tolerance for fragment matching
            max_fragment_charge: Maximum fragment charge state to consider
            min_intensity: Minimum intensity for a peak to be considered
        """
        if not IMSPY_AVAILABLE:
            raise ImportError(
                "imspy-core is required for spectrum annotation. "
                "Ensure it's installed and on the Python path."
            )

        self.config = MatchConfig(
            mz_tolerance_ppm=mz_tolerance_ppm,
            max_fragment_charge=max_fragment_charge,
            min_intensity=min_intensity,
            ion_types=['b', 'y'],
        )
        self.matcher = FragmentMatcher(self.config)

    def annotate(
        self,
        sequence: str,
        charge: int,
        fragment_mz: np.ndarray,
        fragment_intensity: np.ndarray,
    ) -> AnnotatedSpectrum:
        """
        Annotate a single spectrum with b/y ion identifications.

        Args:
            sequence: Peptide sequence (Sage format with [+mass] mods supported)
            charge: Precursor charge state
            fragment_mz: Array of fragment m/z values
            fragment_intensity: Array of fragment intensities

        Returns:
            AnnotatedSpectrum with parallel arrays and summary metrics
        """
        fragment_mz = np.asarray(fragment_mz, dtype=np.float64)
        fragment_intensity = np.asarray(fragment_intensity, dtype=np.float64)
        n_peaks = len(fragment_mz)

        # Handle empty spectrum
        if n_peaks == 0:
            return self._empty_annotation(is_annotated=bool(sequence))

        # Handle missing sequence
        if not sequence or not isinstance(sequence, str):
            return self._unannotated_spectrum(fragment_mz, fragment_intensity)

        # Match using FragmentMatcher
        try:
            match_result = self.matcher.match_spectrum(
                sequence=sequence,
                precursor_charge=charge,
                experimental_mz=fragment_mz,
                experimental_intensity=fragment_intensity,
            )
        except Exception as e:
            # If matching fails (e.g., invalid sequence), return unannotated
            print(f"Warning: Annotation failed for {sequence}: {e}")
            return self._unannotated_spectrum(fragment_mz, fragment_intensity)

        # Build parallel arrays from matches
        ion_type = ["other"] * n_peaks
        ion_number = np.zeros(n_peaks, dtype=np.int32)
        ion_charge = np.zeros(n_peaks, dtype=np.int32)
        theoretical_mz = np.zeros(n_peaks, dtype=np.float64)
        error_ppm = np.zeros(n_peaks, dtype=np.float64)

        # Create lookup from experimental m/z to peak index
        # (FragmentMatcher returns matched experimental m/z values)
        mz_to_idx = {mz: i for i, mz in enumerate(fragment_mz)}

        for match in match_result.matches:
            # Find the peak index for this match
            # Use closest m/z since floating point comparison may differ
            exp_mz = match.experimental_mz
            best_idx = None
            best_diff = float('inf')

            for mz, idx in mz_to_idx.items():
                diff = abs(mz - exp_mz)
                if diff < best_diff:
                    best_diff = diff
                    best_idx = idx

            if best_idx is not None and best_diff < 0.001:  # Within 1 mDa
                ion_type[best_idx] = match.ion_type
                ion_number[best_idx] = match.ion_number
                ion_charge[best_idx] = match.charge
                theoretical_mz[best_idx] = match.theoretical_mz
                error_ppm[best_idx] = match.mass_error_ppm

        return AnnotatedSpectrum(
            ion_type=ion_type,
            ion_number=ion_number,
            ion_charge=ion_charge,
            theoretical_mz=theoretical_mz,
            error_ppm=error_ppm,
            is_annotated=True,
            n_peaks=n_peaks,
            n_matched_peaks=match_result.n_matched,
            intensity_explained=match_result.intensity_explained,
            sequence_coverage_b=match_result.coverage_b,
            sequence_coverage_y=match_result.coverage_y,
            match_result=match_result,
        )

    def annotate_batch(
        self,
        sequences: List[str],
        charges: List[int],
        spectra: List[Tuple[np.ndarray, np.ndarray]],
    ) -> List[AnnotatedSpectrum]:
        """
        Annotate a batch of spectra.

        Args:
            sequences: List of peptide sequences
            charges: List of precursor charges
            spectra: List of (mz_array, intensity_array) tuples

        Returns:
            List of AnnotatedSpectrum objects
        """
        results = []
        for seq, charge, (mz, intensity) in zip(sequences, charges, spectra):
            result = self.annotate(seq, charge, mz, intensity)
            results.append(result)
        return results

    def annotate_dataframe(
        self,
        df: pd.DataFrame,
        sequence_col: str = 'sage_peptide',
        charge_col: str = 'charge',
        mz_col: str = 'fragment_mz',
        intensity_col: str = 'fragment_intensity',
    ) -> pd.DataFrame:
        """
        Annotate spectra from a DataFrame and add annotation columns.

        Args:
            df: DataFrame with spectrum data
            sequence_col: Column containing peptide sequences
            charge_col: Column containing charges
            mz_col: Column containing fragment m/z arrays
            intensity_col: Column containing fragment intensity arrays

        Returns:
            DataFrame with added annotation columns
        """
        # Prepare output columns
        annotations = []

        for _, row in df.iterrows():
            seq = row.get(sequence_col)
            charge = int(row.get(charge_col, 2))
            mz = row.get(mz_col)
            intensity = row.get(intensity_col)

            # Handle missing/empty arrays
            if mz is None or len(mz) == 0:
                result = self._empty_annotation(is_annotated=bool(seq))
            else:
                result = self.annotate(
                    sequence=seq if pd.notna(seq) else None,
                    charge=charge,
                    fragment_mz=np.array(mz),
                    fragment_intensity=np.array(intensity),
                )

            annotations.append(result.to_dict())

        # Add columns to DataFrame
        result_df = df.copy()
        for key in annotations[0].keys():
            result_df[key] = [a[key] for a in annotations]

        return result_df

    def _empty_annotation(self, is_annotated: bool = False) -> AnnotatedSpectrum:
        """Create annotation for empty spectrum."""
        return AnnotatedSpectrum(
            ion_type=[],
            ion_number=np.array([], dtype=np.int32),
            ion_charge=np.array([], dtype=np.int32),
            theoretical_mz=np.array([], dtype=np.float64),
            error_ppm=np.array([], dtype=np.float64),
            is_annotated=is_annotated,
            n_peaks=0,
            n_matched_peaks=0,
            intensity_explained=0.0,
            sequence_coverage_b=0.0,
            sequence_coverage_y=0.0,
        )

    def _unannotated_spectrum(
        self,
        fragment_mz: np.ndarray,
        fragment_intensity: np.ndarray,
    ) -> AnnotatedSpectrum:
        """Create annotation for spectrum without sequence (all peaks unmatched)."""
        n_peaks = len(fragment_mz)
        return AnnotatedSpectrum(
            ion_type=["other"] * n_peaks,
            ion_number=np.zeros(n_peaks, dtype=np.int32),
            ion_charge=np.zeros(n_peaks, dtype=np.int32),
            theoretical_mz=np.zeros(n_peaks, dtype=np.float64),
            error_ppm=np.zeros(n_peaks, dtype=np.float64),
            is_annotated=False,
            n_peaks=n_peaks,
            n_matched_peaks=0,
            intensity_explained=0.0,
            sequence_coverage_b=0.0,
            sequence_coverage_y=0.0,
        )


# ============================================================================
# Demo / Testing
# ============================================================================

def demo():
    """Demonstrate spectrum annotation."""
    print("Spectrum Annotation Demo")
    print("=" * 60)

    if not IMSPY_AVAILABLE:
        print("ERROR: imspy-core is not available. Cannot run demo.")
        return

    # Create annotator
    annotator = SpectrumAnnotator(mz_tolerance_ppm=20.0)

    # Example peptide
    sequence = "PEPTIDEK"
    charge = 2

    print(f"\nPeptide: {sequence}")
    print(f"Charge: {charge}+")

    # Generate theoretical fragments for demo
    theoretical = annotator.matcher.generate_theoretical_fragments(sequence, charge)

    # Create simulated spectrum with some b and y ions
    np.random.seed(42)

    # Take some theoretical peaks
    b_mz = [mz for _, _, mz, _ in theoretical['b'][:4]]
    y_mz = [mz for _, _, mz, _ in theoretical['y'][:4]]

    # Add small mass errors
    fragment_mz = np.array(b_mz + y_mz)
    fragment_mz = fragment_mz * (1 + np.random.uniform(-10, 10, len(fragment_mz)) * 1e-6)

    # Add noise peaks
    noise_mz = np.random.uniform(200, 800, 10)
    all_mz = np.concatenate([fragment_mz, noise_mz])

    # Intensities
    signal_int = np.random.uniform(500, 1000, len(fragment_mz))
    noise_int = np.random.uniform(10, 100, len(noise_mz))
    all_int = np.concatenate([signal_int, noise_int])

    # Sort by m/z
    order = np.argsort(all_mz)
    all_mz = all_mz[order]
    all_int = all_int[order]

    print(f"\nSpectrum: {len(all_mz)} peaks ({len(fragment_mz)} signal, {len(noise_mz)} noise)")

    # Annotate
    result = annotator.annotate(
        sequence=sequence,
        charge=charge,
        fragment_mz=all_mz,
        fragment_intensity=all_int,
    )

    print("\nAnnotation Results:")
    print("-" * 40)
    print(f"Annotated: {result.is_annotated}")
    print(f"Peaks: {result.n_peaks}")
    print(f"Matched: {result.n_matched_peaks}")
    print(f"Intensity explained: {result.intensity_explained:.1%}")
    print(f"Coverage b: {result.sequence_coverage_b:.1%}")
    print(f"Coverage y: {result.sequence_coverage_y:.1%}")

    print("\nParallel arrays (first 10 peaks):")
    print(f"  ion_type:    {result.ion_type[:10]}")
    print(f"  ion_number:  {result.ion_number[:10].tolist()}")
    print(f"  ion_charge:  {result.ion_charge[:10].tolist()}")
    print(f"  theo_mz:     {np.round(result.theoretical_mz[:10], 2).tolist()}")
    print(f"  error_ppm:   {np.round(result.error_ppm[:10], 1).tolist()}")

    # Test unannotated spectrum
    print("\n" + "=" * 60)
    print("Unannotated spectrum (no sequence):")
    unannotated = annotator.annotate(
        sequence=None,
        charge=2,
        fragment_mz=all_mz,
        fragment_intensity=all_int,
    )
    print(f"  is_annotated: {unannotated.is_annotated}")
    print(f"  n_matched_peaks: {unannotated.n_matched_peaks}")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Spectrum annotation")
    parser.add_argument("--demo", action="store_true", help="Run demo")

    args = parser.parse_args()

    if args.demo:
        demo()
    else:
        parser.print_help()
