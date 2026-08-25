#!/usr/bin/env python3
"""
Fragment Ion Matching using imspy

Generates theoretical fragment ions from peptide sequences and matches them
against experimental spectra. Provides search-engine-independent fragment
annotation.

Uses imspy-core's PeptideSequence for theoretical fragment generation.
"""

import os
import sys
from dataclasses import dataclass, field
from typing import List, Dict, Tuple, Optional, Union
from pathlib import Path

import re
import numpy as np
import pandas as pd

# Add imspy-core to path
# Local checkout of https://github.com/theGreatHerrLebert/rustims.
# Override with RUSTIMS_ROOT; defaults to a sibling of this repository.
RUSTIMS_ROOT = Path(
    os.environ.get("RUSTIMS_ROOT", Path(__file__).resolve().parents[1] / ".." / "rustims")
).expanduser().resolve()
IMSPY_CORE_PATH = RUSTIMS_ROOT / "packages" / "imspy-core" / "src"
if str(IMSPY_CORE_PATH) not in sys.path:
    sys.path.insert(0, str(IMSPY_CORE_PATH))

try:
    from imspy_core.data.peptide import PeptideSequence, PeptideProductIon
    from imspy_core.data.spectrum import MzSpectrum
    IMSPY_AVAILABLE = True
except ImportError as e:
    print(f"Warning: imspy-core not available: {e}")
    IMSPY_AVAILABLE = False


# ============================================================================
# Sequence Format Conversion
# ============================================================================

# Common mass -> UNIMOD mapping (for converting Sage [+mass] to [UNIMOD:ID])
MASS_TO_UNIMOD = {
    57.021465: 4,    # Carbamidomethyl (C)
    57.021464: 4,    # Carbamidomethyl (C) - variant
    57.02146: 4,     # Carbamidomethyl (C) - rounded
    15.994915: 35,   # Oxidation (M)
    15.9949: 35,     # Oxidation (M) - rounded
    42.010565: 1,    # Acetyl (N-term)
    42.0106: 1,      # Acetyl - rounded
    -17.026549: 28,  # Gln->pyro-Glu (Q)
    -17.026548: 28,  # Gln->pyro-Glu (Q) - variant
    -18.010565: 27,  # Glu->pyro-Glu (E)
    79.966331: 21,   # Phospho (STY)
}


def mass_delta_to_unimod(mass: float, tolerance: float = 0.001) -> Optional[int]:
    """Convert a mass delta to UNIMOD ID if known."""
    for known_mass, unimod_id in MASS_TO_UNIMOD.items():
        if abs(mass - known_mass) < tolerance:
            return unimod_id
    return None


def sage_to_imspy_sequence(sage_seq: str) -> str:
    """
    Convert Sage peptide notation to imspy format.

    Sage format: PEPTIDE[+15.994915]K or C[+57.021465]
    imspy format: PEPTIDE[UNIMOD:35]K or C[UNIMOD:4]
    """
    def replace_mass(match):
        mass_str = match.group(1)
        mass = float(mass_str)

        unimod_id = mass_delta_to_unimod(mass)
        if unimod_id is not None:
            return f'[UNIMOD:{unimod_id}]'

        # If not found, keep as mass delta (imspy might support it)
        return f'[{mass_str}]'

    # Replace [+mass] or [-mass] patterns
    result = re.sub(r'\[([+-]?\d+\.?\d*)\]', replace_mass, sage_seq)
    return result


def fragpipe_to_imspy_sequence(fp_seq: str) -> str:
    """
    Convert FragPipe modified peptide notation to imspy format.

    FragPipe format: n[+42.0106]PEPTIDEK or PEPTIDEK[+15.9949]
    imspy format: [UNIMOD:1]PEPTIDEK or PEPTIDEK[UNIMOD:35]
    """
    # Similar conversion, but also handle n-term notation
    result = fp_seq

    # Handle n-terminal mods: n[+42.0106] -> [UNIMOD:1]
    result = re.sub(r'n\[([+-]?\d+\.?\d*)\]', lambda m: f'[UNIMOD:{mass_delta_to_unimod(float(m.group(1))) or m.group(1)}]', result)

    # Handle regular mods
    def replace_mass(match):
        mass_str = match.group(1)
        mass = float(mass_str)

        unimod_id = mass_delta_to_unimod(mass)
        if unimod_id is not None:
            return f'[UNIMOD:{unimod_id}]'

        return f'[{mass_str}]'

    result = re.sub(r'\[([+-]?\d+\.?\d*)\]', replace_mass, result)
    return result


@dataclass
class FragmentMatch:
    """A matched fragment ion."""
    ion_type: str           # 'b' or 'y'
    ion_number: int         # Position (e.g., 3 for b3)
    charge: int             # Fragment charge
    theoretical_mz: float   # Calculated m/z
    experimental_mz: float  # Observed m/z
    intensity: float        # Observed intensity
    mass_error_ppm: float   # ppm error
    sequence: str           # Fragment sequence (e.g., "PEP" for b3 of PEPTIDE)


@dataclass
class FragmentMatchResult:
    """Result of fragment matching for a single PSM."""
    peptide: str
    charge: int
    matches: List[FragmentMatch] = field(default_factory=list)

    # Summary stats
    n_theoretical_b: int = 0
    n_theoretical_y: int = 0
    n_matched_b: int = 0
    n_matched_y: int = 0

    # Intensity coverage
    matched_intensity: float = 0.0
    total_intensity: float = 0.0

    @property
    def n_matched(self) -> int:
        return len(self.matches)

    @property
    def coverage_b(self) -> float:
        return self.n_matched_b / self.n_theoretical_b if self.n_theoretical_b > 0 else 0

    @property
    def coverage_y(self) -> float:
        return self.n_matched_y / self.n_theoretical_y if self.n_theoretical_y > 0 else 0

    @property
    def intensity_explained(self) -> float:
        return self.matched_intensity / self.total_intensity if self.total_intensity > 0 else 0

    def to_dict(self) -> Dict:
        return {
            'peptide': self.peptide,
            'charge': self.charge,
            'n_theoretical_b': self.n_theoretical_b,
            'n_theoretical_y': self.n_theoretical_y,
            'n_matched_b': self.n_matched_b,
            'n_matched_y': self.n_matched_y,
            'n_matched_total': self.n_matched,
            'coverage_b': self.coverage_b,
            'coverage_y': self.coverage_y,
            'intensity_explained': self.intensity_explained,
        }


@dataclass
class MatchConfig:
    """Configuration for fragment matching."""
    # Mass tolerance for matching
    mz_tolerance_ppm: float = 20.0
    mz_tolerance_da: Optional[float] = None  # If set, overrides ppm

    # Ion types to match
    ion_types: List[str] = field(default_factory=lambda: ['b', 'y'])

    # Fragment charge states to consider (relative to precursor)
    max_fragment_charge: int = 2

    # Minimum intensity for a peak to be considered
    min_intensity: float = 0.0

    # Whether to use isotope patterns
    use_isotopes: bool = False
    isotope_tolerance_da: float = 0.5


class FragmentMatcher:
    """
    Matches theoretical fragments against experimental spectra.

    Uses imspy's PeptideSequence for theoretical fragment generation.
    """

    def __init__(self, config: Optional[MatchConfig] = None):
        if not IMSPY_AVAILABLE:
            raise ImportError("imspy-core is required for fragment matching")
        self.config = config or MatchConfig()

    def generate_theoretical_fragments(
        self,
        sequence: str,
        precursor_charge: int,
    ) -> Dict[str, List[Tuple[int, int, float, str]]]:
        """
        Generate theoretical fragment ions for a peptide.

        Args:
            sequence: Peptide sequence (with modifications in any supported format)
            precursor_charge: Precursor charge state

        Returns:
            Dict mapping ion_type to list of (ion_number, charge, mz, fragment_sequence)
            Ion numbers follow standard proteomics convention:
            - b ions: numbered from N-terminus (b1, b2, ...)
            - y ions: numbered from C-terminus (y1, y2, ...)
        """
        # Convert sequence to imspy format if needed (handles Sage [+mass] notation)
        imspy_seq = sage_to_imspy_sequence(sequence)
        pep = PeptideSequence(imspy_seq)

        fragments = {}

        # Determine fragment charge states to consider
        max_frag_charge = min(self.config.max_fragment_charge, precursor_charge - 1)
        max_frag_charge = max(1, max_frag_charge)

        for ion_type in self.config.ion_types:
            fragments[ion_type] = []

            for frag_charge in range(1, max_frag_charge + 1):
                # Get N-terminal (b/a/c) and C-terminal (y/x/z) ions
                n_ions, c_ions = pep.calculate_product_ion_series(
                    charge=frag_charge,
                    fragment_type=ion_type
                )

                # Add ions based on type with proper numbering
                if ion_type in ['b', 'a', 'c']:
                    # N-terminal ions: numbered 1, 2, 3, ... from N-terminus
                    for i, ion in enumerate(n_ions, 1):
                        fragments[ion_type].append((i, frag_charge, ion.mz, ion.sequence or ""))
                elif ion_type in ['y', 'x', 'z']:
                    # C-terminal ions: imspy returns them in reverse order
                    # imspy y[0] is the longest (almost full peptide)
                    # imspy y[-1] is the shortest (just C-terminal residue)
                    # Standard convention: y1 is shortest, yn-1 is longest
                    n_y = len(c_ions)
                    for i, ion in enumerate(c_ions, 1):
                        # Convert imspy index to standard y ion number
                        # imspy index 1 -> standard y(n_y)
                        # imspy index n_y -> standard y(1)
                        standard_y_number = n_y + 1 - i
                        fragments[ion_type].append((standard_y_number, frag_charge, ion.mz, ion.sequence or ""))

        return fragments

    def _mz_within_tolerance(self, theoretical_mz: float, experimental_mz: float) -> Tuple[bool, float]:
        """
        Check if two m/z values are within tolerance.

        Returns:
            (is_match, error_ppm)
        """
        if self.config.mz_tolerance_da is not None:
            error_da = abs(theoretical_mz - experimental_mz)
            error_ppm = error_da / theoretical_mz * 1e6
            return error_da <= self.config.mz_tolerance_da, error_ppm
        else:
            error_ppm = abs(theoretical_mz - experimental_mz) / theoretical_mz * 1e6
            return error_ppm <= self.config.mz_tolerance_ppm, error_ppm

    def match_spectrum(
        self,
        sequence: str,
        precursor_charge: int,
        experimental_mz: np.ndarray,
        experimental_intensity: np.ndarray,
    ) -> FragmentMatchResult:
        """
        Match theoretical fragments against an experimental spectrum.

        Args:
            sequence: Peptide sequence
            precursor_charge: Precursor charge state
            experimental_mz: Array of experimental m/z values
            experimental_intensity: Array of experimental intensities

        Returns:
            FragmentMatchResult with all matches
        """
        result = FragmentMatchResult(
            peptide=sequence,
            charge=precursor_charge,
        )

        # Filter by minimum intensity
        mask = experimental_intensity >= self.config.min_intensity
        exp_mz = experimental_mz[mask]
        exp_int = experimental_intensity[mask]

        result.total_intensity = float(np.sum(exp_int))

        # Generate theoretical fragments
        try:
            theoretical = self.generate_theoretical_fragments(sequence, precursor_charge)
        except Exception as e:
            print(f"Warning: Could not generate fragments for {sequence}: {e}")
            return result

        # Track which experimental peaks have been matched
        matched_exp_indices = set()

        for ion_type, ions in theoretical.items():
            if ion_type in ['b', 'a', 'c']:
                result.n_theoretical_b += len(ions)
            else:
                result.n_theoretical_y += len(ions)

            # ions is now list of (ion_number, charge, mz, sequence)
            for ion_number, ion_charge, theo_mz, ion_seq in ions:

                # Find best matching experimental peak
                best_match_idx = None
                best_error_ppm = float('inf')

                for i, (mz, intensity) in enumerate(zip(exp_mz, exp_int)):
                    if i in matched_exp_indices:
                        continue

                    is_match, error_ppm = self._mz_within_tolerance(theo_mz, mz)

                    if is_match and error_ppm < best_error_ppm:
                        best_match_idx = i
                        best_error_ppm = error_ppm

                if best_match_idx is not None:
                    matched_exp_indices.add(best_match_idx)

                    match = FragmentMatch(
                        ion_type=ion_type,
                        ion_number=ion_number,
                        charge=ion_charge,
                        theoretical_mz=theo_mz,
                        experimental_mz=float(exp_mz[best_match_idx]),
                        intensity=float(exp_int[best_match_idx]),
                        mass_error_ppm=best_error_ppm,
                        sequence=ion_seq,
                    )
                    result.matches.append(match)
                    result.matched_intensity += exp_int[best_match_idx]

                    if ion_type in ['b', 'a', 'c']:
                        result.n_matched_b += 1
                    else:
                        result.n_matched_y += 1

        return result

    def match_psm_batch(
        self,
        psm_df: pd.DataFrame,
        spectra: Dict[int, Tuple[np.ndarray, np.ndarray]],
        sequence_col: str = 'peptide',
        charge_col: str = 'charge',
        spectrum_id_col: str = 'psm_id',
    ) -> pd.DataFrame:
        """
        Match fragments for a batch of PSMs.

        Args:
            psm_df: DataFrame with PSMs
            spectra: Dict mapping spectrum_id -> (mz_array, intensity_array)
            sequence_col: Column name for peptide sequence
            charge_col: Column name for charge
            spectrum_id_col: Column name for spectrum identifier

        Returns:
            DataFrame with match results
        """
        results = []

        for _, row in psm_df.iterrows():
            spec_id = row[spectrum_id_col]

            if spec_id not in spectra:
                continue

            mz, intensity = spectra[spec_id]

            result = self.match_spectrum(
                sequence=row[sequence_col],
                precursor_charge=int(row[charge_col]),
                experimental_mz=mz,
                experimental_intensity=intensity,
            )

            result_dict = result.to_dict()
            result_dict[spectrum_id_col] = spec_id
            results.append(result_dict)

        return pd.DataFrame(results)


# Aliases for backwards compatibility
convert_sage_to_imspy_sequence = sage_to_imspy_sequence
convert_fragpipe_to_imspy_sequence = fragpipe_to_imspy_sequence


# ============================================================================
# Main / Demo
# ============================================================================

def demo_fragment_matching():
    """Demonstrate fragment matching on a simple example."""

    if not IMSPY_AVAILABLE:
        print("ERROR: imspy-core is not available. Cannot run demo.")
        return

    print("Fragment Matching Demo")
    print("=" * 60)

    # Example peptide
    sequence = "PEPTIDEK"
    charge = 2

    print(f"\nPeptide: {sequence}")
    print(f"Charge: {charge}+")

    # Create matcher
    config = MatchConfig(
        mz_tolerance_ppm=20.0,
        ion_types=['b', 'y'],
        max_fragment_charge=1,
    )
    matcher = FragmentMatcher(config)

    # Generate theoretical fragments using the matcher (handles y ion numbering correctly)
    print("\nTheoretical Fragments (standard proteomics convention):")
    print("-" * 40)

    theoretical = matcher.generate_theoretical_fragments(sequence, charge)

    print("\nb ions:")
    for ion_num, ion_charge, mz, seq in sorted(theoretical['b']):
        print(f"  b{ion_num}: m/z = {mz:.4f}, seq = {seq}")

    print("\ny ions:")
    for ion_num, ion_charge, mz, seq in sorted(theoretical['y']):
        print(f"  y{ion_num}: m/z = {mz:.4f}, seq = {seq}")

    # Simulate an experimental spectrum with some noise
    print("\n" + "=" * 60)
    print("Simulating experimental spectrum...")

    # Take some theoretical peaks and add noise
    np.random.seed(42)

    # Get m/z values from theoretical fragments (first 4 b and first 4 y)
    b_mz = sorted([mz for _, _, mz, _ in theoretical['b']])[:4]
    y_mz = sorted([mz for _, _, mz, _ in theoretical['y']])[:4]
    theoretical_mz = np.array(b_mz + y_mz)

    # Add mass error (up to 10 ppm)
    experimental_mz = theoretical_mz * (1 + np.random.uniform(-10, 10, len(theoretical_mz)) * 1e-6)

    # Add some noise peaks
    noise_mz = np.random.uniform(200, 800, 20)
    noise_int = np.random.uniform(10, 100, 20)

    all_mz = np.concatenate([experimental_mz, noise_mz])
    all_int = np.concatenate([np.random.uniform(500, 1000, len(experimental_mz)), noise_int])

    # Sort by m/z
    order = np.argsort(all_mz)
    all_mz = all_mz[order]
    all_int = all_int[order]

    print(f"Spectrum: {len(all_mz)} peaks ({len(experimental_mz)} signal, {len(noise_mz)} noise)")

    # Match
    result = matcher.match_spectrum(
        sequence=sequence,
        precursor_charge=charge,
        experimental_mz=all_mz,
        experimental_intensity=all_int,
    )

    print("\nMatching Results:")
    print("-" * 40)
    print(f"Matched: {result.n_matched} / {result.n_theoretical_b + result.n_theoretical_y} theoretical")
    print(f"  b ions: {result.n_matched_b} / {result.n_theoretical_b} ({result.coverage_b:.1%})")
    print(f"  y ions: {result.n_matched_y} / {result.n_theoretical_y} ({result.coverage_y:.1%})")
    print(f"Intensity explained: {result.intensity_explained:.1%}")

    print("\nIndividual matches:")
    for m in result.matches[:10]:
        print(f"  {m.ion_type}{m.ion_number}+{m.charge}: theo={m.theoretical_mz:.4f}, "
              f"exp={m.experimental_mz:.4f}, err={m.mass_error_ppm:.1f}ppm, int={m.intensity:.0f}")


def match_sage_results(
    sage_parquet: str,
    matched_fragments_parquet: str,
    output_parquet: str,
    mz_tolerance_ppm: float = 20.0,
):
    """
    Re-match Sage PSMs using imspy and compare with Sage's matching.

    Args:
        sage_parquet: Path to results.sage.parquet
        matched_fragments_parquet: Path to matched_fragments.sage.parquet
        output_parquet: Path for output comparison
        mz_tolerance_ppm: Matching tolerance
    """
    import pyarrow.parquet as pq

    print(f"Loading Sage results from {sage_parquet}...")
    sage_df = pq.read_table(sage_parquet).to_pandas()

    print(f"Loading Sage matched fragments from {matched_fragments_parquet}...")
    sage_frags = pq.read_table(matched_fragments_parquet).to_pandas()

    print(f"\nSage PSMs: {len(sage_df):,}")
    print(f"Sage fragment matches: {len(sage_frags):,}")

    # Configure matcher
    config = MatchConfig(
        mz_tolerance_ppm=mz_tolerance_ppm,
        ion_types=['b', 'y'],
        max_fragment_charge=2,
    )
    matcher = FragmentMatcher(config)

    # For each PSM, compare Sage's fragments with imspy theoretical
    # This validates that the fragment generation is consistent

    print("\nComparing fragment generation (sample)...")

    sample_psms = sage_df.head(10)

    for _, psm in sample_psms.iterrows():
        psm_id = psm['psm_id']
        peptide = psm['peptide']
        charge = psm['charge']

        # Get Sage's matched fragments for this PSM
        sage_matches = sage_frags[sage_frags['psm_id'] == psm_id]

        # Generate theoretical with imspy
        try:
            seq = convert_sage_to_imspy_sequence(peptide)
            theoretical = matcher.generate_theoretical_fragments(seq, charge)

            n_theo = sum(len(ions) for ions in theoretical.values())
            n_sage = len(sage_matches)

            print(f"\nPSM {psm_id}: {peptide}")
            print(f"  imspy theoretical: {n_theo} fragments")
            print(f"  Sage matched: {n_sage} fragments")

        except Exception as e:
            print(f"\nPSM {psm_id}: {peptide} - ERROR: {e}")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Fragment ion matching using imspy")
    parser.add_argument("--demo", action="store_true", help="Run demo")
    parser.add_argument("--sage", type=str, help="Path to Sage results.sage.parquet")
    parser.add_argument("--sage-fragments", type=str, help="Path to Sage matched_fragments.sage.parquet")
    parser.add_argument("--output", type=str, help="Output parquet path")
    parser.add_argument("--tolerance", type=float, default=20.0, help="m/z tolerance in ppm")

    args = parser.parse_args()

    if args.demo:
        demo_fragment_matching()
    elif args.sage and args.sage_fragments:
        match_sage_results(
            sage_parquet=args.sage,
            matched_fragments_parquet=args.sage_fragments,
            output_parquet=args.output or "imspy_fragments.parquet",
            mz_tolerance_ppm=args.tolerance,
        )
    else:
        parser.print_help()
