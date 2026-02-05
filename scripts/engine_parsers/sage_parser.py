"""
Sage Results Parser

Parses Sage results.sage.parquet files into standardized format.

Key features:
- m/z CALCULATION: Sage reports experimental mass, we calculate m/z from mass and charge
  Formula: mz = (expmass + charge * 1.007276) / charge
  This is intentional and documented behavior.
- RT CONVERSION: Sage reports RT in MINUTES, we convert to SECONDS at parse time
- scannr from mzML conversion is NOT the timsTOF precursor_id
- Standardizes [+mass] format to [UNIMOD:X]
- Filters decoy hits
"""

from pathlib import Path
from typing import Optional
import pandas as pd
import numpy as np

import sys
sys.path.insert(0, str(Path(__file__).parent.parent))
from sequence_utils import standardize_sage_sequence

from .base import BaseParser


# Proton mass for m/z calculation
PROTON_MASS = 1.007276


class SageParser(BaseParser):
    """Parser for Sage results.sage.parquet results.

    Sage uses mzML files converted from timsTOF .d files. The scannr field
    comes from mzML conversion and does NOT correspond to timsTOF precursor_id.
    Matching requires sequence or coordinate matching.
    """

    @property
    def engine_name(self) -> str:
        return "sage"

    @property
    def has_precursor_id(self) -> bool:
        """Sage scannr is NOT the timsTOF precursor_id."""
        return False

    def find_result_files(self, base_dir: Path, accession: str) -> list[Path]:
        """Find Sage results.sage.parquet file."""
        sage_path = base_dir / accession / "sage" / "results.sage.parquet"
        if sage_path.exists():
            return [sage_path]
        return []

    def parse(
        self,
        base_dir: Path,
        accession: str,
        raw_file_filter: Optional[str] = None,
    ) -> pd.DataFrame:
        """Parse Sage results into standardized DataFrame.

        IMPORTANT:
        - m/z is calculated from experimental mass and charge.
          This is intentional - Sage reports mass, not m/z directly.
        - RT is converted from minutes to seconds at parse time.

        Args:
            base_dir: Base directory for processed data
            accession: PRIDE accession
            raw_file_filter: Optional filter for specific raw file

        Returns:
            DataFrame with sage_* prefixed columns
        """
        sage_path = base_dir / accession / "sage" / "results.sage.parquet"

        if not sage_path.exists():
            print(f"  Sage results not found: {sage_path}")
            return pd.DataFrame()

        df = pd.read_parquet(sage_path)

        # Filter decoys
        if "is_decoy" in df.columns:
            df = df[~df["is_decoy"]]

        # Filter by raw file if specified
        if raw_file_filter and "filename" in df.columns:
            df = df[df["filename"].str.contains(raw_file_filter, case=False, na=False)]

        if df.empty:
            return pd.DataFrame()

        df = df.copy()

        # INTENTIONAL: Calculate precursor m/z from experimental mass and charge
        # Sage reports expmass (neutral mass), we need m/z for matching
        # Formula: m/z = (M + z * H+) / z where H+ = 1.007276 Da
        df["calculated_mz"] = (df["expmass"] + df["charge"] * PROTON_MASS) / df["charge"]

        # Standardize Sage modified sequence to UNIMOD format
        # Sage "peptide" column contains [+mass] format, e.g., "PEPTIDE[+15.9949]K"
        peptide_col = df.get("peptide", pd.Series(dtype=str))
        df["modified_std"] = peptide_col.apply(standardize_sage_sequence)

        # Strip modifications to get plain peptide sequence
        import re
        def strip_mods(seq):
            if pd.isna(seq) or not seq:
                return ""
            # Remove anything in brackets: [+57.021465] or [UNIMOD:4]
            return re.sub(r'\[[^\]]+\]', '', str(seq))
        df["stripped_peptide"] = peptide_col.apply(strip_mods)

        # Extract raw_file from filename
        df["raw_file"] = df.get("filename", "").apply(
            lambda x: Path(x).stem.replace(".d", "") if pd.notna(x) else ""
        )

        # Convert posterior_error to PEP (it's stored as log)
        pep_values = np.nan
        if "posterior_error" in df.columns:
            pep_values = np.exp(df["posterior_error"])

        # CRITICAL: Convert RT from minutes to seconds
        # Sage reports RT in minutes (confirmed: data range 0-120 for 2h runs)
        rt_minutes = df.get("rt")
        rt_seconds = rt_minutes * 60.0 if rt_minutes is not None else None

        # Build result DataFrame with standardized columns
        result = pd.DataFrame({
            "raw_file": df["raw_file"],
            "sage_psm_id": df.get("psm_id"),  # Links to matched_fragments.sage.parquet
            "sage_scannr": df.get("scannr"),  # Keep for reference, NOT for matching
            "sage_peptide": df.get("stripped_peptide"),
            "sage_modified": df["modified_std"],
            "sage_protein": df.get("proteins"),
            "sage_charge": df.get("charge"),
            "sage_mz": df["calculated_mz"],  # Calculated from mass
            "sage_qvalue": df.get("spectrum_q"),
            "sage_peptide_qvalue": df.get("peptide_q"),
            "sage_protein_qvalue": df.get("protein_q"),
            "sage_pep": pep_values,
            "sage_hyperscore": df.get("hyperscore"),
            "sage_rt": rt_seconds,  # Converted from minutes to seconds
            "sage_mobility": df.get("ion_mobility"),
        })

        return result
