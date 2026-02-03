"""
FragPipe PSM Parser

Parses FragPipe psm.tsv files into standardized format.

Key features:
- Direct precursor_id extraction from Spectrum column (format: rawfile.scannum.scannum.charge)
- Handles both "Assigned Modifications" and "Modified Peptide" columns
- Standardizes modifications to UNIMOD format
- RT already in seconds (no conversion needed)
"""

from pathlib import Path
from typing import Optional
import pandas as pd
import numpy as np

import sys
sys.path.insert(0, str(Path(__file__).parent.parent))
from sequence_utils import (
    standardize_fragpipe_sequence,
    standardize_fragpipe_modified_peptide,
)

from .base import BaseParser


def parse_spectrum_id(spectrum: str) -> tuple[Optional[str], Optional[int], Optional[int]]:
    """Parse FragPipe spectrum ID to extract raw file, precursor ID, and charge.

    FragPipe Spectrum format: rawfile.scannum.scannum.charge
    Example: "sample01.12345.12345.2"

    Args:
        spectrum: FragPipe Spectrum column value

    Returns:
        Tuple of (raw_file, precursor_id, charge)
    """
    if pd.isna(spectrum) or not spectrum:
        return None, None, None

    parts = str(spectrum).rsplit(".", 3)
    if len(parts) >= 4:
        raw_file = parts[0]
        try:
            precursor_id = int(parts[1])
            charge = int(parts[3])
            return raw_file, precursor_id, charge
        except ValueError:
            return parts[0], None, None

    return None, None, None


class FragPipeParser(BaseParser):
    """Parser for FragPipe psm.tsv results.

    FragPipe provides direct precursor_id mapping through its Spectrum column,
    making it the anchor for merging with raw data and other engines.
    """

    @property
    def engine_name(self) -> str:
        return "fragpipe"

    @property
    def has_precursor_id(self) -> bool:
        """FragPipe has direct precursor_id from Spectrum column."""
        return True

    def find_result_files(self, base_dir: Path, accession: str) -> list[Path]:
        """Find all psm.tsv files under the accession directory."""
        processed_dir = base_dir / accession
        return list(processed_dir.rglob("psm.tsv"))

    def parse(
        self,
        base_dir: Path,
        accession: str,
        raw_file_filter: Optional[str] = None,
    ) -> pd.DataFrame:
        """Parse FragPipe PSM files into standardized DataFrame.

        Args:
            base_dir: Base directory for processed data
            accession: PRIDE accession
            raw_file_filter: Optional filter for specific raw file

        Returns:
            DataFrame with fragpipe_* prefixed columns
        """
        psm_files = self.find_result_files(base_dir, accession)

        if not psm_files:
            print(f"  No FragPipe PSM files found in {base_dir / accession}")
            return pd.DataFrame()

        dfs = []
        for psm_file in psm_files:
            df = pd.read_csv(psm_file, sep="\t")

            # Parse spectrum ID to get raw_file, precursor_id, charge
            parsed = df["Spectrum"].apply(parse_spectrum_id)
            df["raw_file"] = parsed.apply(lambda x: x[0])
            df["precursor_id"] = parsed.apply(lambda x: x[1])
            df["extracted_charge"] = parsed.apply(lambda x: x[2])

            # Filter by raw file if specified
            if raw_file_filter:
                df = df[df["raw_file"].str.contains(raw_file_filter, case=False, na=False)]

            if not df.empty:
                dfs.append(df)

        if not dfs:
            return pd.DataFrame()

        all_psms = pd.concat(dfs, ignore_index=True)

        # Keep best PSM per precursor (highest probability)
        all_psms = all_psms.sort_values("Probability", ascending=False)
        all_psms = all_psms.drop_duplicates(subset=["raw_file", "precursor_id"], keep="first")

        # Standardize modifications to UNIMOD format
        all_psms["modified_std"] = self._standardize_modifications(all_psms)

        # Build result DataFrame with standardized columns
        result = pd.DataFrame({
            "raw_file": all_psms["raw_file"],
            "precursor_id": pd.to_numeric(all_psms["precursor_id"], errors="coerce").astype("Int64"),
            "fragpipe_peptide": all_psms["Peptide"],
            "fragpipe_modified": all_psms["modified_std"],
            "fragpipe_protein": all_psms["Protein"],
            "fragpipe_probability": all_psms["Probability"],
            "fragpipe_pep": 1.0 - all_psms["Probability"],
            "fragpipe_hyperscore": all_psms.get("Hyperscore", np.nan),
            "fragpipe_qvalue": all_psms.get("Qvalue", np.nan),
            "fragpipe_mz": all_psms.get("Calibrated Observed M/Z", all_psms.get("Observed M/Z")),
            "fragpipe_rt": all_psms.get("Retention"),  # Already in seconds
            "fragpipe_mobility": all_psms.get("Ion Mobility"),
            "fragpipe_charge": all_psms.get("Charge", all_psms["extracted_charge"]),
        })

        return result

    def _standardize_modifications(self, df: pd.DataFrame) -> pd.Series:
        """Standardize FragPipe modifications to UNIMOD format.

        Priority:
        1. Use Peptide + Assigned Modifications (includes fixed mods like carbamidomethyl)
        2. Fallback to Modified Peptide column (may miss fixed modifications)
        3. Use plain Peptide if no modifications available
        """
        if "Assigned Modifications" in df.columns:
            return df.apply(
                lambda r: standardize_fragpipe_sequence(
                    r["Peptide"],
                    r.get("Assigned Modifications", "") if pd.notna(r.get("Assigned Modifications")) else ""
                ),
                axis=1
            )
        elif "Modified Peptide" in df.columns:
            return df["Modified Peptide"].apply(standardize_fragpipe_modified_peptide)
        else:
            return df["Peptide"]
