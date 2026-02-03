"""
DIA-NN Report Parser

Parses DIA-NN report.parquet files into standardized format.

Key features:
- RT conversion: DIA-NN reports RT in MINUTES, we convert to SECONDS at parse time
- No direct precursor_id - requires coordinate or sequence matching
- Standardizes (UniMod:X) format to [UNIMOD:X]
- Includes CCS values from DIA-NN predictions
"""

from pathlib import Path
from typing import Optional
import pandas as pd

import sys
sys.path.insert(0, str(Path(__file__).parent.parent))
from sequence_utils import standardize_diann_sequence

from .base import BaseParser


class DiannParser(BaseParser):
    """Parser for DIA-NN report.parquet results.

    DIA-NN does not provide direct precursor_id mapping. Matching to raw data
    requires either:
    1. Sequence + charge matching (if identified by FragPipe/Sage first)
    2. Coordinate matching (m/z + RT + IM)
    """

    @property
    def engine_name(self) -> str:
        return "diann"

    @property
    def has_precursor_id(self) -> bool:
        """DIA-NN does not have direct precursor_id."""
        return False

    def find_result_files(self, base_dir: Path, accession: str) -> list[Path]:
        """Find DIA-NN report.parquet file."""
        report_path = base_dir / accession / "diann" / "report.parquet"
        if report_path.exists():
            return [report_path]
        return []

    def parse(
        self,
        base_dir: Path,
        accession: str,
        raw_file_filter: Optional[str] = None,
    ) -> pd.DataFrame:
        """Parse DIA-NN report into standardized DataFrame.

        IMPORTANT: RT is converted from minutes to seconds at parse time.
        This is the canonical place for this conversion.

        Args:
            base_dir: Base directory for processed data
            accession: PRIDE accession
            raw_file_filter: Optional filter for specific raw file

        Returns:
            DataFrame with diann_* prefixed columns
        """
        report_path = base_dir / accession / "diann" / "report.parquet"

        if not report_path.exists():
            print(f"  DIA-NN report not found: {report_path}")
            return pd.DataFrame()

        df = pd.read_parquet(report_path)

        # Filter by raw file if specified
        if raw_file_filter and "Run" in df.columns:
            df = df[df["Run"].str.contains(raw_file_filter, case=False, na=False)]

        if df.empty:
            return pd.DataFrame()

        # Standardize DIA-NN modified sequence to UNIMOD format
        modified_col = df.get("Modified.Sequence", pd.Series(dtype=str))
        df["modified_std"] = modified_col.apply(standardize_diann_sequence)

        # Extract raw_file from Run column
        df["raw_file"] = df.get("Run", "").apply(
            lambda x: Path(x).stem.replace(".d", "") if pd.notna(x) else ""
        )

        # CRITICAL: Convert RT from minutes to seconds
        # DIA-NN reports RT in minutes, we standardize to seconds
        rt_minutes = df.get("RT")
        rt_seconds = rt_minutes * 60.0 if rt_minutes is not None else None

        # Build result DataFrame with standardized columns
        result = pd.DataFrame({
            "raw_file": df["raw_file"],
            "diann_peptide": df.get("Stripped.Sequence", df.get("Sequence", "")),
            "diann_modified": df["modified_std"],
            "diann_protein": df.get("Protein.Ids", ""),
            "diann_charge": df.get("Precursor.Charge"),
            "diann_mz": df.get("Precursor.Mz"),
            "diann_qvalue": df.get("Q.Value"),
            "diann_pep": df.get("PEP"),
            "diann_global_qvalue": df.get("Global.Q.Value"),
            "diann_pg_qvalue": df.get("PG.Q.Value"),
            "diann_rt": rt_seconds,  # Now in SECONDS
            "diann_mobility": df.get("IM"),
            "diann_ccs": df.get("CCS"),
        })

        return result
