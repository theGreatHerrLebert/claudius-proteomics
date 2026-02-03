"""
Base Engine Parser Protocol

Defines the interface and data structures for search engine result parsers.
All parsers produce standardized output for downstream merging.
"""

from abc import ABC, abstractmethod
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Protocol, runtime_checkable
import pandas as pd


@dataclass
class StandardizedPSM:
    """Standardized PSM representation across all engines.

    This dataclass defines the canonical fields that all engine parsers
    must produce. Engine-specific fields are prefixed with the engine name.
    """
    # Identity
    raw_file: str
    precursor_id: Optional[int]  # Only FragPipe has direct mapping

    # Sequence
    sequence: str               # Plain sequence
    modified_sequence: str      # Sequence with [UNIMOD:X] format
    protein: str

    # Coordinates
    charge: int
    mz: float
    rt_seconds: float           # Normalized to seconds (DIA-NN reports minutes)
    mobility: Optional[float]

    # Quality scores (engine-specific)
    probability: Optional[float]  # FragPipe
    qvalue: Optional[float]       # All engines
    pep: Optional[float]          # Posterior error probability
    hyperscore: Optional[float]   # FragPipe, Sage


@runtime_checkable
class EngineParser(Protocol):
    """Protocol for search engine result parsers.

    Each parser implements this interface to provide standardized
    loading and parsing of search engine results.
    """

    @property
    def engine_name(self) -> str:
        """Return the engine name (lowercase): 'fragpipe', 'diann', 'sage'."""
        ...

    @property
    def has_precursor_id(self) -> bool:
        """Whether this engine provides direct precursor_id mapping."""
        ...

    def find_result_files(self, base_dir: Path, accession: str) -> list[Path]:
        """Find all result files for this engine under the accession directory."""
        ...

    def parse(
        self,
        base_dir: Path,
        accession: str,
        raw_file_filter: Optional[str] = None,
    ) -> pd.DataFrame:
        """Parse engine results into standardized DataFrame.

        Args:
            base_dir: Base directory for processed data
            accession: PRIDE accession
            raw_file_filter: Optional filter for specific raw file (partial match)

        Returns:
            DataFrame with standardized columns:
            - raw_file: Source file name
            - precursor_id: Direct mapping (FragPipe only, None for others)
            - {engine}_peptide: Plain sequence
            - {engine}_modified: Modified sequence in UNIMOD format
            - {engine}_protein: Protein accession(s)
            - {engine}_charge: Precursor charge
            - {engine}_mz: Precursor m/z
            - {engine}_rt: RT in seconds (normalized)
            - {engine}_mobility: Ion mobility
            - {engine}_probability/qvalue/pep: Quality scores
        """
        ...


class BaseParser(ABC):
    """Abstract base class for engine parsers.

    Provides common functionality and enforces the EngineParser protocol.
    """

    @property
    @abstractmethod
    def engine_name(self) -> str:
        """Return the engine name (lowercase)."""
        pass

    @property
    def has_precursor_id(self) -> bool:
        """Whether this engine provides direct precursor_id mapping.

        Only FragPipe has direct precursor_id from the Spectrum column.
        Override in FragPipeParser to return True.
        """
        return False

    @abstractmethod
    def find_result_files(self, base_dir: Path, accession: str) -> list[Path]:
        """Find all result files for this engine."""
        pass

    @abstractmethod
    def parse(
        self,
        base_dir: Path,
        accession: str,
        raw_file_filter: Optional[str] = None,
    ) -> pd.DataFrame:
        """Parse engine results into standardized DataFrame."""
        pass

    def _prefix_columns(self, df: pd.DataFrame, exclude: list[str] = None) -> pd.DataFrame:
        """Prefix all columns with engine name except excluded ones.

        Args:
            df: DataFrame to modify
            exclude: Columns to keep without prefix (e.g., ['raw_file', 'precursor_id'])

        Returns:
            DataFrame with prefixed columns
        """
        exclude = exclude or ['raw_file', 'precursor_id']
        prefix = f"{self.engine_name}_"

        new_columns = {}
        for col in df.columns:
            if col in exclude:
                new_columns[col] = col
            elif not col.startswith(prefix):
                new_columns[col] = f"{prefix}{col}"
            else:
                new_columns[col] = col

        return df.rename(columns=new_columns)
