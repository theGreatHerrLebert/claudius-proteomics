"""
Engine Parsers Package

Provides standardized parsing of search engine results (FragPipe, DIA-NN, Sage)
into a common format for precursor merging.

Each parser produces a DataFrame with standardized columns:
- raw_file: Source file name (without .d extension)
- precursor_id: Direct mapping if available (FragPipe only)
- sequence: Plain peptide sequence
- modified_sequence: Sequence with UNIMOD modifications
- protein: Protein accession(s)
- charge: Precursor charge state
- mz: Precursor m/z
- rt_seconds: Retention time in seconds (all engines normalized to seconds)
- mobility: Ion mobility (1/K0)
- probability/qvalue/pep: Engine-specific quality scores
"""

from .base import EngineParser, StandardizedPSM
from .fragpipe_parser import FragPipeParser
from .diann_parser import DiannParser
from .sage_parser import SageParser

__all__ = [
    'EngineParser',
    'StandardizedPSM',
    'FragPipeParser',
    'DiannParser',
    'SageParser',
]
