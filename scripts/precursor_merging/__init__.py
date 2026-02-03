"""
Precursor Merging Package

Provides matching and consensus logic for merging search engine results
with raw timsTOF precursor data.

Strategy (FragPipe-anchored):
1. FragPipe has direct precursor_id from Spectrum column - use for direct raw data join
2. DIA-NN/Sage match to FragPipe via (raw_file, sequence_unimod, charge)
3. Coordinate matching as fallback for unmatched entries (m/z + RT + IM)
"""

from .config import MatchConfig, MatchTier, MatchResult
from .sequence_matcher import SequenceMatcher
from .coordinate_matcher import CoordinateMatcher
from .merger import FragPipeAnchoredMerger
from .consensus import calculate_consensus, count_engines

__all__ = [
    'MatchConfig',
    'MatchTier',
    'MatchResult',
    'SequenceMatcher',
    'CoordinateMatcher',
    'FragPipeAnchoredMerger',
    'calculate_consensus',
    'count_engines',
]
