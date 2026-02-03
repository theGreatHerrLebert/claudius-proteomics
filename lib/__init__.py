"""
San José Pipeline Library

Shared utilities for the 6-step pipeline:
- sequence_utils: UNIMOD standardization across search engines
- precursor_matching: Tolerance-based precursor matching
- quality_metrics: Gaussian fits, isotope cosine similarity
- parquet_utils: Parquet read/write helpers
"""

from lib.sequence_utils import (
    standardize_diann_sequence,
    standardize_sage_sequence,
    standardize_fragpipe_sequence,
    standardize_fragpipe_modified_peptide,
    standardize_modified_sequence,
    remove_modifications,
    normalize_sequence_il,
    sequences_match,
    add_carbamidomethyl_to_cysteine,
    mass_to_unimod,
    MASS_TO_UNIMOD,
    ABSOLUTE_MASS_TO_UNIMOD,
)

from lib.precursor_matching import (
    MatchTier,
    MatchConfig,
    MatchResult,
    PrecursorMatcher,
    normalize_sequence_for_matching,
    compute_mz_diff_ppm,
    match_precursors,
)

from lib.quality_metrics import (
    calculate_moments,
    fit_gaussian,
    compute_isotope_cosine_similarity,
    gaussian,
)

__all__ = [
    # sequence_utils
    "standardize_diann_sequence",
    "standardize_sage_sequence",
    "standardize_fragpipe_sequence",
    "standardize_fragpipe_modified_peptide",
    "standardize_modified_sequence",
    "remove_modifications",
    "normalize_sequence_il",
    "sequences_match",
    "add_carbamidomethyl_to_cysteine",
    "mass_to_unimod",
    "MASS_TO_UNIMOD",
    "ABSOLUTE_MASS_TO_UNIMOD",
    # precursor_matching
    "MatchTier",
    "MatchConfig",
    "MatchResult",
    "PrecursorMatcher",
    "normalize_sequence_for_matching",
    "compute_mz_diff_ppm",
    "match_precursors",
    # quality_metrics
    "calculate_moments",
    "fit_gaussian",
    "compute_isotope_cosine_similarity",
    "gaussian",
]
