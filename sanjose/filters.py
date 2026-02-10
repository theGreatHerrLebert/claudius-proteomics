"""Filter builder and column group constants for PyArrow queries."""

from typing import Optional, List, Tuple, Any

# Column group constants — matching precursor_store.parquet schema

IDENTITY_COLUMNS = [
    'precursor_id', 'raw_file', 'charge', 'mono_mz', 'mz',
    'isolation_mz', 'isolation_width', 'precursor_intensity',
    'rt_seconds', 'mobility',
]

ENGINE_COLUMNS = [
    'fragpipe_modified', 'fragpipe_peptide', 'fragpipe_protein', 'fragpipe_probability',
    'diann_modified', 'diann_peptide', 'diann_protein', 'diann_qvalue',
    'sage_modified', 'sage_peptide', 'sage_protein', 'sage_qvalue',
]

CONSENSUS_COLUMNS = [
    'n_engines', 'consensus_peptide', 'sequence_normalized',
    'confidence_weight', 'is_high_quality',
]

QUALITY_COLUMNS = [
    'ms1_rt_apex', 'ms1_rt_fwhm', 'ms1_rt_skew', 'ms1_rt_sigma', 'ms1_rt_r2',
    'ms1_im_apex', 'ms1_im_fwhm', 'ms1_im_skew', 'ms1_im_sigma', 'ms1_im_r2',
    'ms1_total_intensity', 'isotope_cosim',
    'ms1_iso_0', 'ms1_iso_1', 'ms1_iso_2', 'ms1_iso_3', 'ms1_iso_4',
]

FRAGMENT_SUMMARY_COLUMNS = [
    'n_fragments_merged', 'n_peaks', 'fragment_total_intensity',
    'fragment_mz_mean', 'fragment_mz_var',
    'fragment_im_mean', 'fragment_im_var', 'fragment_im_apex', 'fragment_im_fwhm',
]

BLOB_COLUMNS = [
    'blob_offset', 'blob_size',
]

# All scalar columns (no list/array columns) — suitable for batch iteration
ALL_SCALAR_COLUMNS = (
    IDENTITY_COLUMNS + FRAGMENT_SUMMARY_COLUMNS + QUALITY_COLUMNS
    + BLOB_COLUMNS + ENGINE_COLUMNS + CONSENSUS_COLUMNS
)

# Lightweight summary for quick browsing
SUMMARY_COLUMNS = [
    'precursor_id', 'raw_file', 'charge', 'mono_mz', 'rt_seconds', 'mobility',
    'n_engines', 'consensus_peptide', 'is_high_quality',
]


def build_filters(
    *,
    sequence: Optional[str] = None,
    charge: Optional[int] = None,
    min_engines: Optional[int] = None,
    max_engines: Optional[int] = None,
    raw_file: Optional[str] = None,
    is_high_quality: Optional[bool] = None,
    min_rt_r2: Optional[float] = None,
    min_im_r2: Optional[float] = None,
    min_isotope_cosim: Optional[float] = None,
    precursor_id: Optional[int] = None,
) -> Optional[List[Tuple[str, str, Any]]]:
    """Build PyArrow filter tuples for predicate pushdown.

    Returns None if no filters specified (read everything).
    """
    filters = []

    if sequence is not None:
        filters.append(('consensus_peptide', '=', sequence))

    if charge is not None:
        filters.append(('charge', '=', charge))

    if min_engines is not None:
        filters.append(('n_engines', '>=', min_engines))

    if max_engines is not None:
        filters.append(('n_engines', '<=', max_engines))

    if raw_file is not None:
        filters.append(('raw_file', '=', raw_file))

    if is_high_quality is not None:
        filters.append(('is_high_quality', '=', is_high_quality))

    if min_rt_r2 is not None:
        filters.append(('ms1_rt_r2', '>=', min_rt_r2))

    if min_im_r2 is not None:
        filters.append(('ms1_im_r2', '>=', min_im_r2))

    if min_isotope_cosim is not None:
        filters.append(('isotope_cosim', '>=', min_isotope_cosim))

    if precursor_id is not None:
        filters.append(('precursor_id', '=', float(precursor_id)))

    return filters if filters else None
