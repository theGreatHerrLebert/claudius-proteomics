"""
FastAPI backend for Precursor Browser

Serves precursor data from Parquet store with efficient columnar reads.
Supports both single-dataset mode and collection mode.

Single-dataset mode:
    python main.py --store data/merged/PXD019086/precursor_store.parquet

Collection mode:
    python main.py --collection /path/to/san_jose_collection/
"""

import json
import re
import struct
from pathlib import Path
from typing import Optional, List, Dict, Any, Set
import numpy as np
import pandas as pd

from fastapi import FastAPI, HTTPException, Query
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel
import pyarrow.parquet as pq
import yaml

import gzip
import re

# Try to import zstd for newer blob format
try:
    import zstandard as zstd
    HAS_ZSTD = True
except ImportError:
    HAS_ZSTD = False


# Mass delta to UNIMOD mapping (for standardizing modification display)
MASS_TO_UNIMOD = {
    57: "[UNIMOD:4]",   # Carbamidomethyl (C)
    58: "[UNIMOD:4]",
    16: "[UNIMOD:35]",  # Oxidation (M)
    15: "[UNIMOD:35]",
    42: "[UNIMOD:1]",   # Acetyl (N-term)
    43: "[UNIMOD:1]",
    80: "[UNIMOD:21]",  # Phospho (S/T/Y)
    79: "[UNIMOD:21]",
}


def standardize_modified_sequence(sequence: str) -> str:
    """Convert any modification format to unified UNIMOD format.

    Handles:
    - Sage: [+57.021465] -> [UNIMOD:4]
    - DIA-NN: (UniMod:35) -> [UNIMOD:35]
    - FragPipe: M[147.0354] -> M[UNIMOD:35]
    """
    if not sequence or pd.isna(sequence):
        return sequence

    result = str(sequence)

    # Sage format: [+X.XXX] or [-X.XXX]
    def replace_sage_mass(match):
        try:
            mass = abs(float(match.group(1)))
            rounded = int(round(mass))
            if rounded in MASS_TO_UNIMOD:
                return MASS_TO_UNIMOD[rounded]
        except ValueError:
            pass
        return match.group(0)

    result = re.sub(r'\[([+-]?\d+\.?\d*)\]', replace_sage_mass, result)

    # DIA-NN format: (UniMod:X) -> [UNIMOD:X]
    result = re.sub(r'\(UniMod:(\d+)\)', r'[UNIMOD:\1]', result, flags=re.IGNORECASE)

    # FragPipe absolute mass: [147.XXX] (oxidized M), [160.XXX] (carbamidomethyl C)
    def replace_fragpipe_mass(match):
        try:
            mass = float(match.group(1))
            rounded = int(round(mass))
            # Common absolute masses
            if 146 <= rounded <= 148:  # Oxidized M
                return "[UNIMOD:35]"
            if 159 <= rounded <= 161:  # Carbamidomethyl C
                return "[UNIMOD:4]"
        except ValueError:
            pass
        return match.group(0)

    result = re.sub(r'\[(\d+\.?\d*)\]', replace_fragpipe_mass, result)

    return result

app = FastAPI(title="Precursor Browser API", version="1.0.0")

# CORS for local development
app.add_middleware(
    CORSMiddleware,
    allow_origins=["http://localhost:5173", "http://localhost:3000"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Global store reference (single dataset mode)
store_path: Optional[Path] = None
parquet_file: Optional[pq.ParquetFile] = None
blob_dir: Optional[Path] = None  # Directory containing extracted/{raw_file}.d/blobs.bin

# Blob file sizes for raw data availability checking
blob_file_sizes: Dict[str, int] = {}  # raw_file_name (without .d) -> blobs.bin file size

# Sage fragment data (loaded from engines/sage/)
sage_fragments_by_psm: Dict[int, List[Dict]] = {}  # psm_id -> list of fragment dicts
sage_data_loaded: bool = False

# Global collection reference (collection mode)
collection_path: Optional[Path] = None
collection_manifest: Optional[Dict[str, Any]] = None
studies_config: Optional[Dict[str, Any]] = None
active_dataset: Optional[str] = None  # Currently loaded dataset accession

# Cache derived overlap identifiers per dataset accession
dataset_identity_cache: Dict[str, Dict[str, Any]] = {}


def load_sage_fragments(sage_dir: Path) -> bool:
    """Load Sage matched fragments for fragment annotation.

    Args:
        sage_dir: Path to directory containing Sage output files

    Looks for:
    - matched_fragments.sage.parquet

    Builds lookup dict:
    - sage_fragments_by_psm: psm_id -> list of fragment dicts

    The precursor_store must have sage_psm_id column for exact matching.
    """
    global sage_fragments_by_psm, sage_data_loaded

    fragments_path = sage_dir / "matched_fragments.sage.parquet"

    if not fragments_path.exists():
        print(f"  Sage fragments not found at {fragments_path}")
        return False

    try:
        # Load fragments
        print(f"  Loading Sage fragments from {fragments_path}")
        fragments_df = pq.read_table(str(fragments_path)).to_pandas()
        print(f"    Read {len(fragments_df)} fragment rows, grouping by psm_id...")

        # Vectorized groupby: build list-of-dicts per psm_id without iterrows()
        fragments_df = fragments_df.sort_values('psm_id')
        psm_ids = fragments_df['psm_id'].values
        types = fragments_df['fragment_type'].values
        ordinals = fragments_df['fragment_ordinals'].values.astype(int)
        charges = fragments_df['fragment_charge'].values.astype(int)
        mz_exp = fragments_df['fragment_mz_experimental'].values.astype(float)
        mz_calc = fragments_df['fragment_mz_calculated'].values.astype(float)
        intensities = fragments_df['fragment_intensity'].values.astype(float)

        # Find group boundaries using numpy
        boundaries = np.where(np.diff(psm_ids) != 0)[0] + 1
        starts = np.concatenate([[0], boundaries])
        ends = np.concatenate([boundaries, [len(psm_ids)]])
        unique_psm_ids = psm_ids[starts]

        sage_fragments_by_psm = {}
        for i in range(len(unique_psm_ids)):
            s, e = starts[i], ends[i]
            sage_fragments_by_psm[unique_psm_ids[i]] = [
                {
                    'fragment_type': types[j],
                    'ion_number': ordinals[j],
                    'charge': charges[j],
                    'mz_experimental': mz_exp[j],
                    'mz_calculated': mz_calc[j],
                    'intensity': intensities[j],
                }
                for j in range(s, e)
            ]
        del fragments_df  # Free memory
        print(f"    Loaded fragments for {len(sage_fragments_by_psm)} PSMs")

        sage_data_loaded = True
        return True

    except Exception as e:
        print(f"  Error loading Sage data: {e}")
        return False


class SageMatchedFragment(BaseModel):
    """A matched b/y ion from Sage search results."""
    fragment_type: str      # "b" or "y"
    ion_number: int         # fragment_ordinals (1, 2, 3...)
    charge: int             # fragment_charge
    mz_experimental: float
    mz_calculated: float
    intensity: float


def compute_spectrum_cosine(
    obs_mz: List[float], obs_int: List[float],
    sage_fragments: List[dict],
    ppm_tol: float = 20.0,
) -> Optional[float]:
    """Cosine similarity between observed (blob) and Sage-matched fragment spectra.

    For each Sage fragment, finds the closest observed peak within ppm_tol
    and builds matched intensity vectors for cosine computation.
    """
    if not obs_mz or not sage_fragments:
        return None

    obs_mz_arr = np.array(obs_mz)
    obs_int_arr = np.array(obs_int)

    # Sort observed spectrum by m/z for efficient searching
    sort_idx = np.argsort(obs_mz_arr)
    obs_mz_sorted = obs_mz_arr[sort_idx]
    obs_int_sorted = obs_int_arr[sort_idx]

    matched_obs = []
    matched_sage = []

    for frag in sage_fragments:
        sage_mz = frag['mz_experimental']
        sage_int = frag['intensity']
        tol = sage_mz * ppm_tol / 1e6

        # Binary search for closest peak
        idx = np.searchsorted(obs_mz_sorted, sage_mz)
        best_idx = None
        best_delta = tol

        for candidate in [idx - 1, idx]:
            if 0 <= candidate < len(obs_mz_sorted):
                delta = abs(obs_mz_sorted[candidate] - sage_mz)
                if delta < best_delta:
                    best_delta = delta
                    best_idx = candidate

        if best_idx is not None:
            matched_obs.append(obs_int_sorted[best_idx])
            matched_sage.append(sage_int)

    if len(matched_obs) < 2:
        return None

    a = np.array(matched_obs)
    b = np.array(matched_sage)
    norm_a = np.linalg.norm(a)
    norm_b = np.linalg.norm(b)

    if norm_a == 0 or norm_b == 0:
        return None

    return float(np.dot(a, b) / (norm_a * norm_b))


def one_over_k0_to_ccs(
    one_over_k0: float, mz: float, charge: int,
    mass_gas: float = 28.013, temp: float = 31.85, t_diff: float = 273.15,
) -> float:
    """Convert reduced ion mobility (1/K0) to CCS using Mason-Schamp equation.

    Uses the same constants as imspy/rustims (mscore/src/chemistry/formulas.rs).
    """
    if one_over_k0 <= 0 or charge <= 0 or mz <= 0:
        return 0.0
    summary_constant = 18509.8632163405
    reduced_mobility = 1.0 / one_over_k0
    ion_mass = mz * charge
    reduced_mass = (ion_mass * mass_gas) / (ion_mass + mass_gas)
    return summary_constant * charge / (np.sqrt(reduced_mass * (temp + t_diff)) * reduced_mobility)


class PrecursorSummary(BaseModel):
    precursor_id: int
    raw_file: str
    mz: float
    charge: int
    rt_seconds: float
    mobility: float
    ccs: Optional[float] = None
    n_engines: int
    confidence_weight: Optional[float] = None
    sage_cosine: Optional[float] = None
    # FragPipe (10 columns)
    fragpipe_peptide: Optional[str] = None
    fragpipe_modified: Optional[str] = None
    fragpipe_protein: Optional[str] = None
    fragpipe_probability: Optional[float] = None
    fragpipe_pep: Optional[float] = None
    fragpipe_hyperscore: Optional[float] = None
    fragpipe_qvalue: Optional[float] = None
    fragpipe_rt: Optional[float] = None
    fragpipe_mz: Optional[float] = None
    fragpipe_mobility: Optional[float] = None
    # Sage (12 columns)
    sage_peptide: Optional[str] = None
    sage_modified: Optional[str] = None
    sage_protein: Optional[str] = None
    sage_qvalue: Optional[float] = None
    sage_pep: Optional[float] = None
    sage_hyperscore: Optional[float] = None
    sage_peptide_qvalue: Optional[float] = None
    sage_protein_qvalue: Optional[float] = None
    sage_rt: Optional[float] = None
    sage_mz: Optional[float] = None
    sage_mobility: Optional[float] = None
    sage_match_tier: Optional[str] = None
    # DIA-NN (13 columns)
    diann_peptide: Optional[str] = None
    diann_modified: Optional[str] = None
    diann_protein: Optional[str] = None
    diann_qvalue: Optional[float] = None
    diann_pep: Optional[float] = None
    diann_global_qvalue: Optional[float] = None
    diann_pg_qvalue: Optional[float] = None
    diann_rt: Optional[float] = None
    diann_mz: Optional[float] = None
    diann_mobility: Optional[float] = None
    diann_ccs: Optional[float] = None
    diann_match_tier: Optional[str] = None
    diann_match_score: Optional[float] = None
    # Raw
    collision_energy: Optional[float] = None
    raw_intensity_meta: Optional[float] = None
    frame_id: Optional[int] = None
    isolation_mz: Optional[float] = None
    # Quality metrics
    ms1_rt_sigma: Optional[float] = None
    ms1_rt_r2: Optional[float] = None
    ms1_im_sigma: Optional[float] = None
    ms1_im_r2: Optional[float] = None
    isotope_cosim: Optional[float] = None


class PrecursorDetail(BaseModel):
    precursor_id: int
    mz: float
    charge: int
    rt_seconds: float
    mobility: float
    n_engines: int
    fragpipe_peptide: Optional[str] = None
    fragpipe_modified: Optional[str] = None
    sage_peptide: Optional[str] = None
    sage_modified: Optional[str] = None
    diann_peptide: Optional[str] = None
    diann_modified: Optional[str] = None

    # Fragment spectrum
    fragment_mz: List[float]
    fragment_intensity: List[float]
    fragment_mobility: List[float]
    fragment_scan: List[int]

    # Sage matched b/y ions (if available)
    sage_matched_fragments: Optional[List[SageMatchedFragment]] = None
    sage_cosine: Optional[float] = None  # Cosine similarity: observed vs Sage-matched spectrum

    # MS1 projections
    xic_rt: List[float]
    xic_intensity: List[float]
    mobilogram_im: List[float]
    mobilogram_intensity: List[float]
    isotope_mz: List[float]
    isotope_intensity: List[float]

    # Raw 4D data
    raw_rt: List[float]
    raw_mz: List[float]
    raw_mobility: List[float]
    raw_intensity: List[float]


class StoreInfo(BaseModel):
    path: str
    num_precursors: int
    num_row_groups: int
    columns: List[str]


# Collection mode models
class DatasetQuality(BaseModel):
    rt_r2_median: Optional[float] = None
    im_r2_median: Optional[float] = None
    pct_high_quality: Optional[float] = None


class DatasetSummary(BaseModel):
    accession: str
    version: str
    study_id: str
    path: str
    n_precursors: int
    n_all_three: Optional[int] = None
    n_at_least_two: Optional[int] = None
    quality: Optional[DatasetQuality] = None
    added_at: Optional[str] = None


class StudySummary(BaseModel):
    id: str
    title: str
    organism: Optional[str] = None
    publication: Optional[str] = None
    description: Optional[str] = None
    n_datasets: int
    n_total_precursors: int
    datasets: List[str]


class CollectionInfo(BaseModel):
    version: str
    updated_at: str
    n_studies: int
    n_datasets: int
    n_total_precursors: int


class DatasetOverlapSummary(BaseModel):
    dataset_a: str
    dataset_b: str
    precursors_a: int
    precursors_b: int
    shared_precursors: int
    unique_precursors_a: int
    unique_precursors_b: int
    precursor_jaccard: float
    peptides_a: int
    peptides_b: int
    shared_peptides: int
    unique_peptides_a: int
    unique_peptides_b: int
    peptide_jaccard: float


def _clean_str(value: Any) -> Optional[str]:
    """Return a stripped string value or None for empty/NaN."""
    if value is None:
        return None
    try:
        if pd.isna(value):
            return None
    except Exception:
        pass
    text = str(value).strip()
    if not text or text.lower() == "nan":
        return None
    return text


def _to_float(value: Any) -> Optional[float]:
    """Convert numeric-like value to float, returning None on failure/NaN."""
    if value is None:
        return None
    try:
        if pd.isna(value):
            return None
    except Exception:
        pass
    try:
        return float(value)
    except Exception:
        return None


def normalize_peptide_sequence(sequence: str) -> str:
    """Normalize a modified peptide string to bare amino-acid sequence."""
    if not sequence:
        return ""
    seq = standardize_modified_sequence(sequence)
    # Remove UNIMOD and other bracketed/parenthesized modification annotations.
    seq = re.sub(r"\[UNIMOD:\d+\]", "", seq, flags=re.IGNORECASE)
    seq = re.sub(r"\[[^\]]*\]", "", seq)
    seq = re.sub(r"\([^)]*\)", "", seq)
    return "".join(ch for ch in seq if "A" <= ch <= "Z")


IDENTITY_SEQUENCE_COLUMNS = [
    "fragpipe_modified",
    "sage_modified",
    "diann_modified",
    "fragpipe_peptide",
    "sage_peptide",
    "diann_peptide",
]


def _select_sequence_value(row: Any, sequence_cols: List[str]) -> Optional[str]:
    """Pick first available sequence value from prioritized columns."""
    for col in sequence_cols:
        candidate = _clean_str(row.get(col))
        if candidate:
            return candidate
    return None


def _build_precursor_identity_key(
    row: Any,
    sequence_cols: List[str],
    charge_col: Optional[str],
    mz_col: Optional[str],
    rt_col: Optional[str],
    im_col: Optional[str],
) -> Optional[str]:
    """Build a stable precursor identity key for overlap matching."""
    charge_val = _to_float(row.get(charge_col)) if charge_col else None
    charge = int(charge_val) if charge_val is not None and charge_val > 0 else None

    seq_value = _select_sequence_value(row, sequence_cols)
    if seq_value and charge is not None:
        seq_key = standardize_modified_sequence(seq_value)
        if seq_key:
            return f"{seq_key}|z{charge}"

    # Fallback identity for rows without peptide strings.
    if charge is not None and mz_col:
        mz_val = _to_float(row.get(mz_col))
        if mz_val is not None:
            rt_val = _to_float(row.get(rt_col)) if rt_col else None
            im_val = _to_float(row.get(im_col)) if im_col else None
            rt_key = f"{rt_val:.1f}" if rt_val is not None else "na"
            im_key = f"{im_val:.3f}" if im_val is not None else "na"
            return f"coord:{mz_val:.4f}|z{charge}|rt:{rt_key}|im:{im_key}"

    return None


def _get_collection_dataset_entry(accession: str) -> Dict[str, Any]:
    """Resolve dataset manifest entry in collection mode."""
    if collection_manifest is None:
        raise HTTPException(400, "No collection loaded. Start with --collection flag.")
    dataset = next((d for d in collection_manifest.get("datasets", []) if d.get("accession") == accession), None)
    if dataset is None:
        raise HTTPException(404, f"Dataset not found: {accession}")
    return dataset


def _resolve_dataset_store_path(accession: str) -> Path:
    """Get the precursor_store.parquet path for a collection dataset accession."""
    if collection_path is None:
        raise HTTPException(400, "No collection loaded. Start with --collection flag.")
    dataset = _get_collection_dataset_entry(accession)
    store = collection_path / dataset.get("path", "") / "precursor_store.parquet"
    if not store.exists():
        raise HTTPException(404, f"Precursor store not found: {store}")
    return store


def _extract_dataset_identity_sets(dataset_store: Path) -> Dict[str, Set[str]]:
    """Build overlap-ready precursor and peptide identity sets for one dataset."""
    pf = pq.ParquetFile(str(dataset_store))
    available = set(pf.schema_arrow.names)

    sequence_cols = [c for c in IDENTITY_SEQUENCE_COLUMNS if c in available]
    charge_col = "charge" if "charge" in available else ("raw_charge" if "raw_charge" in available else None)
    mz_col = "mz" if "mz" in available else ("raw_mz" if "raw_mz" in available else None)
    rt_col = "rt_seconds" if "rt_seconds" in available else ("raw_rt_seconds" if "raw_rt_seconds" in available else None)
    im_col = "mobility" if "mobility" in available else ("raw_mobility" if "raw_mobility" in available else None)

    cols_to_read = [c for c in sequence_cols if c in available]
    for col in [charge_col, mz_col, rt_col, im_col]:
        if col and col not in cols_to_read:
            cols_to_read.append(col)

    if not cols_to_read:
        return {"precursor_keys": set(), "peptide_keys": set()}

    table = pq.read_table(str(dataset_store), columns=cols_to_read)
    df = table.to_pandas()

    precursor_keys: Set[str] = set()
    peptide_keys: Set[str] = set()

    for _, row in df.iterrows():
        seq_value = _select_sequence_value(row, sequence_cols)

        if seq_value:
            pep_key = normalize_peptide_sequence(seq_value)
            if pep_key:
                peptide_keys.add(pep_key)

        precursor_key = _build_precursor_identity_key(
            row=row,
            sequence_cols=sequence_cols,
            charge_col=charge_col,
            mz_col=mz_col,
            rt_col=rt_col,
            im_col=im_col,
        )
        if precursor_key:
            precursor_keys.add(precursor_key)

    return {"precursor_keys": precursor_keys, "peptide_keys": peptide_keys}


def _get_dataset_identity(accession: str) -> Dict[str, Set[str]]:
    """Get or compute cached overlap identity sets for one dataset accession."""
    if accession in dataset_identity_cache:
        return dataset_identity_cache[accession]
    dataset_store = _resolve_dataset_store_path(accession)
    identity = _extract_dataset_identity_sets(dataset_store)
    dataset_identity_cache[accession] = identity
    return identity


@app.get("/")
async def root():
    return {
        "status": "ok",
        "mode": "collection" if collection_path else "single",
        "store_loaded": parquet_file is not None,
        "active_dataset": active_dataset,
        "collection_loaded": collection_manifest is not None,
    }


@app.post("/load")
async def load_store(path: str):
    """Load a Parquet store."""
    global store_path, parquet_file

    p = Path(path)
    if not p.exists():
        raise HTTPException(404, f"Store not found: {path}")

    try:
        parquet_file = pq.ParquetFile(str(p))
        store_path = p
        return {"status": "loaded", "path": str(p), "num_precursors": parquet_file.metadata.num_rows}
    except Exception as e:
        raise HTTPException(500, f"Failed to load store: {e}")


@app.get("/info", response_model=StoreInfo)
async def get_info():
    """Get store metadata."""
    if parquet_file is None:
        raise HTTPException(400, "No store loaded")

    return StoreInfo(
        path=str(store_path),
        num_precursors=parquet_file.metadata.num_rows,
        num_row_groups=parquet_file.metadata.num_row_groups,
        columns=parquet_file.schema_arrow.names,
    )


@app.get("/precursors", response_model=List[PrecursorSummary])
async def list_precursors(
    offset: int = Query(0, ge=0),
    limit: int = Query(100, ge=1, le=1000),
    min_engines: int = Query(0, ge=0, le=3),
    max_engines: Optional[int] = Query(None, ge=0, le=3, description="Maximum number of engines (use 0 for unidentified)"),
    charge: Optional[int] = Query(None, ge=1, le=5),
    raw_file: Optional[str] = Query(None, description="Filter by raw file name"),
    has_ms1: bool = Query(False, description="Only show precursors with MS1 signal data"),
    has_raw_data: bool = Query(False, description="Only show precursors with readable raw blob data"),
    sort_by: str = Query("quality", pattern="^(quality|raw_intensity_meta|precursor_intensity|n_engines|mz|rt_seconds|precursor_id|mobility|fragpipe_probability|fragpipe_hyperscore|sage_hyperscore|sage_qvalue|diann_qvalue|diann_match_tier|ms1_rt_r2|ms1_im_r2|isotope_cosim|sage_cosine)$"),
    sort_desc: bool = Query(True),
    overlap_mode: Optional[str] = Query(None, pattern="^(shared|unique_a)$", description="Overlap filter mode relative to overlap_dataset"),
    overlap_dataset: Optional[str] = Query(None, description="Comparison dataset accession for overlap filtering"),
):
    """List precursors with pagination and filtering."""
    if parquet_file is None:
        raise HTTPException(400, "No store loaded")

    # Detect column naming convention by trying to read
    # Try non-prefixed first (newer format)
    try:
        test_table = pq.read_table(str(store_path), columns=['mz'], filters=[('precursor_id', '=', 1)])
        use_raw_prefix = False
    except Exception:
        use_raw_prefix = True

    # Map to actual column names
    mz_col = 'raw_mz' if use_raw_prefix else 'mz'
    charge_col = 'raw_charge' if use_raw_prefix else 'charge'
    rt_col = 'raw_rt_seconds' if use_raw_prefix else 'rt_seconds'
    mobility_col = 'raw_mobility' if use_raw_prefix else 'mobility'
    # Try different intensity column names
    available_cols_check = parquet_file.schema_arrow.names
    if 'precursor_intensity' in available_cols_check:
        intensity_col = 'precursor_intensity'
    elif 'ms1_total_intensity' in available_cols_check:
        intensity_col = 'ms1_total_intensity'
    else:
        intensity_col = 'raw_intensity_meta'

    # Build column list - base columns that must exist
    base_columns = ['precursor_id', 'raw_file', mz_col, charge_col, rt_col, mobility_col, 'n_engines']

    # Optional columns - read what exists
    optional_columns = [
        'confidence_weight',
        # FragPipe (10 columns)
        'fragpipe_peptide', 'fragpipe_modified', 'fragpipe_protein',
        'fragpipe_probability', 'fragpipe_pep', 'fragpipe_hyperscore', 'fragpipe_qvalue',
        'fragpipe_rt', 'fragpipe_mz', 'fragpipe_mobility',
        # Sage (12 columns)
        'sage_peptide', 'sage_modified', 'sage_protein',
        'sage_qvalue', 'sage_pep', 'sage_hyperscore', 'sage_peptide_qvalue', 'sage_protein_qvalue',
        'sage_rt', 'sage_mz', 'sage_mobility', 'sage_match_tier', 'sage_psm_id',
        # DIA-NN (13 columns)
        'diann_peptide', 'diann_modified', 'diann_protein',
        'diann_qvalue', 'diann_pep', 'diann_global_qvalue', 'diann_pg_qvalue',
        'diann_rt', 'diann_mz', 'diann_mobility', 'diann_ccs', 'diann_match_tier', 'diann_match_score',
        # Raw (intensity_col is auto-detected: precursor_intensity, ms1_total_intensity, or raw_intensity_meta)
        'collision_energy', intensity_col, 'frame_id', 'isolation_mz',
        # Quality metrics
        'ms1_rt_sigma', 'ms1_rt_r2', 'ms1_im_sigma', 'ms1_im_r2', 'isotope_cosim',
        'sage_cosine',
        # Blob metadata (needed for has_raw_data filter)
        'blob_offset', 'blob_size',
    ]

    # Get available columns from schema
    available_cols = parquet_file.schema_arrow.names
    columns = base_columns + [c for c in optional_columns if c in available_cols]

    # Build filters
    filters = []
    if min_engines > 0:
        filters.append(('n_engines', '>=', min_engines))
    if max_engines is not None:
        filters.append(('n_engines', '<=', max_engines))
    if charge is not None:
        filters.append((charge_col, '=', charge))
    if raw_file is not None:
        filters.append(('raw_file', '=', raw_file))

    # Read with filters
    table = pq.read_table(str(store_path), columns=columns, filters=filters if filters else None)
    df = table.to_pandas()

    # Filter for precursors with MS1 data
    if has_ms1:
        df = df[df[rt_col].notna()]

    # Filter for precursors with readable raw blob data
    if has_raw_data and blob_file_sizes and 'blob_offset' in df.columns and 'blob_size' in df.columns:
        raw_clean = df['raw_file'].astype(str).str.replace('.d', '', regex=False)
        file_size_series = raw_clean.map(blob_file_sizes).fillna(0)
        blob_end = df['blob_offset'].fillna(0) + df['blob_size'].fillna(0)
        df = df[blob_end <= file_size_series]

    # Optional overlap-based filtering (collection mode only)
    if overlap_mode is not None:
        if overlap_dataset is None:
            raise HTTPException(400, "overlap_dataset is required when overlap_mode is set.")
        if collection_manifest is None or active_dataset is None:
            raise HTTPException(400, "Overlap filtering requires collection mode with an active dataset.")

        compare_identity = _get_dataset_identity(overlap_dataset)
        compare_precursor_keys = compare_identity["precursor_keys"]
        sequence_cols_present = [c for c in IDENTITY_SEQUENCE_COLUMNS if c in df.columns]

        row_keys = [
            _build_precursor_identity_key(
                row=row,
                sequence_cols=sequence_cols_present,
                charge_col=charge_col,
                mz_col=mz_col,
                rt_col=rt_col,
                im_col=mobility_col,
            )
            for _, row in df.iterrows()
        ]
        shared_mask = pd.Series(
            [(k in compare_precursor_keys) if k else False for k in row_keys],
            index=df.index,
        )

        if overlap_mode == "shared":
            df = df[shared_mask]
        else:  # unique_a
            df = df[~shared_mask]

    # Sort - when sorting by n_engines, also sort by intensity within each group
    # Map sort column to actual column name
    sort_col = sort_by
    if sort_by == 'raw_intensity_meta' or sort_by == 'precursor_intensity':
        sort_col = intensity_col
    elif sort_by == 'mz':
        sort_col = mz_col
    elif sort_by == 'rt_seconds':
        sort_col = rt_col
    elif sort_by == 'mobility':
        sort_col = mobility_col

    if sort_by == 'quality':
        # Composite quality score: agreement + spectral quality + intensity
        # Higher is better for all components
        quality_score = pd.Series(0.0, index=df.index)

        # Engine agreement (n_engines/3, weight=0.35)
        if 'n_engines' in df.columns:
            quality_score += (df['n_engines'].fillna(0) / 3.0) * 0.35

        # Confidence weight (already 0-1, weight=0.25)
        if 'confidence_weight' in df.columns:
            quality_score += df['confidence_weight'].fillna(0) * 0.25

        # Isotope cosine similarity (already 0-1, weight=0.20)
        if 'isotope_cosim' in df.columns:
            quality_score += df['isotope_cosim'].fillna(0) * 0.20

        # Intensity (log-normalized, weight=0.20)
        if intensity_col in df.columns:
            int_vals = df[intensity_col].fillna(0).clip(lower=1)
            int_log = np.log10(int_vals)
            int_normalized = (int_log - int_log.min()) / (int_log.max() - int_log.min() + 1e-10)
            quality_score += int_normalized * 0.20

        df['_quality_score'] = quality_score
        df = df.sort_values('_quality_score', ascending=not sort_desc, na_position='last')
        df = df.drop(columns=['_quality_score'])
    elif sort_by == 'n_engines' and intensity_col in df.columns:
        # Special case: sort by n_engines with intensity as secondary
        df = df.sort_values(
            ['n_engines', intensity_col],
            ascending=[not sort_desc, False],
            na_position='last'
        )
    elif sort_col in df.columns:
        df = df.sort_values(sort_col, ascending=not sort_desc, na_position='last')

    # Paginate
    df = df.iloc[offset:offset + limit]

    # Helper to safely get float
    def safe_float(val):
        if val is None or (isinstance(val, float) and np.isnan(val)):
            return None
        try:
            return float(val)
        except (ValueError, TypeError):
            return None

    def safe_str(val):
        if val is None or (isinstance(val, float) and np.isnan(val)):
            return None
        return str(val) if val else None

    def safe_int(val):
        if val is None or (isinstance(val, float) and np.isnan(val)):
            return None
        try:
            return int(val)
        except (ValueError, TypeError):
            return None

    # Convert to response
    results = []
    for idx, row in df.iterrows():
        # Handle NaN precursor_id - use DataFrame index as fallback
        prec_id = row['precursor_id']
        if pd.isna(prec_id):
            prec_id = idx if isinstance(idx, int) else hash(str(idx)) % 1000000000

        # Get mz with fallback to engine-specific values
        mz_val = row.get(mz_col)
        if pd.isna(mz_val):
            for fallback in ['fragpipe_mz', 'sage_mz', 'diann_mz']:
                if fallback in row and pd.notna(row.get(fallback)):
                    mz_val = row[fallback]
                    break

        # Get charge with fallback
        charge_val = row.get(charge_col)
        if pd.isna(charge_val):
            charge_val = 0

        # Get n_engines safely
        n_engines_val = row.get('n_engines', 0)
        if pd.isna(n_engines_val):
            n_engines_val = 0

        # Use precomputed sage_cosine from parquet if available,
        # otherwise fall back to on-the-fly blob-read computation
        row_sage_cosine = safe_float(row.get('sage_cosine'))
        if row_sage_cosine is None:
            sage_psm_val = row.get('sage_psm_id')
            if (sage_data_loaded and pd.notna(sage_psm_val)
                    and pd.notna(row.get('blob_offset')) and pd.notna(row.get('blob_size'))):
                psm_id = int(sage_psm_val)
                if psm_id in sage_fragments_by_psm:
                    blob = read_blob(
                        str(row['raw_file']) if pd.notna(row.get('raw_file')) else '',
                        int(row['blob_offset']),
                        int(row['blob_size']),
                    )
                    if blob:
                        row_sage_cosine = compute_spectrum_cosine(
                            blob.get('fragment_mz', []),
                            blob.get('fragment_intensity', []),
                            sage_fragments_by_psm[psm_id],
                        )

        # Compute CCS from 1/K0, m/z, and charge (Mason-Schamp equation)
        row_mz = float(mz_val) if pd.notna(mz_val) else 0.0
        row_mobility = float(row[mobility_col]) if pd.notna(row.get(mobility_col)) else 0.0
        row_charge = int(charge_val)
        row_ccs = one_over_k0_to_ccs(row_mobility, row_mz, row_charge) if row_mobility > 0 and row_charge > 0 and row_mz > 0 else None

        results.append(PrecursorSummary(
            precursor_id=int(prec_id),
            raw_file=str(row['raw_file']) if pd.notna(row.get('raw_file')) else '',
            mz=row_mz,
            charge=row_charge,
            rt_seconds=float(row[rt_col]) if pd.notna(row.get(rt_col)) else 0.0,
            mobility=row_mobility,
            ccs=row_ccs,
            n_engines=int(n_engines_val),
            confidence_weight=safe_float(row.get('confidence_weight')),
            sage_cosine=row_sage_cosine,
            # FragPipe (10 columns)
            fragpipe_peptide=safe_str(row.get('fragpipe_peptide')),
            fragpipe_modified=standardize_modified_sequence(safe_str(row.get('fragpipe_modified'))),
            fragpipe_protein=safe_str(row.get('fragpipe_protein')),
            fragpipe_probability=safe_float(row.get('fragpipe_probability')),
            fragpipe_pep=safe_float(row.get('fragpipe_pep')),
            fragpipe_hyperscore=safe_float(row.get('fragpipe_hyperscore')),
            fragpipe_qvalue=safe_float(row.get('fragpipe_qvalue')),
            fragpipe_rt=safe_float(row.get('fragpipe_rt')),
            fragpipe_mz=safe_float(row.get('fragpipe_mz')),
            fragpipe_mobility=safe_float(row.get('fragpipe_mobility')),
            # Sage (12 columns)
            sage_peptide=safe_str(row.get('sage_peptide')),
            sage_modified=standardize_modified_sequence(safe_str(row.get('sage_modified'))),
            sage_protein=safe_str(row.get('sage_protein')),
            sage_qvalue=safe_float(row.get('sage_qvalue')),
            sage_pep=safe_float(row.get('sage_pep')),
            sage_hyperscore=safe_float(row.get('sage_hyperscore')),
            sage_peptide_qvalue=safe_float(row.get('sage_peptide_qvalue')),
            sage_protein_qvalue=safe_float(row.get('sage_protein_qvalue')),
            sage_rt=safe_float(row.get('sage_rt')),
            sage_mz=safe_float(row.get('sage_mz')),
            sage_mobility=safe_float(row.get('sage_mobility')),
            sage_match_tier=safe_str(row.get('sage_match_tier')),
            # DIA-NN (13 columns)
            diann_peptide=safe_str(row.get('diann_peptide')),
            diann_modified=standardize_modified_sequence(safe_str(row.get('diann_modified'))),
            diann_protein=safe_str(row.get('diann_protein')),
            diann_qvalue=safe_float(row.get('diann_qvalue')),
            diann_pep=safe_float(row.get('diann_pep')),
            diann_global_qvalue=safe_float(row.get('diann_global_qvalue')),
            diann_pg_qvalue=safe_float(row.get('diann_pg_qvalue')),
            diann_rt=safe_float(row.get('diann_rt')),
            diann_mz=safe_float(row.get('diann_mz')),
            diann_mobility=safe_float(row.get('diann_mobility')),
            diann_ccs=safe_float(row.get('diann_ccs')),
            diann_match_tier=safe_str(row.get('diann_match_tier')),
            diann_match_score=safe_float(row.get('diann_match_score')),
            # Raw
            collision_energy=safe_float(row.get('collision_energy')),
            raw_intensity_meta=safe_float(row.get(intensity_col)),
            frame_id=safe_int(row.get('frame_id')),
            isolation_mz=safe_float(row.get('isolation_mz')),
            # Quality metrics
            ms1_rt_sigma=safe_float(row.get('ms1_rt_sigma')),
            ms1_rt_r2=safe_float(row.get('ms1_rt_r2')),
            ms1_im_sigma=safe_float(row.get('ms1_im_sigma')),
            ms1_im_r2=safe_float(row.get('ms1_im_r2')),
            isotope_cosim=safe_float(row.get('isotope_cosim')),
        ))

    return results


def read_blob(raw_file: str, offset: int, size: int) -> Optional[dict]:
    """Read and decompress a precursor blob from the blobs.bin file.

    Returns dict with xic, mobilogram, isotope, and fragment data.
    """
    if blob_dir is None:
        return None

    # Find blob file for this raw file
    raw_file_clean = raw_file.replace('.d', '')
    blob_path = blob_dir / f"{raw_file_clean}.d" / "blobs.bin"

    if not blob_path.exists():
        # Try with .d suffix
        blob_path = blob_dir / f"{raw_file}.d" / "blobs.bin" if not raw_file.endswith('.d') else blob_dir / raw_file / "blobs.bin"

    if not blob_path.exists():
        return None

    try:
        with open(blob_path, 'rb') as f:
            file_size = f.seek(0, 2)  # Get file size
            if offset + size > file_size:
                return None  # Blob offset beyond file (truncated extraction)
            f.seek(offset)
            compressed = f.read(size)

        if len(compressed) == 0:
            return None

        # Detect compression format and decompress
        # zstd magic: 0x28 0xB5 0x2F 0xFD
        # gzip magic: 0x1F 0x8B
        if compressed[:4] == b'\x28\xb5\x2f\xfd':
            if not HAS_ZSTD:
                print("zstd blob found but zstandard not installed")
                return None
            dctx = zstd.ZstdDecompressor()
            combined = dctx.decompress(compressed)
        elif compressed[:2] == b'\x1f\x8b':
            combined = gzip.decompress(compressed)
        else:
            print(f"Unknown compression format: {compressed[:4].hex()}")
            return None

        # Parse: [4 bytes metadata_len][metadata JSON][npz arrays]
        metadata_len = int.from_bytes(combined[:4], "little")
        metadata_bytes = combined[4:4+metadata_len]
        metadata = json.loads(metadata_bytes.decode("utf-8"))

        # Extract arrays
        import io
        npz_bytes = combined[4+metadata_len:]
        npz_buffer = io.BytesIO(npz_bytes)
        arrays = np.load(npz_buffer)

        # Build result from actual blob structure
        result = {
            # Fragment spectrum (frag_* arrays)
            "fragment_mz": arrays["frag_mz"].tolist() if "frag_mz" in arrays else [],
            "fragment_intensity": arrays["frag_intensity"].tolist() if "frag_intensity" in arrays else [],
            "fragment_mobility": arrays["frag_mobility"].tolist() if "frag_mobility" in arrays else [],
            "fragment_scan": arrays["frag_scan"].tolist() if "frag_scan" in arrays else [],
            # XIC (ms1_rt_* arrays)
            "xic_rt": arrays["ms1_rt_coords"].tolist() if "ms1_rt_coords" in arrays else [],
            "xic_intensity": arrays["ms1_rt_intensities"].tolist() if "ms1_rt_intensities" in arrays else [],
            # Mobilogram (ms1_im_* arrays)
            "mobilogram_im": arrays["ms1_im_coords"].tolist() if "ms1_im_coords" in arrays else [],
            "mobilogram_intensity": arrays["ms1_im_intensities"].tolist() if "ms1_im_intensities" in arrays else [],
            # Isotope envelope (from metadata)
            "isotope_mz": [],  # Not stored separately, would need to compute from mono_mz
            "isotope_intensity": metadata.get("ms1_isotope_intensities", []),
            # Raw 4D MS1 point cloud (for RT vs IM heatmap)
            "raw_rt": arrays["raw_rt"].tolist() if "raw_rt" in arrays else [],
            "raw_mz": arrays["raw_mz"].tolist() if "raw_mz" in arrays else [],
            "raw_mobility": arrays["raw_mobility"].tolist() if "raw_mobility" in arrays else [],
            "raw_intensity": arrays["raw_intensity"].tolist() if "raw_intensity" in arrays else [],
        }

        return result

    except Exception as e:
        print(f"Error reading blob: {e}")
        return None


@app.get("/precursor/{precursor_id}", response_model=PrecursorDetail)
async def get_precursor(
    precursor_id: int,
    raw_file: Optional[str] = Query(None, description="Specific raw file to load data from"),
):
    """Get full detail for a single precursor."""
    if parquet_file is None:
        raise HTTPException(400, "No store loaded")

    # Read with filter
    filters = [('precursor_id', '=', precursor_id)]
    if raw_file:
        filters.append(('raw_file', '=', raw_file))
    table = pq.read_table(str(store_path), filters=filters)

    if len(table) == 0:
        raise HTTPException(404, f"Precursor not found: {precursor_id}")

    df = table.to_pandas()

    def to_list(val):
        if val is None:
            return []
        if isinstance(val, (list, np.ndarray)):
            return [float(x) for x in val]
        return []

    def to_int_list(val):
        if val is None:
            return []
        if isinstance(val, (list, np.ndarray)):
            return [int(x) for x in val]
        return []

    def safe_str(val):
        if val is None or (isinstance(val, float) and np.isnan(val)):
            return None
        return str(val) if val else None

    # When multiple rows exist (same precursor across raw files),
    # prefer a row with readable blob data
    blob_data = None
    row = df.iloc[0]
    if blob_file_sizes and 'blob_offset' in df.columns and 'blob_size' in df.columns:
        for _, candidate in df.iterrows():
            rf = str(candidate.get('raw_file', '')).replace('.d', '')
            fsize = blob_file_sizes.get(rf, 0)
            off = candidate.get('blob_offset')
            sz = candidate.get('blob_size')
            if pd.notna(off) and pd.notna(sz) and int(off) + int(sz) <= fsize:
                row = candidate
                break

    if pd.notna(row.get('blob_offset')) and pd.notna(row.get('blob_size')) and pd.notna(row.get('raw_file')):
        blob_data = read_blob(
            str(row['raw_file']),
            int(row['blob_offset']),
            int(row['blob_size'])
        )

    # Get precursor properties for isotope m/z calculation
    precursor_mz = float(row['mz']) if pd.notna(row.get('mz')) else 0.0
    precursor_charge = int(row['charge']) if pd.notna(row.get('charge')) and int(row['charge']) > 0 else 1

    # Use blob data if available, otherwise try parquet columns, otherwise empty
    if blob_data:
        fragment_mz = blob_data.get("fragment_mz", [])
        fragment_intensity = blob_data.get("fragment_intensity", [])
        fragment_mobility = blob_data.get("fragment_mobility", [])
        fragment_scan = blob_data.get("fragment_scan", [])
        xic_rt = blob_data.get("xic_rt", [])
        xic_intensity = blob_data.get("xic_intensity", [])
        mobilogram_im = blob_data.get("mobilogram_im", [])
        mobilogram_intensity = blob_data.get("mobilogram_intensity", [])
        isotope_intensity = blob_data.get("isotope_intensity", [])
        # Compute isotope m/z from precursor m/z and charge
        # Delta mass between C13 and C12 is ~1.003355 Da
        delta_mass = 1.003355
        isotope_mz = [precursor_mz + (i * delta_mass / precursor_charge) for i in range(len(isotope_intensity))]
        # Raw 4D MS1 point cloud (for RT vs IM heatmap)
        raw_rt = blob_data.get("raw_rt", [])
        raw_mz = blob_data.get("raw_mz", [])
        raw_mobility = blob_data.get("raw_mobility", [])
        raw_intensity = blob_data.get("raw_intensity", [])
    else:
        fragment_mz = to_list(row.get('fragment_mz'))
        fragment_intensity = to_list(row.get('fragment_intensity'))
        fragment_mobility = to_list(row.get('fragment_mobility'))
        fragment_scan = to_int_list(row.get('fragment_scan')) if 'fragment_scan' in row else []
        xic_rt = to_list(row.get('xic_rt'))
        xic_intensity = to_list(row.get('xic_intensity'))
        mobilogram_im = to_list(row.get('mobilogram_im'))
        mobilogram_intensity = to_list(row.get('mobilogram_intensity'))
        isotope_mz = to_list(row.get('isotope_mz'))
        isotope_intensity = to_list(row.get('isotope_intensity'))
        # Raw 4D from parquet if available
        raw_rt = to_list(row.get('raw_rt'))
        raw_mz = to_list(row.get('raw_mz'))
        raw_mobility = to_list(row.get('raw_mobility'))
        raw_intensity = to_list(row.get('raw_intensity'))

    # Get Sage modified sequence (standardize format)
    sage_modified = standardize_modified_sequence(safe_str(row.get('sage_modified')))

    # Get Sage matched fragments if available (exact psm_id match only)
    sage_matched_fragments = None
    sage_cosine = None
    sage_psm_id = row.get('sage_psm_id')
    if sage_data_loaded and pd.notna(sage_psm_id):
        psm_id = int(sage_psm_id)
        if psm_id in sage_fragments_by_psm:
            frag_dicts = sage_fragments_by_psm[psm_id]
            sage_matched_fragments = [
                SageMatchedFragment(**frag) for frag in frag_dicts
            ]
            # Compute cosine similarity between blob spectrum and Sage-matched peaks
            sage_cosine = compute_spectrum_cosine(
                fragment_mz, fragment_intensity, frag_dicts
            )

    return PrecursorDetail(
        precursor_id=int(row['precursor_id']),
        mz=float(row['mz']) if pd.notna(row.get('mz')) else 0.0,
        charge=int(row['charge']) if pd.notna(row.get('charge')) else 0,
        rt_seconds=float(row['rt_seconds']) if pd.notna(row.get('rt_seconds')) else 0.0,
        mobility=float(row['mobility']) if pd.notna(row.get('mobility')) else 0.0,
        n_engines=int(row['n_engines']) if pd.notna(row.get('n_engines')) else 0,
        fragpipe_peptide=safe_str(row.get('fragpipe_peptide')),
        fragpipe_modified=standardize_modified_sequence(safe_str(row.get('fragpipe_modified'))),
        sage_peptide=safe_str(row.get('sage_peptide')),
        sage_modified=sage_modified,
        diann_peptide=safe_str(row.get('diann_peptide')),
        diann_modified=standardize_modified_sequence(safe_str(row.get('diann_modified'))),
        fragment_mz=fragment_mz,
        fragment_intensity=fragment_intensity,
        fragment_mobility=fragment_mobility,
        fragment_scan=fragment_scan,
        sage_matched_fragments=sage_matched_fragments,
        sage_cosine=sage_cosine,
        xic_rt=xic_rt,
        xic_intensity=xic_intensity,
        mobilogram_im=mobilogram_im,
        mobilogram_intensity=mobilogram_intensity,
        isotope_mz=isotope_mz,
        isotope_intensity=isotope_intensity,
        raw_rt=raw_rt,
        raw_mz=raw_mz,
        raw_mobility=raw_mobility,
        raw_intensity=raw_intensity,
    )


@app.get("/stats")
async def get_stats():
    """Get summary statistics."""
    if parquet_file is None:
        raise HTTPException(400, "No store loaded")

    # Read all potential columns and use what exists
    try:
        table = pq.read_table(str(store_path), columns=['n_engines', 'charge', 'mz', 'raw_file'])
        charge_col, mz_col = 'charge', 'mz'
    except Exception:
        table = pq.read_table(str(store_path), columns=['n_engines', 'raw_charge', 'raw_mz', 'raw_file'])
        charge_col, mz_col = 'raw_charge', 'raw_mz'

    df = table.to_pandas()

    return {
        "total_precursors": len(df),
        "by_engines": df['n_engines'].value_counts().to_dict(),
        "by_charge": df[charge_col].value_counts().to_dict(),
        "mz_range": [float(df[mz_col].min()), float(df[mz_col].max())],
        "raw_files": df['raw_file'].value_counts().to_dict(),
    }


@app.get("/raw_files")
async def list_raw_files():
    """Get list of available raw files with counts."""
    if parquet_file is None:
        raise HTTPException(400, "No store loaded")

    table = pq.read_table(str(store_path), columns=['raw_file'])
    df = table.to_pandas()

    files = []
    for name, count in df['raw_file'].value_counts().items():
        files.append({
            "name": str(name),
            "count": int(count),
        })

    return files


# =============================================================================
# Collection Mode Endpoints
# =============================================================================

@app.get("/collection", response_model=CollectionInfo)
async def get_collection_info():
    """Get collection metadata (collection mode only)."""
    if collection_manifest is None:
        raise HTTPException(400, "No collection loaded. Start with --collection flag.")

    datasets = collection_manifest.get("datasets", [])
    studies = collection_manifest.get("studies", [])

    return CollectionInfo(
        version=collection_manifest.get("version", "1.0"),
        updated_at=collection_manifest.get("updated_at", ""),
        n_studies=len(studies),
        n_datasets=len(datasets),
        n_total_precursors=sum(d.get("n_precursors", 0) for d in datasets),
    )


@app.get("/studies", response_model=List[StudySummary])
async def list_studies():
    """List all studies in the collection (collection mode only)."""
    if collection_manifest is None:
        raise HTTPException(400, "No collection loaded. Start with --collection flag.")

    studies = []
    studies_from_manifest = collection_manifest.get("studies", [])

    # Merge with studies.yaml for full details
    studies_yaml = studies_config.get("studies", []) if studies_config else []
    yaml_by_id = {s["id"]: s for s in studies_yaml}

    for study in studies_from_manifest:
        yaml_info = yaml_by_id.get(study["id"], {})
        studies.append(StudySummary(
            id=study["id"],
            title=study.get("title", study["id"]),
            organism=yaml_info.get("organism") or study.get("organism"),
            publication=yaml_info.get("publication"),
            description=yaml_info.get("description"),
            n_datasets=study.get("n_datasets", 0),
            n_total_precursors=study.get("n_total_precursors", 0),
            datasets=study.get("datasets", []),
        ))

    return studies


@app.get("/studies/{study_id}/datasets", response_model=List[DatasetSummary])
async def list_study_datasets(study_id: str):
    """List all datasets in a study (collection mode only)."""
    if collection_manifest is None:
        raise HTTPException(400, "No collection loaded. Start with --collection flag.")

    datasets = collection_manifest.get("datasets", [])
    study_datasets = [d for d in datasets if d.get("study_id") == study_id]

    if not study_datasets:
        # Check if study exists
        studies = collection_manifest.get("studies", [])
        if not any(s["id"] == study_id for s in studies):
            raise HTTPException(404, f"Study not found: {study_id}")

    result = []
    for d in study_datasets:
        quality = None
        if d.get("quality"):
            quality = DatasetQuality(**d["quality"])

        result.append(DatasetSummary(
            accession=d["accession"],
            version=d.get("version", "1.0"),
            study_id=d["study_id"],
            path=d.get("path", ""),
            n_precursors=d.get("n_precursors", 0),
            n_all_three=d.get("n_all_three"),
            n_at_least_two=d.get("n_at_least_two"),
            quality=quality,
            added_at=d.get("added_at"),
        ))

    return result


@app.get("/datasets/{accession}/info", response_model=DatasetSummary)
async def get_dataset_info(accession: str):
    """Get dataset info without loading it (collection mode only)."""
    if collection_manifest is None:
        raise HTTPException(400, "No collection loaded. Start with --collection flag.")

    datasets = collection_manifest.get("datasets", [])
    dataset = next((d for d in datasets if d["accession"] == accession), None)

    if dataset is None:
        raise HTTPException(404, f"Dataset not found: {accession}")

    quality = None
    if dataset.get("quality"):
        quality = DatasetQuality(**dataset["quality"])

    return DatasetSummary(
        accession=dataset["accession"],
        version=dataset.get("version", "1.0"),
        study_id=dataset["study_id"],
        path=dataset.get("path", ""),
        n_precursors=dataset.get("n_precursors", 0),
        n_all_three=dataset.get("n_all_three"),
        n_at_least_two=dataset.get("n_at_least_two"),
        quality=quality,
        added_at=dataset.get("added_at"),
    )


@app.get("/datasets/overlap", response_model=DatasetOverlapSummary)
async def get_dataset_overlap(
    dataset_a: str = Query(..., description="Primary dataset accession"),
    dataset_b: str = Query(..., description="Comparison dataset accession"),
):
    """Compute cross-dataset overlap for precursor and peptide identities."""
    if collection_manifest is None or collection_path is None:
        raise HTTPException(400, "No collection loaded. Start with --collection flag.")
    if dataset_a == dataset_b:
        raise HTTPException(400, "dataset_a and dataset_b must be different.")

    identity_a = _get_dataset_identity(dataset_a)
    identity_b = _get_dataset_identity(dataset_b)

    precursor_a = identity_a["precursor_keys"]
    precursor_b = identity_b["precursor_keys"]
    peptide_a = identity_a["peptide_keys"]
    peptide_b = identity_b["peptide_keys"]

    shared_precursors = len(precursor_a & precursor_b)
    unique_precursors_a = len(precursor_a - precursor_b)
    unique_precursors_b = len(precursor_b - precursor_a)
    precursor_union = len(precursor_a | precursor_b)
    precursor_jaccard = (shared_precursors / precursor_union) if precursor_union > 0 else 0.0

    shared_peptides = len(peptide_a & peptide_b)
    unique_peptides_a = len(peptide_a - peptide_b)
    unique_peptides_b = len(peptide_b - peptide_a)
    peptide_union = len(peptide_a | peptide_b)
    peptide_jaccard = (shared_peptides / peptide_union) if peptide_union > 0 else 0.0

    return DatasetOverlapSummary(
        dataset_a=dataset_a,
        dataset_b=dataset_b,
        precursors_a=len(precursor_a),
        precursors_b=len(precursor_b),
        shared_precursors=shared_precursors,
        unique_precursors_a=unique_precursors_a,
        unique_precursors_b=unique_precursors_b,
        precursor_jaccard=precursor_jaccard,
        peptides_a=len(peptide_a),
        peptides_b=len(peptide_b),
        shared_peptides=shared_peptides,
        unique_peptides_a=unique_peptides_a,
        unique_peptides_b=unique_peptides_b,
        peptide_jaccard=peptide_jaccard,
    )


@app.post("/datasets/{accession}/load")
async def load_dataset(accession: str):
    """Load a dataset from the collection (collection mode only)."""
    global store_path, parquet_file, active_dataset, blob_dir, blob_file_sizes
    global sage_fragments_by_psm, sage_data_loaded

    if collection_path is None:
        raise HTTPException(400, "No collection loaded. Start with --collection flag.")

    if collection_manifest is None:
        raise HTTPException(400, "Collection manifest not loaded")

    # Find dataset in manifest
    datasets = collection_manifest.get("datasets", [])
    dataset = next((d for d in datasets if d["accession"] == accession), None)

    if dataset is None:
        raise HTTPException(404, f"Dataset not found: {accession}")

    # Build path to precursor store
    dataset_dir = collection_path / dataset["path"]
    dataset_path = dataset_dir / "precursor_store.parquet"

    if not dataset_path.exists():
        raise HTTPException(404, f"Precursor store not found: {dataset_path}")

    try:
        parquet_file = pq.ParquetFile(str(dataset_path))
        store_path = dataset_path
        active_dataset = accession

        # Set up blob directory for raw data reading
        import os
        extracted_dir = dataset_dir / "extracted"
        blob_file_sizes = {}
        if extracted_dir.exists():
            blob_dir = extracted_dir
            for d in sorted(extracted_dir.iterdir()):
                if d.is_dir() and (d / "blobs.bin").exists():
                    raw_name = d.name.replace('.d', '')
                    blob_file_sizes[raw_name] = os.path.getsize(d / "blobs.bin")
            print(f"  Blob directory: {blob_dir} ({len(blob_file_sizes)} files)")
        else:
            blob_dir = None
            print(f"  No blob directory found at {extracted_dir}")

        # Load Sage fragment data
        sage_fragments_by_psm = {}
        sage_data_loaded = False
        sage_dir = dataset_dir / "engines" / "sage"
        if sage_dir.exists():
            load_sage_fragments(sage_dir)

        return {
            "status": "loaded",
            "accession": accession,
            "path": str(dataset_path),
            "num_precursors": parquet_file.metadata.num_rows,
            "blob_files": len(blob_file_sizes),
            "sage_fragments": sage_data_loaded,
        }
    except Exception as e:
        raise HTTPException(500, f"Failed to load dataset: {e}")


@app.get("/datasets/active")
async def get_active_dataset():
    """Get the currently active dataset (collection mode only)."""
    return {
        "accession": active_dataset,
        "path": str(store_path) if store_path else None,
        "loaded": parquet_file is not None,
    }


def load_collection(path: Path) -> None:
    """Load a collection from disk."""
    global collection_path, collection_manifest, studies_config, dataset_identity_cache

    collection_path = path
    dataset_identity_cache = {}

    # Load collection manifest
    manifest_path = path / "collection_manifest.json"
    if manifest_path.exists():
        with open(manifest_path) as f:
            collection_manifest = json.load(f)
        print(f"Loaded collection manifest: {len(collection_manifest.get('datasets', []))} datasets")
    else:
        collection_manifest = {"version": "1.0", "studies": [], "datasets": []}
        print("Warning: No collection_manifest.json found")

    # Load studies config
    studies_path = path / "studies.yaml"
    if studies_path.exists():
        with open(studies_path) as f:
            studies_config = yaml.safe_load(f)
        print(f"Loaded studies config: {len(studies_config.get('studies', []))} studies")
    else:
        studies_config = {"studies": []}
        print("Warning: No studies.yaml found")


if __name__ == "__main__":
    import uvicorn
    import argparse

    parser = argparse.ArgumentParser(
        description="Precursor Browser API Server",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Modes:
    Single dataset:
        python main.py --store data/merged/PXD019086/precursor_store.parquet

    Collection mode:
        python main.py --collection /path/to/san_jose_collection/
        """,
    )
    parser.add_argument("--store", help="Path to Parquet store (single dataset mode)")
    parser.add_argument("--blob-dir", help="Path to extracted/ directory containing blobs")
    parser.add_argument("--sage-dir", help="Path to engines/sage/ directory containing matched_fragments")
    parser.add_argument("--collection", help="Path to collection root (collection mode)")
    parser.add_argument("--port", type=int, default=8000)
    parser.add_argument("--host", default="127.0.0.1")
    args = parser.parse_args()

    if args.collection:
        # Collection mode
        p = Path(args.collection)
        if not p.exists():
            print(f"Error: Collection path not found: {p}")
            exit(1)
        load_collection(p)
        print(f"Running in collection mode: {p}")

    elif args.store:
        # Single dataset mode
        p = Path(args.store)
        if p.exists():
            parquet_file = pq.ParquetFile(str(p))
            store_path = p
            print(f"Loaded store: {p} ({parquet_file.metadata.num_rows} precursors)")

            # Set blob directory - try to auto-detect from store path
            if args.blob_dir:
                blob_dir = Path(args.blob_dir)
            else:
                # Auto-detect: if store is in PXD019086/, look for PXD019086/extracted/
                blob_dir = p.parent / "extracted"

            if blob_dir.exists():
                print(f"Blob directory: {blob_dir}")
                print("  Blob reading enabled (auto-detect gzip/zstd)")
                # Index blob file sizes for raw data availability filtering
                import os
                for d in sorted(blob_dir.iterdir()):
                    if d.is_dir() and (d / "blobs.bin").exists():
                        raw_name = d.name.replace('.d', '')
                        blob_file_sizes[raw_name] = os.path.getsize(d / "blobs.bin")
                if blob_file_sizes:
                    print(f"  Indexed {len(blob_file_sizes)} blob files for raw data filtering")
            else:
                blob_dir = None
                print(f"Warning: Blob directory not found: {p.parent / 'extracted'}")

            # Try to load Sage fragment data
            if args.sage_dir:
                sage_dir = Path(args.sage_dir)
            else:
                # Auto-detect: try parent/engines/sage/ or sibling PXD*/engines/sage/
                sage_dir = p.parent / "engines" / "sage"
            sage_loaded = load_sage_fragments(sage_dir)
            if sage_loaded:
                print("  Sage fragment annotation enabled")
        else:
            print(f"Warning: Store not found: {p}")

    uvicorn.run(app, host=args.host, port=args.port)
