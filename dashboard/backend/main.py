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
from typing import Optional, List, Dict, Any
import numpy as np
import pandas as pd

from fastapi import FastAPI, HTTPException, Query
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel
import pyarrow.parquet as pq
import yaml

import gzip
import re


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

# Sage fragment data (loaded from engines/sage/)
sage_fragments_by_psm: Dict[int, List[Dict]] = {}  # psm_id -> list of fragment dicts
sage_peptide_to_psm: Dict[str, int] = {}  # "peptide_unimod_charge" -> psm_id (fallback lookup)
sage_data_loaded: bool = False

# Global collection reference (collection mode)
collection_path: Optional[Path] = None
collection_manifest: Optional[Dict[str, Any]] = None
studies_config: Optional[Dict[str, Any]] = None
active_dataset: Optional[str] = None  # Currently loaded dataset accession


def normalize_sage_to_unimod(peptide: str) -> str:
    """Convert Sage mass delta format [+57.021465] to UNIMOD format [UNIMOD:4]."""
    if not peptide:
        return peptide
    result = peptide
    result = re.sub(r'\[\+57\.02\d*\]', '[UNIMOD:4]', result)   # Carbamidomethyl
    result = re.sub(r'\[\+15\.99\d*\]', '[UNIMOD:35]', result)  # Oxidation
    result = re.sub(r'\[\+42\.01\d*\]', '[UNIMOD:1]', result)   # Acetyl
    result = re.sub(r'\[\+79\.96\d*\]', '[UNIMOD:21]', result)  # Phospho
    return result


def load_sage_fragments(sage_dir: Path) -> bool:
    """Load Sage matched fragments for fragment annotation.

    Args:
        sage_dir: Path to directory containing Sage output files

    Looks for:
    - matched_fragments.sage.parquet
    - results.sage.parquet (for fallback peptide+charge lookup)

    Builds lookup dicts:
    - sage_fragments_by_psm: psm_id -> list of fragment dicts
    - sage_peptide_to_psm: "peptide_unimod_charge" -> psm_id (fallback)

    The precursor_store.parquet should have sage_psm_id column linking
    to these fragments. If not (older data), falls back to peptide+charge matching.
    """
    global sage_fragments_by_psm, sage_peptide_to_psm, sage_data_loaded

    fragments_path = sage_dir / "matched_fragments.sage.parquet"
    results_path = sage_dir / "results.sage.parquet"

    if not fragments_path.exists():
        print(f"  Sage fragments not found at {fragments_path}")
        return False

    try:
        # Load fragments
        print(f"  Loading Sage fragments from {fragments_path}")
        fragments_df = pq.read_table(str(fragments_path)).to_pandas()

        # Group fragments by psm_id
        sage_fragments_by_psm = {}
        for psm_id, group in fragments_df.groupby('psm_id'):
            fragments = []
            for _, row in group.iterrows():
                fragments.append({
                    'fragment_type': row['fragment_type'],
                    'ion_number': int(row['fragment_ordinals']),
                    'charge': int(row['fragment_charge']),
                    'mz_experimental': float(row['fragment_mz_experimental']),
                    'mz_calculated': float(row['fragment_mz_calculated']),
                    'intensity': float(row['fragment_intensity']),
                })
            sage_fragments_by_psm[psm_id] = fragments
        print(f"    Loaded fragments for {len(sage_fragments_by_psm)} PSMs")

        # Build fallback peptide+charge -> psm_id lookup from results
        if results_path.exists():
            print(f"  Building fallback peptide lookup from {results_path}")
            results_df = pq.read_table(str(results_path)).to_pandas()

            # Normalize peptide format to UNIMOD
            results_df['peptide_norm'] = results_df['peptide'].apply(normalize_sage_to_unimod)

            # Create key: peptide_norm + "_" + charge
            results_df['lookup_key'] = results_df['peptide_norm'] + '_' + results_df['charge'].astype(str)

            # For each unique key, take the PSM with lowest posterior_error (best match)
            best_psms = results_df.sort_values('posterior_error').drop_duplicates('lookup_key', keep='first')
            sage_peptide_to_psm = dict(zip(best_psms['lookup_key'], best_psms['psm_id']))
            print(f"    Built fallback lookup for {len(sage_peptide_to_psm)} peptide+charge combos")

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


class PrecursorSummary(BaseModel):
    precursor_id: int
    raw_file: str
    mz: float
    charge: int
    rt_seconds: float
    mobility: float
    n_engines: int
    confidence_weight: Optional[float] = None
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
    sage_peptide: Optional[str] = None
    sage_modified: Optional[str] = None
    diann_peptide: Optional[str] = None

    # Fragment spectrum
    fragment_mz: List[float]
    fragment_intensity: List[float]
    fragment_mobility: List[float]
    fragment_scan: List[int]

    # Sage matched b/y ions (if available)
    sage_matched_fragments: Optional[List[SageMatchedFragment]] = None

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
    charge: Optional[int] = Query(None),
    raw_file: Optional[str] = Query(None, description="Filter by raw file name"),
    has_ms1: bool = Query(False, description="Only show precursors with MS1 signal data"),
    sort_by: str = Query("n_engines", pattern="^(raw_intensity_meta|precursor_intensity|n_engines|mz|rt_seconds|precursor_id|mobility|fragpipe_probability|fragpipe_hyperscore|sage_hyperscore|sage_qvalue|diann_qvalue|ms1_rt_r2|ms1_im_r2|isotope_cosim)$"),
    sort_desc: bool = Query(True),
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
        'sage_rt', 'sage_mz', 'sage_mobility', 'sage_match_tier',
        # DIA-NN (13 columns)
        'diann_peptide', 'diann_modified', 'diann_protein',
        'diann_qvalue', 'diann_pep', 'diann_global_qvalue', 'diann_pg_qvalue',
        'diann_rt', 'diann_mz', 'diann_mobility', 'diann_ccs', 'diann_match_tier', 'diann_match_score',
        # Raw (intensity_col is auto-detected: precursor_intensity, ms1_total_intensity, or raw_intensity_meta)
        intensity_col, 'frame_id', 'isolation_mz',
        # Quality metrics
        'ms1_rt_sigma', 'ms1_rt_r2', 'ms1_im_sigma', 'ms1_im_r2', 'isotope_cosim'
    ]

    # Get available columns from schema
    available_cols = parquet_file.schema_arrow.names
    columns = base_columns + [c for c in optional_columns if c in available_cols]

    # Build filters
    filters = []
    if min_engines > 0:
        filters.append(('n_engines', '>=', min_engines))
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

    if sort_by == 'n_engines' and intensity_col in df.columns:
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

        results.append(PrecursorSummary(
            precursor_id=int(prec_id),
            raw_file=str(row['raw_file']) if pd.notna(row.get('raw_file')) else '',
            mz=float(mz_val) if pd.notna(mz_val) else 0.0,
            charge=int(charge_val),
            rt_seconds=float(row[rt_col]) if pd.notna(row.get(rt_col)) else 0.0,
            mobility=float(row[mobility_col]) if pd.notna(row.get(mobility_col)) else 0.0,
            n_engines=int(n_engines_val),
            confidence_weight=safe_float(row.get('confidence_weight')),
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
            f.seek(offset)
            compressed = f.read(size)

        # Decompress with gzip
        combined = gzip.decompress(compressed)

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
        }

        return result

    except Exception as e:
        print(f"Error reading blob: {e}")
        return None


@app.get("/precursor/{precursor_id}", response_model=PrecursorDetail)
async def get_precursor(precursor_id: int):
    """Get full detail for a single precursor."""
    if parquet_file is None:
        raise HTTPException(400, "No store loaded")

    # Read with filter
    table = pq.read_table(str(store_path), filters=[('precursor_id', '=', precursor_id)])

    if len(table) == 0:
        raise HTTPException(404, f"Precursor not found: {precursor_id}")

    df = table.to_pandas()
    row = df.iloc[0]

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

    # Try to read blob data if available
    blob_data = None
    if pd.notna(row.get('blob_offset')) and pd.notna(row.get('blob_size')) and pd.notna(row.get('raw_file')):
        blob_data = read_blob(
            str(row['raw_file']),
            int(row['blob_offset']),
            int(row['blob_size'])
        )

    # Get precursor properties for isotope m/z calculation
    precursor_mz = float(row['mz']) if pd.notna(row.get('mz')) else 0.0
    precursor_charge = int(row['charge']) if pd.notna(row.get('charge')) else 1

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

    # Get Sage matched fragments if available
    sage_matched_fragments = None
    sage_modified = safe_str(row.get('sage_modified'))
    sage_psm_id = row.get('sage_psm_id')
    if sage_data_loaded:
        psm_id = None
        # Try direct psm_id lookup first (new data format)
        if pd.notna(sage_psm_id):
            psm_id = int(sage_psm_id)
        # Fallback to peptide+charge lookup (old data format)
        elif sage_modified and sage_peptide_to_psm:
            lookup_key = f"{sage_modified}_{precursor_charge}"
            psm_id = sage_peptide_to_psm.get(lookup_key)

        if psm_id is not None and psm_id in sage_fragments_by_psm:
            sage_matched_fragments = [
                SageMatchedFragment(**frag) for frag in sage_fragments_by_psm[psm_id]
            ]

    return PrecursorDetail(
        precursor_id=int(row['precursor_id']),
        mz=float(row['mz']) if pd.notna(row.get('mz')) else 0.0,
        charge=int(row['charge']) if pd.notna(row.get('charge')) else 0,
        rt_seconds=float(row['rt_seconds']) if pd.notna(row.get('rt_seconds')) else 0.0,
        mobility=float(row['mobility']) if pd.notna(row.get('mobility')) else 0.0,
        n_engines=int(row['n_engines']) if pd.notna(row.get('n_engines')) else 0,
        fragpipe_peptide=safe_str(row.get('fragpipe_peptide')),
        sage_peptide=safe_str(row.get('sage_peptide')),
        sage_modified=sage_modified,
        diann_peptide=safe_str(row.get('diann_peptide')),
        fragment_mz=fragment_mz,
        fragment_intensity=fragment_intensity,
        fragment_mobility=fragment_mobility,
        fragment_scan=fragment_scan,
        sage_matched_fragments=sage_matched_fragments,
        xic_rt=xic_rt,
        xic_intensity=xic_intensity,
        mobilogram_im=mobilogram_im,
        mobilogram_intensity=mobilogram_intensity,
        isotope_mz=isotope_mz,
        isotope_intensity=isotope_intensity,
        raw_rt=[],  # Not stored in blob format
        raw_mz=[],
        raw_mobility=[],
        raw_intensity=[],
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


@app.post("/datasets/{accession}/load")
async def load_dataset(accession: str):
    """Load a dataset from the collection (collection mode only)."""
    global store_path, parquet_file, active_dataset

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
    dataset_path = collection_path / dataset["path"] / "precursor_store.parquet"

    if not dataset_path.exists():
        raise HTTPException(404, f"Precursor store not found: {dataset_path}")

    try:
        parquet_file = pq.ParquetFile(str(dataset_path))
        store_path = dataset_path
        active_dataset = accession

        return {
            "status": "loaded",
            "accession": accession,
            "path": str(dataset_path),
            "num_precursors": parquet_file.metadata.num_rows,
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
    global collection_path, collection_manifest, studies_config

    collection_path = path

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
                print("  Blob reading enabled (gzip format)")
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
