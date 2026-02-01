"""
FastAPI backend for Precursor Browser

Serves precursor data from Parquet store with efficient columnar reads.
"""

from pathlib import Path
from typing import Optional, List
import numpy as np
import pandas as pd

from fastapi import FastAPI, HTTPException, Query
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel
import pyarrow.parquet as pq

app = FastAPI(title="Precursor Browser API", version="1.0.0")

# CORS for local development
app.add_middleware(
    CORSMiddleware,
    allow_origins=["http://localhost:5173", "http://localhost:3000"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Global store reference
store_path: Optional[Path] = None
parquet_file: Optional[pq.ParquetFile] = None


class PrecursorSummary(BaseModel):
    precursor_id: int
    raw_file: str
    mz: float
    charge: int
    rt_seconds: float
    mobility: float
    n_engines: int
    confidence_weight: Optional[float] = None
    # FragPipe
    fragpipe_peptide: Optional[str] = None
    fragpipe_modified: Optional[str] = None
    fragpipe_probability: Optional[float] = None
    fragpipe_qvalue: Optional[float] = None
    # Sage
    sage_peptide: Optional[str] = None
    sage_modified: Optional[str] = None
    sage_qvalue: Optional[float] = None
    sage_match_tier: Optional[str] = None
    # DIA-NN
    diann_peptide: Optional[str] = None
    diann_modified: Optional[str] = None
    diann_qvalue: Optional[float] = None
    diann_ccs: Optional[float] = None
    diann_match_tier: Optional[str] = None
    diann_match_score: Optional[float] = None
    # Raw
    raw_intensity_meta: Optional[float] = None
    frame_id: Optional[int] = None
    isolation_mz: Optional[float] = None


class PrecursorDetail(BaseModel):
    precursor_id: int
    mz: float
    charge: int
    rt_seconds: float
    mobility: float
    n_engines: int
    fragpipe_peptide: Optional[str] = None
    sage_peptide: Optional[str] = None
    diann_peptide: Optional[str] = None

    # Fragment spectrum
    fragment_mz: List[float]
    fragment_intensity: List[float]
    fragment_mobility: List[float]
    fragment_scan: List[int]

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


@app.get("/")
async def root():
    return {"status": "ok", "store_loaded": parquet_file is not None}


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
    sort_by: str = Query("n_engines", pattern="^(raw_intensity_meta|n_engines|mz|rt_seconds|precursor_id)$"),
    sort_desc: bool = Query(True),
):
    """List precursors with pagination and filtering."""
    if parquet_file is None:
        raise HTTPException(400, "No store loaded")

    # Read only scalar columns for listing
    # Use raw_* columns from precursor index
    columns = [
        'precursor_id', 'raw_file', 'raw_mz', 'raw_charge', 'raw_rt_seconds', 'raw_mobility',
        'n_engines', 'confidence_weight',
        # FragPipe
        'fragpipe_peptide', 'fragpipe_modified', 'fragpipe_probability', 'fragpipe_qvalue',
        # Sage
        'sage_peptide', 'sage_modified', 'sage_qvalue', 'sage_match_tier',
        # DIA-NN
        'diann_peptide', 'diann_modified', 'diann_qvalue', 'diann_ccs', 'diann_match_tier', 'diann_match_score',
        # Raw
        'raw_intensity', 'frame_id', 'isolation_mz'
    ]

    # Build filters
    filters = []
    if min_engines > 0:
        filters.append(('n_engines', '>=', min_engines))
    if charge is not None:
        filters.append(('raw_charge', '=', charge))
    if raw_file is not None:
        filters.append(('raw_file', '=', raw_file))

    # Read with filters
    table = pq.read_table(str(store_path), columns=columns, filters=filters if filters else None)
    df = table.to_pandas()

    # Filter for precursors with MS1 data (check if raw_rt_seconds is not null)
    if has_ms1:
        df = df[df['raw_rt_seconds'].notna()]

    # Sort - when sorting by n_engines, also sort by intensity within each group
    if sort_by == 'n_engines' and 'raw_intensity' in df.columns:
        # Multi-column sort: engines first, then intensity
        df = df.sort_values(
            ['n_engines', 'raw_intensity'],
            ascending=[not sort_desc, False],  # engines: user choice, intensity: always descending
            na_position='last'
        )
    elif sort_by == 'mz' and 'raw_mz' in df.columns:
        df = df.sort_values('raw_mz', ascending=not sort_desc, na_position='last')
    elif sort_by == 'rt_seconds' and 'raw_rt_seconds' in df.columns:
        df = df.sort_values('raw_rt_seconds', ascending=not sort_desc, na_position='last')
    elif sort_by in df.columns:
        df = df.sort_values(sort_by, ascending=not sort_desc, na_position='last')

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
    for _, row in df.iterrows():
        results.append(PrecursorSummary(
            precursor_id=int(row['precursor_id']),
            raw_file=str(row['raw_file']) if row['raw_file'] else '',
            mz=float(row['raw_mz']) if pd.notna(row.get('raw_mz')) else 0.0,
            charge=int(row['raw_charge']) if pd.notna(row.get('raw_charge')) else 0,
            rt_seconds=float(row['raw_rt_seconds']) if pd.notna(row.get('raw_rt_seconds')) else 0.0,
            mobility=float(row['raw_mobility']) if pd.notna(row.get('raw_mobility')) else 0.0,
            n_engines=int(row['n_engines']),
            confidence_weight=safe_float(row.get('confidence_weight')),
            # FragPipe
            fragpipe_peptide=safe_str(row.get('fragpipe_peptide')),
            fragpipe_modified=safe_str(row.get('fragpipe_modified')),
            fragpipe_probability=safe_float(row.get('fragpipe_probability')),
            fragpipe_qvalue=safe_float(row.get('fragpipe_qvalue')),
            # Sage
            sage_peptide=safe_str(row.get('sage_peptide')),
            sage_modified=safe_str(row.get('sage_modified')),
            sage_qvalue=safe_float(row.get('sage_qvalue')),
            sage_match_tier=safe_str(row.get('sage_match_tier')),
            # DIA-NN
            diann_peptide=safe_str(row.get('diann_peptide')),
            diann_modified=safe_str(row.get('diann_modified')),
            diann_qvalue=safe_float(row.get('diann_qvalue')),
            diann_ccs=safe_float(row.get('diann_ccs')),
            diann_match_tier=safe_str(row.get('diann_match_tier')),
            diann_match_score=safe_float(row.get('diann_match_score')),
            # Raw
            raw_intensity_meta=safe_float(row.get('raw_intensity')),
            frame_id=safe_int(row.get('frame_id')),
            isolation_mz=safe_float(row.get('isolation_mz')),
        ))

    return results


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

    # Handle missing fragment_scan column for older stores
    fragment_scan = to_int_list(row.get('fragment_scan')) if 'fragment_scan' in row else []

    return PrecursorDetail(
        precursor_id=int(row['precursor_id']),
        mz=float(row['mz']),
        charge=int(row['charge']),
        rt_seconds=float(row['rt_seconds']),
        mobility=float(row['mobility']),
        n_engines=int(row['n_engines']),
        fragpipe_peptide=row['fragpipe_peptide'] if row['fragpipe_peptide'] else None,
        sage_peptide=row['sage_peptide'] if row['sage_peptide'] else None,
        diann_peptide=row['diann_peptide'] if row['diann_peptide'] else None,
        fragment_mz=to_list(row['fragment_mz']),
        fragment_intensity=to_list(row['fragment_intensity']),
        fragment_mobility=to_list(row['fragment_mobility']),
        fragment_scan=fragment_scan,
        xic_rt=to_list(row['xic_rt']),
        xic_intensity=to_list(row['xic_intensity']),
        mobilogram_im=to_list(row['mobilogram_im']),
        mobilogram_intensity=to_list(row['mobilogram_intensity']),
        isotope_mz=to_list(row['isotope_mz']),
        isotope_intensity=to_list(row['isotope_intensity']),
        raw_rt=to_list(row['raw_rt']),
        raw_mz=to_list(row['raw_mz']),
        raw_mobility=to_list(row['raw_mobility']),
        raw_intensity=to_list(row['raw_intensity']),
    )


@app.get("/stats")
async def get_stats():
    """Get summary statistics."""
    if parquet_file is None:
        raise HTTPException(400, "No store loaded")

    # Read minimal columns for stats
    table = pq.read_table(str(store_path), columns=['n_engines', 'raw_charge', 'raw_mz', 'raw_file'])
    df = table.to_pandas()

    return {
        "total_precursors": len(df),
        "by_engines": df['n_engines'].value_counts().to_dict(),
        "by_charge": df['raw_charge'].value_counts().to_dict(),
        "mz_range": [float(df['raw_mz'].min()), float(df['raw_mz'].max())],
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


if __name__ == "__main__":
    import uvicorn
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--store", help="Path to Parquet store to load on startup")
    parser.add_argument("--port", type=int, default=8000)
    parser.add_argument("--host", default="127.0.0.1")
    args = parser.parse_args()

    if args.store:
        p = Path(args.store)
        if p.exists():
            parquet_file = pq.ParquetFile(str(p))
            store_path = p
            print(f"Loaded store: {p} ({parquet_file.metadata.num_rows} precursors)")

    uvicorn.run(app, host=args.host, port=args.port)
