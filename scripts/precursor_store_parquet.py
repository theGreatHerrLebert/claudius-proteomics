#!/usr/bin/env python3
"""
Parquet-based Precursor Store

Training-optimized storage using Parquet with list columns for variable-length data.
Each row is one precursor with all associated signal data.

Benefits over Zarr ragged arrays:
- Row groups provide natural batching for training
- Columnar reads - only load columns you need
- Predicate pushdown - filter before reading
- Native PyArrow/Polars/DuckDB support

Storage layout:
  precursors.parquet
    - One row per precursor
    - Scalar columns: precursor_id, mz, charge, rt, mobility, peptides, etc.
    - List columns: fragment_mz[], fragment_intensity[], raw_rt[], etc.

Usage:
    # Create store
    PrecursorStoreParquet.create_from_index_and_raw(
        index_path="index.parquet",
        raw_data_path="data.d",
        output_path="precursors.parquet",
    )

    # Load for training
    import pyarrow.parquet as pq
    table = pq.read_table("precursors.parquet",
                          columns=['precursor_id', 'raw_rt', 'raw_mz', 'raw_mobility', 'raw_intensity'],
                          filters=[('charge', '=', 2)])

    # Iterate batches
    for batch in pq.ParquetFile("precursors.parquet").iter_batches(batch_size=1000):
        ...
"""

import gc
from pathlib import Path
from dataclasses import dataclass
from typing import Optional, List, Iterator
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq


@dataclass
class PrecursorData:
    """Complete data for a single precursor."""
    # Identifiers
    precursor_id: int
    raw_file: str

    # Properties
    mz: float
    charge: int
    rt_seconds: float
    mobility: float

    # Search engine IDs (may be None)
    fragpipe_peptide: Optional[str] = None
    diann_peptide: Optional[str] = None
    sage_peptide: Optional[str] = None
    consensus_peptide: Optional[str] = None
    n_engines: int = 0

    # Fragment spectrum (variable length)
    fragment_mz: Optional[np.ndarray] = None
    fragment_intensity: Optional[np.ndarray] = None
    fragment_mobility: Optional[np.ndarray] = None

    # MS1 signals - 1D projections (variable length)
    xic_rt: Optional[np.ndarray] = None
    xic_intensity: Optional[np.ndarray] = None
    mobilogram_im: Optional[np.ndarray] = None
    mobilogram_intensity: Optional[np.ndarray] = None
    isotope_mz: Optional[np.ndarray] = None
    isotope_intensity: Optional[np.ndarray] = None

    # Raw 4D MS1 data (variable length)
    raw_rt: Optional[np.ndarray] = None
    raw_mz: Optional[np.ndarray] = None
    raw_mobility: Optional[np.ndarray] = None
    raw_intensity: Optional[np.ndarray] = None

    # Quality metrics
    fragpipe_probability: Optional[float] = None
    diann_qvalue: Optional[float] = None
    sage_qvalue: Optional[float] = None

    # Fragment annotations (parallel arrays, same length as fragment_mz)
    # Added by build_training_spectra.py
    fragment_intensity_norm: Optional[np.ndarray] = None  # Normalized intensities
    fragment_ion_type: Optional[List[str]] = None         # "b", "y", or "other"
    fragment_ion_number: Optional[np.ndarray] = None      # Ion position (0=unmatched)
    fragment_ion_charge: Optional[np.ndarray] = None      # Fragment charge
    fragment_theoretical_mz: Optional[np.ndarray] = None  # Theoretical m/z
    fragment_error_ppm: Optional[np.ndarray] = None       # Mass error in ppm

    # Annotation summary metrics
    is_annotated: bool = False                 # Has sequence for annotation
    n_matched_peaks: int = 0                   # Peaks matched to theoretical
    intensity_explained: float = 0.0           # Fraction of intensity matched
    sequence_coverage_b: float = 0.0           # b ion series coverage
    sequence_coverage_y: float = 0.0           # y ion series coverage

    # Normalization metadata
    normalization_factor: Optional[float] = None    # For reversibility
    normalization_method: Optional[str] = None      # Method used


class PrecursorStoreParquet:
    """Parquet-based precursor store optimized for batch training."""

    def __init__(self, store_path: str):
        """Open existing store."""
        self.store_path = Path(store_path)
        self.parquet_file = pq.ParquetFile(str(self.store_path))

        # Load lightweight index (scalar columns only)
        self._load_index()

    def _load_index(self):
        """Load scalar columns as index DataFrame."""
        # Read only scalar columns for fast indexing
        scalar_cols = [
            'precursor_id', 'raw_file', 'mz', 'charge', 'rt_seconds', 'mobility',
            'fragpipe_peptide', 'diann_peptide', 'sage_peptide', 'consensus_peptide',
            'n_engines', 'fragpipe_probability', 'diann_qvalue', 'sage_qvalue',
            'raw_intensity_meta',  # metadata intensity from timsTOF
        ]
        # Filter to columns that exist
        available_cols = self.parquet_file.schema.names
        cols_to_read = [c for c in scalar_cols if c in available_cols]

        self.index_df = self.parquet_file.read(columns=cols_to_read).to_pandas()

        # Build precursor_id -> row index lookup
        self._id_to_idx = {
            pid: idx for idx, pid in enumerate(self.index_df['precursor_id'])
        }

    def __len__(self):
        return len(self.index_df)

    def get_precursor(self, precursor_id: int) -> Optional[PrecursorData]:
        """Get complete data for a precursor by ID."""
        if precursor_id not in self._id_to_idx:
            return None

        idx = self._id_to_idx[precursor_id]

        # Read single row with all columns
        # Use row group filtering for efficiency
        row_group_idx = idx // self.parquet_file.metadata.row_group(0).num_rows

        # Read the full row
        table = self.parquet_file.read_row_group(row_group_idx)
        local_idx = idx % self.parquet_file.metadata.row_group(0).num_rows

        row = table.slice(local_idx, 1).to_pandas().iloc[0]

        return self._row_to_precursor_data(row)

    def _row_to_precursor_data(self, row) -> PrecursorData:
        """Convert a DataFrame row to PrecursorData."""
        def to_array(val):
            if val is None or (isinstance(val, float) and np.isnan(val)):
                return None
            if isinstance(val, (list, np.ndarray)):
                arr = np.array(val)
                return arr if len(arr) > 0 else None
            return None

        return PrecursorData(
            precursor_id=int(row['precursor_id']),
            raw_file=str(row.get('raw_file', '')),
            mz=float(row.get('mz', 0)),
            charge=int(row.get('charge', 2)),
            rt_seconds=float(row.get('rt_seconds', 0)),
            mobility=float(row.get('mobility', 0)),
            fragpipe_peptide=row.get('fragpipe_peptide'),
            diann_peptide=row.get('diann_peptide'),
            sage_peptide=row.get('sage_peptide'),
            consensus_peptide=row.get('consensus_peptide'),
            n_engines=int(row.get('n_engines', 0)),
            fragment_mz=to_array(row.get('fragment_mz')),
            fragment_intensity=to_array(row.get('fragment_intensity')),
            fragment_mobility=to_array(row.get('fragment_mobility')),
            xic_rt=to_array(row.get('xic_rt')),
            xic_intensity=to_array(row.get('xic_intensity')),
            mobilogram_im=to_array(row.get('mobilogram_im')),
            mobilogram_intensity=to_array(row.get('mobilogram_intensity')),
            isotope_mz=to_array(row.get('isotope_mz')),
            isotope_intensity=to_array(row.get('isotope_intensity')),
            raw_rt=to_array(row.get('raw_rt')),
            raw_mz=to_array(row.get('raw_mz')),
            raw_mobility=to_array(row.get('raw_mobility')),
            raw_intensity=to_array(row.get('raw_intensity')),
            fragpipe_probability=row.get('fragpipe_probability'),
            diann_qvalue=row.get('diann_qvalue'),
            sage_qvalue=row.get('sage_qvalue'),
            # Annotation fields (added by build_training_spectra.py)
            fragment_intensity_norm=to_array(row.get('fragment_intensity_norm')),
            fragment_ion_type=row.get('fragment_ion_type') if isinstance(row.get('fragment_ion_type'), list) else None,
            fragment_ion_number=to_array(row.get('fragment_ion_number')),
            fragment_ion_charge=to_array(row.get('fragment_ion_charge')),
            fragment_theoretical_mz=to_array(row.get('fragment_theoretical_mz')),
            fragment_error_ppm=to_array(row.get('fragment_error_ppm')),
            is_annotated=bool(row.get('is_annotated', False)),
            n_matched_peaks=int(row.get('n_matched_peaks', 0)),
            intensity_explained=float(row.get('intensity_explained', 0.0)),
            sequence_coverage_b=float(row.get('sequence_coverage_b', 0.0)),
            sequence_coverage_y=float(row.get('sequence_coverage_y', 0.0)),
            normalization_factor=row.get('normalization_factor'),
            normalization_method=row.get('normalization_method'),
        )

    def get_by_agreement(self, min_engines: int = 1) -> pd.DataFrame:
        """Get precursors sorted by agreement (n_engines descending)."""
        df = self.index_df[self.index_df['n_engines'] >= min_engines].copy()
        sort_cols = ['n_engines']
        if 'fragpipe_probability' in df.columns:
            sort_cols.append('fragpipe_probability')
        return df.sort_values(sort_cols, ascending=[False] * len(sort_cols))

    def iter_batches(
        self,
        batch_size: int = 1000,
        columns: Optional[List[str]] = None,
        filters: Optional[List] = None,
    ) -> Iterator[pa.RecordBatch]:
        """
        Iterate over precursors in batches (training-optimized).

        Args:
            batch_size: Number of precursors per batch
            columns: Columns to load (None = all)
            filters: PyArrow filter expressions, e.g. [('charge', '=', 2)]

        Yields:
            PyArrow RecordBatch objects
        """
        # Use ParquetFile.iter_batches for memory efficiency
        for batch in self.parquet_file.iter_batches(
            batch_size=batch_size,
            columns=columns,
        ):
            yield batch

    def read_filtered(
        self,
        columns: Optional[List[str]] = None,
        filters: Optional[List] = None,
    ) -> pa.Table:
        """
        Read with column selection and row filtering.

        Args:
            columns: Columns to load
            filters: PyArrow filter expressions

        Returns:
            PyArrow Table
        """
        return pq.read_table(
            str(self.store_path),
            columns=columns,
            filters=filters,
        )

    @classmethod
    def create_from_index_and_raw(
        cls,
        index_path: str,
        raw_data_path: str,
        output_path: str,
        num_threads: int = 16,
        rt_window_sec: float = 30.0,
        mz_tol_ppm: float = 20.0,
        im_window: float = 0.1,
        batch_size: int = 2000,
        limit: Optional[int] = None,
        row_group_size: int = 10000,
        use_calibration: bool = True,
    ) -> 'PrecursorStoreParquet':
        """
        Create Parquet store from unified index + raw data extraction.
        Uses streaming writes to handle large datasets without OOM.

        Args:
            index_path: Path to precursor_index.parquet
            raw_data_path: Path to .d folder
            output_path: Where to create .parquet store
            num_threads: Threads for extraction
            rt_window_sec: RT window for MS1 extraction
            mz_tol_ppm: m/z tolerance for MS1
            im_window: Ion mobility window
            batch_size: Batch size for extraction and writing
            limit: Limit number of precursors (for testing)
            row_group_size: Parquet row group size
            use_calibration: Use pre-computed IM calibration for accurate values

        Returns:
            PrecursorStoreParquet instance
        """
        from imspy_core.timstof import TimsDatasetDDA
        from imspy_connector import py_dda
        from imspy_connector.py_dataset import PyTimsDataset

        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        # Load index
        print(f"Loading index: {index_path}")
        index_df = pd.read_parquet(index_path)

        # Filter to this raw file
        raw_name = Path(raw_data_path).stem
        file_index = index_df[index_df['raw_file'] == raw_name].copy().reset_index(drop=True)
        print(f"  {len(file_index)} precursors for {raw_name}")

        if len(file_index) == 0:
            raise ValueError(f"No precursors found for raw file: {raw_name}")

        # Apply limit for testing
        if limit is not None and limit < len(file_index):
            print(f"  Limiting to {limit} precursors")
            file_index = file_index.head(limit)

        # Load dataset with optional calibration
        print(f"Loading dataset: {raw_data_path}")
        raw_data_path_obj = Path(raw_data_path)

        if use_calibration:
            # Import calibration utilities
            import sys
            scripts_dir = Path(__file__).parent
            if str(scripts_dir) not in sys.path:
                sys.path.insert(0, str(scripts_dir))
            from extract_calibration import ensure_calibration, get_calibration_path

            # Check for cached calibration or extract it
            cal_path = get_calibration_path(str(raw_data_path))
            if cal_path.exists():
                print(f"  Using cached IM calibration: {cal_path}")
                im_lookup = np.load(cal_path)
            else:
                print(f"  Extracting IM calibration (one-time, using Bruker SDK)...")
                im_lookup = ensure_calibration(str(raw_data_path), verbose=True)

            # Create DDA dataset with calibration (fast + accurate + thread-safe)
            from imspy_connector.py_dda import PyTimsDatasetDDA as PyTimsDatasetDDARust
            rust_dataset = PyTimsDatasetDDARust.with_calibration(
                str(raw_data_path), False, im_lookup.tolist()
            )
            print(f"  Dataset loaded with calibrated IM (thread-safe parallel extraction)")
        else:
            # Fallback: use linear interpolation (fast but inaccurate)
            print(f"  WARNING: Using linear IM interpolation (may be inaccurate)")
            dataset = TimsDatasetDDA(str(raw_data_path), in_memory=False, use_bruker_sdk=False)
            rust_dataset = dataset.get_py_ptr()

        # Helper function to extract and merge fragments for a batch of precursor IDs
        def extract_fragments_for_batch(precursor_ids: List[int]) -> dict:
            """
            Extract fragments for specific precursor IDs using batched Rust API.
            Merges re-targeted precursors (same precursor fragmented multiple times).
            Returns dict: precursor_id -> {mz, intensity, mobility, scan, rt}
            """
            if not precursor_ids:
                return {}

            # Use the new batched API - only loads frames for requested precursors
            fragments_list = rust_dataset.get_pasef_fragments_for_precursors(
                [int(pid) for pid in precursor_ids],
                num_threads,
            )

            # Merge re-targeted precursors (multiple PASEF events for same precursor)
            fragment_lookup = {}
            for frag in fragments_list:
                pid = int(frag.precursor_id)
                frame = frag.selected_fragment

                if pid not in fragment_lookup:
                    fragment_lookup[pid] = {
                        'mz': list(frame.mz),
                        'intensity': list(frame.intensity),
                        'mobility': list(frame.mobility),
                        'scan': [int(s) for s in frame.scan],
                        'rt': frame.retention_time,
                    }
                else:
                    # Merge re-captured frames
                    existing = fragment_lookup[pid]
                    existing['mz'].extend(frame.mz)
                    existing['intensity'].extend(frame.intensity)
                    existing['mobility'].extend(frame.mobility)
                    existing['scan'].extend([int(s) for s in frame.scan])

            return fragment_lookup

        print("Using batched fragment extraction (memory-efficient)...")

        # Create PyArrow schema
        schema = pa.schema([
            ('precursor_id', pa.int64()),
            ('raw_file', pa.string()),
            ('mz', pa.float64()),
            ('charge', pa.int32()),
            ('rt_seconds', pa.float64()),
            ('mobility', pa.float64()),

            ('fragpipe_peptide', pa.string()),
            ('fragpipe_modified', pa.string()),
            ('diann_peptide', pa.string()),
            ('diann_modified', pa.string()),
            ('sage_peptide', pa.string()),
            ('sage_modified', pa.string()),
            ('consensus_peptide', pa.string()),
            ('n_engines', pa.int32()),

            ('fragpipe_probability', pa.float64()),
            ('diann_qvalue', pa.float64()),
            ('sage_qvalue', pa.float64()),
            ('raw_intensity_meta', pa.float64()),

            ('fragment_mz', pa.list_(pa.float64())),
            ('fragment_intensity', pa.list_(pa.float64())),
            ('fragment_mobility', pa.list_(pa.float64())),
            ('fragment_scan', pa.list_(pa.int32())),  # Scan indices for heatmap visualization

            ('xic_rt', pa.list_(pa.float64())),
            ('xic_intensity', pa.list_(pa.float64())),
            ('mobilogram_im', pa.list_(pa.float64())),
            ('mobilogram_intensity', pa.list_(pa.float64())),
            ('isotope_mz', pa.list_(pa.float64())),
            ('isotope_intensity', pa.list_(pa.float64())),

            ('raw_rt', pa.list_(pa.float64())),
            ('raw_mz', pa.list_(pa.float64())),
            ('raw_mobility', pa.list_(pa.float64())),
            ('raw_intensity', pa.list_(pa.float64())),
        ])

        # Process in batches and write incrementally
        n_precursors = len(file_index)
        n_batches = (n_precursors + batch_size - 1) // batch_size

        # Stats accumulators
        total_frag_peaks = 0
        total_xic_points = 0
        total_raw_points = 0

        print(f"Processing {n_precursors} precursors in {n_batches} batches...")

        # Use ParquetWriter for streaming writes
        writer = None

        for batch_idx in range(n_batches):
            batch_start = batch_idx * batch_size
            batch_end = min(batch_start + batch_size, n_precursors)
            batch_df = file_index.iloc[batch_start:batch_end]

            print(f"  Batch {batch_idx + 1}/{n_batches}: precursors {batch_start}-{batch_end}...")

            # Extract fragments for this batch using the batched Rust API
            batch_precursor_ids = batch_df['precursor_id'].tolist()
            fragment_lookup = extract_fragments_for_batch(batch_precursor_ids)

            # Build coords for MS1 extraction
            coords = []
            batch_meta = {}

            for _, row in batch_df.iterrows():
                pid = int(row['precursor_id'])

                # Get RT from fragments if available (RT is already in seconds)
                frag_data = fragment_lookup.get(pid, {})
                if frag_data and frag_data.get('rt'):
                    rt_sec = float(frag_data.get('rt', 0))
                else:
                    # FragPipe RT is also in seconds
                    rt_val = row.get('fragpipe_rt')
                    rt_sec = float(rt_val) if pd.notna(rt_val) else 0.0

                # Get precursor properties - use pd.notna() to handle NaN properly
                # (Python's `or` doesn't work correctly with NaN - np.nan is truthy)
                mz_fp = row.get('fragpipe_mz')
                mz_raw = row.get('raw_mz')
                mz = float(mz_fp) if pd.notna(mz_fp) else (float(mz_raw) if pd.notna(mz_raw) else 0.0)

                charge_fp = row.get('fragpipe_charge')
                charge_raw = row.get('raw_charge')
                charge_val = charge_fp if pd.notna(charge_fp) else (charge_raw if pd.notna(charge_raw) else 2)
                charge = int(charge_val) if pd.notna(charge_val) else 2

                mobility_fp = row.get('fragpipe_mobility')
                mobility_raw = row.get('raw_mobility')
                mobility = float(mobility_fp) if pd.notna(mobility_fp) else (float(mobility_raw) if pd.notna(mobility_raw) else 0.9)

                if mz > 0:
                    coord = py_dda.PyPrecursorCoord(
                        precursor_id=pid,
                        mz=float(mz),
                        rt_seconds=float(rt_sec),
                        mobility=float(mobility),
                        charge=charge,
                    )
                    coords.append(coord)

                batch_meta[pid] = {
                    'rt_seconds': rt_sec,
                    'mz': mz,
                    'charge': charge,
                    'mobility': mobility,
                }

            # Extract MS1 signals for this batch
            ms1_lookup = {}
            if coords:
                batch_signals = rust_dataset.extract_precursor_ms1_signals(
                    coords,
                    rt_window_sec=rt_window_sec,
                    mz_tol_ppm=mz_tol_ppm,
                    im_window=im_window,
                    n_isotopes=5,
                    num_threads=num_threads,
                )

                for sig in batch_signals:
                    ms1_lookup[sig.precursor_id] = {
                        'rt_coords': list(sig.rt_coords),
                        'rt_intensities': list(sig.rt_intensities),
                        'im_coords': list(sig.im_coords),
                        'im_intensities': list(sig.im_intensities),
                        'isotope_mz': list(sig.isotope_mz),
                        'isotope_intensity': list(sig.isotope_intensity),
                        'raw_rt': list(sig.raw_rt) if hasattr(sig, 'raw_rt') else [],
                        'raw_mz': list(sig.raw_mz) if hasattr(sig, 'raw_mz') else [],
                        'raw_mobility': list(sig.raw_mobility) if hasattr(sig, 'raw_mobility') else [],
                        'raw_intensity': list(sig.raw_intensity) if hasattr(sig, 'raw_intensity') else [],
                    }
                del batch_signals

            # Build rows for this batch
            rows = []
            for _, row in batch_df.iterrows():
                pid = int(row['precursor_id'])
                meta = batch_meta.get(pid, {})
                ms1 = ms1_lookup.get(pid, {})
                frag = fragment_lookup.get(pid, {})

                frag_mz = frag.get('mz', [])
                frag_int = frag.get('intensity', [])
                frag_mob = frag.get('mobility', [])
                frag_scan = frag.get('scan', [])

                total_frag_peaks += len(frag_mz)
                total_xic_points += len(ms1.get('rt_coords', []))
                total_raw_points += len(ms1.get('raw_rt', []))

                precursor_row = {
                    'precursor_id': pid,
                    'raw_file': str(row.get('raw_file', '')),
                    'mz': float(meta.get('mz', 0)),
                    'charge': int(meta.get('charge', 2)),
                    'rt_seconds': float(meta.get('rt_seconds', 0)),
                    'mobility': float(meta.get('mobility', 0)),

                    'fragpipe_peptide': row.get('fragpipe_peptide') if pd.notna(row.get('fragpipe_peptide')) else None,
                    'fragpipe_modified': row.get('fragpipe_modified') if pd.notna(row.get('fragpipe_modified')) else None,
                    'diann_peptide': row.get('diann_peptide') if pd.notna(row.get('diann_peptide')) else None,
                    'diann_modified': row.get('diann_modified') if pd.notna(row.get('diann_modified')) else None,
                    'sage_peptide': row.get('sage_peptide') if pd.notna(row.get('sage_peptide')) else None,
                    'sage_modified': row.get('sage_modified') if pd.notna(row.get('sage_modified')) else None,
                    'consensus_peptide': row.get('consensus_peptide') if pd.notna(row.get('consensus_peptide')) else None,
                    'n_engines': int(row.get('n_engines', 0)),

                    'fragpipe_probability': float(row['fragpipe_probability']) if pd.notna(row.get('fragpipe_probability')) else None,
                    'diann_qvalue': float(row['diann_qvalue']) if pd.notna(row.get('diann_qvalue')) else None,
                    'sage_qvalue': float(row['sage_qvalue']) if pd.notna(row.get('sage_qvalue')) else None,
                    'raw_intensity_meta': float(row['raw_intensity']) if pd.notna(row.get('raw_intensity')) else None,

                    'fragment_mz': frag_mz,
                    'fragment_intensity': frag_int,
                    'fragment_mobility': frag_mob,
                    'fragment_scan': frag_scan,

                    'xic_rt': ms1.get('rt_coords', []),
                    'xic_intensity': ms1.get('rt_intensities', []),
                    'mobilogram_im': ms1.get('im_coords', []),
                    'mobilogram_intensity': ms1.get('im_intensities', []),
                    'isotope_mz': ms1.get('isotope_mz', []),
                    'isotope_intensity': ms1.get('isotope_intensity', []),

                    'raw_rt': ms1.get('raw_rt', []),
                    'raw_mz': ms1.get('raw_mz', []),
                    'raw_mobility': ms1.get('raw_mobility', []),
                    'raw_intensity': ms1.get('raw_intensity', []),
                }
                rows.append(precursor_row)

            # Convert batch to table and write
            batch_table = pa.Table.from_pylist(rows, schema=schema)

            if writer is None:
                writer = pq.ParquetWriter(
                    str(output_path),
                    schema,
                    compression='zstd',
                    compression_level=3,
                    use_dictionary=True,
                    write_statistics=True,
                )

            writer.write_table(batch_table)

            # Clean up batch data
            del rows, batch_table, ms1_lookup, batch_meta, coords, fragment_lookup
            gc.collect()

        # Close writer (outside the for loop)
        if writer:
            writer.close()

        print(f"Store created: {output_path}")
        print(f"  Precursors: {n_precursors:,}")
        print(f"  Fragment peaks: {total_frag_peaks:,}")
        print(f"  XIC points: {total_xic_points:,}")
        print(f"  Raw 4D points: {total_raw_points:,}")
        print(f"  File size: {output_path.stat().st_size / 1e6:.1f} MB")

        return cls(str(output_path))


def main():
    """CLI for creating stores."""
    import argparse

    parser = argparse.ArgumentParser(description="Create Parquet precursor store")
    parser.add_argument("--index", required=True, help="Path to precursor_index.parquet")
    parser.add_argument("--raw", required=True, help="Path to .d folder")
    parser.add_argument("--output", required=True, help="Output .parquet path")
    parser.add_argument("--threads", type=int, default=16)
    parser.add_argument("--rt-window", type=float, default=30.0)
    parser.add_argument("--batch-size", type=int, default=2000)
    parser.add_argument("--limit", type=int, default=None, help="Limit precursors (testing)")
    parser.add_argument("--row-group-size", type=int, default=10000)
    parser.add_argument("--no-calibration", action="store_true",
                        help="Skip IM calibration (use linear interpolation - less accurate)")

    args = parser.parse_args()

    PrecursorStoreParquet.create_from_index_and_raw(
        index_path=args.index,
        raw_data_path=args.raw,
        output_path=args.output,
        num_threads=args.threads,
        rt_window_sec=args.rt_window,
        batch_size=args.batch_size,
        limit=args.limit,
        row_group_size=args.row_group_size,
        use_calibration=not args.no_calibration,
    )


if __name__ == "__main__":
    main()
