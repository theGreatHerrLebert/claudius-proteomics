"""PrecursorStore: PyArrow query engine with predicate pushdown.

No in-memory index — uses PyArrow predicate pushdown for all queries.
Startup is instant regardless of dataset size.
"""

from pathlib import Path
from typing import Optional, List, Iterator, Dict, Any

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

from .models import Precursor, RawSignal
from .blob import BlobReader
from .filters import build_filters, ALL_SCALAR_COLUMNS, SUMMARY_COLUMNS, BLOB_COLUMNS


class PrecursorStore:
    """Query engine for precursor_store.parquet with lazy blob loading."""

    def __init__(self, store_path: str, blob_reader: Optional[BlobReader] = None):
        """Open a precursor store.

        Args:
            store_path: Path to precursor_store.parquet
            blob_reader: Optional BlobReader for raw signal access
        """
        self.store_path = Path(store_path)
        self.parquet_file = pq.ParquetFile(str(self.store_path))
        self.blob_reader = blob_reader

        self._schema = self.parquet_file.schema_arrow
        self._columns = self._schema.names
        self._num_rows = self.parquet_file.metadata.num_rows

    @property
    def num_precursors(self) -> int:
        return self._num_rows

    @property
    def columns(self) -> List[str]:
        return list(self._columns)

    def _available(self, cols: List[str]) -> List[str]:
        """Filter column list to those present in schema."""
        return [c for c in cols if c in self._columns]

    def _row_to_precursor(self, row: pd.Series) -> Precursor:
        """Convert a DataFrame row to a Precursor object."""
        def safe_float(val):
            if val is None or (isinstance(val, float) and np.isnan(val)):
                return None
            try:
                return float(val)
            except (TypeError, ValueError):
                return None

        def safe_int(val, default=0):
            if val is None:
                return default
            try:
                f = float(val)
                if np.isnan(f):
                    return default
                return int(f)
            except (TypeError, ValueError):
                return default

        def safe_str(val):
            if val is None or (isinstance(val, float) and np.isnan(val)):
                return None
            return str(val)

        def safe_bool(val):
            if val is None or (isinstance(val, float) and np.isnan(val)):
                return None
            return bool(val)

        blob_offset = safe_float(row.get('blob_offset'))
        blob_size = safe_float(row.get('blob_size'))

        return Precursor(
            precursor_id=safe_int(row.get('precursor_id'), 0),
            raw_file=safe_str(row.get('raw_file')) or '',
            mono_mz=float(safe_float(row.get('mono_mz')) or safe_float(row.get('mz')) or 0),
            charge=safe_int(row.get('charge'), 2),
            rt_seconds=float(safe_float(row.get('rt_seconds')) or 0),
            mobility=float(safe_float(row.get('mobility')) or 0),
            n_engines=safe_int(row.get('n_engines'), 0),
            consensus_peptide=safe_str(row.get('consensus_peptide')),
            fragpipe_peptide=safe_str(row.get('fragpipe_peptide')),
            fragpipe_modified=safe_str(row.get('fragpipe_modified')),
            diann_peptide=safe_str(row.get('diann_peptide')),
            diann_modified=safe_str(row.get('diann_modified')),
            sage_peptide=safe_str(row.get('sage_peptide')),
            sage_modified=safe_str(row.get('sage_modified')),
            is_high_quality=safe_bool(row.get('is_high_quality')),
            ms1_rt_r2=safe_float(row.get('ms1_rt_r2')),
            ms1_im_r2=safe_float(row.get('ms1_im_r2')),
            isotope_cosim=safe_float(row.get('isotope_cosim')),
            fragpipe_probability=safe_float(row.get('fragpipe_probability')),
            diann_qvalue=safe_float(row.get('diann_qvalue')),
            sage_qvalue=safe_float(row.get('sage_qvalue')),
            _blob_offset=int(blob_offset) if blob_offset is not None else None,
            _blob_size=int(blob_size) if blob_size is not None else None,
            _blob_reader=self.blob_reader,
        )

    def get(
        self, precursor_id: int, raw_file: Optional[str] = None,
    ) -> Optional[Precursor]:
        """Get precursor(s) by ID.

        Since precursor_id is not globally unique (it's per-raw-file),
        providing raw_file is recommended.

        Args:
            precursor_id: Bruker precursor ID
            raw_file: Raw file name (optional but recommended)

        Returns:
            Single Precursor if raw_file given, or first match.
            None if not found.
        """
        filters = [('precursor_id', '=', float(precursor_id))]
        if raw_file is not None:
            filters.append(('raw_file', '=', raw_file))

        table = pq.read_table(
            str(self.store_path),
            columns=self._available(ALL_SCALAR_COLUMNS),
            filters=filters,
        )

        if len(table) == 0:
            return None

        df = table.to_pandas()
        return self._row_to_precursor(df.iloc[0])

    def query(
        self,
        *,
        columns: Optional[List[str]] = None,
        limit: Optional[int] = None,
        **filter_kwargs,
    ) -> pd.DataFrame:
        """Query precursors with predicate pushdown.

        Returns scalar columns only (no blob data).

        Args:
            columns: Specific columns to return (default: SUMMARY_COLUMNS)
            limit: Max rows to return
            **filter_kwargs: Passed to build_filters() — sequence, charge,
                min_engines, raw_file, is_high_quality, min_rt_r2, etc.

        Returns:
            pandas DataFrame
        """
        if columns is None:
            columns = SUMMARY_COLUMNS
        columns = self._available(columns)

        filters = build_filters(**filter_kwargs)

        table = pq.read_table(
            str(self.store_path),
            columns=columns,
            filters=filters,
        )

        df = table.to_pandas()

        if limit is not None:
            df = df.head(limit)

        return df

    def iter(
        self,
        *,
        limit: Optional[int] = None,
        **filter_kwargs,
    ) -> Iterator[Precursor]:
        """Iterate over precursors as Precursor objects with lazy blob loading.

        Args:
            limit: Max precursors to yield
            **filter_kwargs: Passed to build_filters()

        Yields:
            Precursor objects
        """
        filters = build_filters(**filter_kwargs)
        columns = self._available(ALL_SCALAR_COLUMNS)

        table = pq.read_table(
            str(self.store_path),
            columns=columns,
            filters=filters,
        )
        df = table.to_pandas()

        if limit is not None:
            df = df.head(limit)

        for _, row in df.iterrows():
            yield self._row_to_precursor(row)

    def batches(
        self,
        batch_size: int = 1000,
        *,
        columns: Optional[List[str]] = None,
        include_signal: bool = False,
        **filter_kwargs,
    ) -> Iterator[pd.DataFrame]:
        """Streaming batch iteration for ML training.

        Args:
            batch_size: Rows per batch
            columns: Columns to include (default: ALL_SCALAR_COLUMNS)
            include_signal: If True, expand blob data into signal columns
                (slower — requires I/O per batch)
            **filter_kwargs: Passed to build_filters()

        Yields:
            pandas DataFrames of batch_size rows
        """
        if columns is None:
            columns = ALL_SCALAR_COLUMNS
        if include_signal:
            # Ensure blob columns are present for expansion
            for col in BLOB_COLUMNS:
                if col not in columns:
                    columns = list(columns) + [col]

        columns = self._available(columns)
        filters = build_filters(**filter_kwargs)

        # Use iter_batches for memory efficiency
        for batch in self.parquet_file.iter_batches(
            batch_size=batch_size,
            columns=columns,
        ):
            df = batch.to_pandas()

            # Apply post-hoc filters (iter_batches doesn't support predicate pushdown)
            if filters:
                df = self._apply_filters(df, filters)

            if len(df) == 0:
                continue

            if include_signal and self.blob_reader is not None:
                df = self._expand_blobs(df)

            yield df

    def _apply_filters(
        self, df: pd.DataFrame, filters: List
    ) -> pd.DataFrame:
        """Apply filter tuples to a DataFrame (post-hoc for iter_batches)."""
        mask = pd.Series(True, index=df.index)
        for col, op, val in filters:
            if col not in df.columns:
                continue
            if op == '=':
                mask &= df[col] == val
            elif op == '>=':
                mask &= df[col] >= val
            elif op == '<=':
                mask &= df[col] <= val
            elif op == '>':
                mask &= df[col] > val
            elif op == '<':
                mask &= df[col] < val
            elif op == '!=':
                mask &= df[col] != val
        return df[mask]

    def _expand_blobs(self, df: pd.DataFrame) -> pd.DataFrame:
        """Expand blob references into signal columns.

        Groups by raw_file for efficient batch reads.
        """
        if self.blob_reader is None:
            return df

        # Initialize signal columns
        frag_mz = [None] * len(df)
        frag_int = [None] * len(df)
        raw_rt = [None] * len(df)
        raw_mz = [None] * len(df)
        raw_mob = [None] * len(df)
        raw_intensity = [None] * len(df)

        # Group by raw_file for batch reads
        for raw_file, group in df.groupby('raw_file'):
            valid = group.dropna(subset=['blob_offset', 'blob_size'])
            if len(valid) == 0:
                continue

            offsets = valid['blob_offset'].astype(int).tolist()
            sizes = valid['blob_size'].astype(int).tolist()
            signals = self.blob_reader.read_batch(raw_file, offsets, sizes)

            for (idx, _), sig in zip(valid.iterrows(), signals):
                pos = df.index.get_loc(idx)
                if sig is not None:
                    frag_mz[pos] = sig.fragment_spectrum.mz
                    frag_int[pos] = sig.fragment_spectrum.intensity
                    raw_rt[pos] = sig.raw_point_cloud.rt
                    raw_mz[pos] = sig.raw_point_cloud.mz
                    raw_mob[pos] = sig.raw_point_cloud.mobility
                    raw_intensity[pos] = sig.raw_point_cloud.intensity

        df = df.copy()
        df['fragment_mz'] = frag_mz
        df['fragment_intensity'] = frag_int
        df['raw_rt'] = raw_rt
        df['raw_mz'] = raw_mz
        df['raw_mobility'] = raw_mob
        df['raw_intensity'] = raw_intensity
        return df

    def raw_files(self) -> List[str]:
        """List unique raw file names."""
        table = pq.read_table(
            str(self.store_path), columns=['raw_file'],
        )
        return table.column('raw_file').to_pylist()

    def raw_files_unique(self) -> List[str]:
        """List unique raw file names (deduplicated)."""
        return sorted(f for f in set(self.raw_files()) if f is not None)

    def summary(self) -> Dict[str, Any]:
        """Quick dataset summary without full scan."""
        raw_files = self.raw_files_unique()

        # Read lightweight columns for stats
        table = pq.read_table(
            str(self.store_path),
            columns=self._available(['n_engines', 'is_high_quality']),
        )
        df = table.to_pandas()

        n_engines = df['n_engines']
        return {
            'num_precursors': self._num_rows,
            'num_raw_files': len(raw_files),
            'raw_files': raw_files,
            'has_blobs': self.blob_reader is not None,
            'engine_agreement': {
                'n_0': int((n_engines == 0).sum()),
                'n_1': int((n_engines == 1).sum()),
                'n_2': int((n_engines == 2).sum()),
                'n_3': int((n_engines == 3).sum()),
            },
            'n_high_quality': int(df['is_high_quality'].sum()) if 'is_high_quality' in df else 0,
        }

    def __len__(self) -> int:
        return self._num_rows

    def __repr__(self) -> str:
        blobs = ", with blobs" if self.blob_reader else ""
        return f"<PrecursorStore: {self._num_rows:,} precursors{blobs}>"
