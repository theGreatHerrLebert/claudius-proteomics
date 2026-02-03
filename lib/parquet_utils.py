#!/usr/bin/env python3
"""
Parquet Read/Write Utilities

Helpers for reading and writing Parquet files with common patterns
used throughout the San José pipeline.
"""

from pathlib import Path
from typing import List, Optional, Dict, Any, Union
import json

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq


def read_parquet_safe(
    path: Union[str, Path],
    columns: Optional[List[str]] = None,
    filters: Optional[List] = None,
) -> pd.DataFrame:
    """Read a parquet file with error handling.

    Args:
        path: Path to parquet file
        columns: Optional list of columns to read
        filters: Optional PyArrow filters for predicate pushdown

    Returns:
        DataFrame (empty if file doesn't exist or read fails)
    """
    path = Path(path)
    if not path.exists():
        return pd.DataFrame()

    try:
        return pd.read_parquet(path, columns=columns, filters=filters)
    except Exception as e:
        print(f"Warning: Failed to read {path}: {e}")
        return pd.DataFrame()


def write_parquet_with_metadata(
    df: pd.DataFrame,
    path: Union[str, Path],
    metadata: Optional[Dict[str, Any]] = None,
    compression: str = "snappy",
) -> None:
    """Write DataFrame to parquet with custom metadata.

    Args:
        df: DataFrame to write
        path: Output path
        metadata: Optional dict of metadata to store in parquet file
        compression: Compression codec (snappy, gzip, zstd)
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    # Convert to PyArrow Table
    table = pa.Table.from_pandas(df)

    # Add custom metadata
    if metadata:
        # Serialize metadata to JSON string
        meta_json = json.dumps(metadata).encode('utf-8')
        existing_meta = table.schema.metadata or {}
        existing_meta[b'san_jose_metadata'] = meta_json
        table = table.replace_schema_metadata(existing_meta)

    # Write with specified compression
    pq.write_table(table, path, compression=compression)


def read_parquet_metadata(path: Union[str, Path]) -> Optional[Dict[str, Any]]:
    """Read custom metadata from a parquet file.

    Args:
        path: Path to parquet file

    Returns:
        Metadata dict or None if not present
    """
    path = Path(path)
    if not path.exists():
        return None

    try:
        parquet_file = pq.ParquetFile(path)
        schema = parquet_file.schema_arrow
        if schema.metadata and b'san_jose_metadata' in schema.metadata:
            meta_json = schema.metadata[b'san_jose_metadata']
            return json.loads(meta_json.decode('utf-8'))
    except Exception:
        pass

    return None


def get_parquet_info(path: Union[str, Path]) -> Dict[str, Any]:
    """Get basic info about a parquet file without reading data.

    Args:
        path: Path to parquet file

    Returns:
        Dict with num_rows, num_columns, columns, size_mb
    """
    path = Path(path)
    if not path.exists():
        return {"exists": False}

    try:
        parquet_file = pq.ParquetFile(path)
        metadata = parquet_file.metadata

        return {
            "exists": True,
            "num_rows": metadata.num_rows,
            "num_columns": metadata.num_columns,
            "columns": [f.name for f in parquet_file.schema_arrow],
            "size_mb": round(path.stat().st_size / (1024 * 1024), 2),
            "compression": metadata.row_group(0).column(0).compression if metadata.num_row_groups > 0 else None,
        }
    except Exception as e:
        return {"exists": True, "error": str(e)}


def concat_parquet_files(
    paths: List[Union[str, Path]],
    output_path: Union[str, Path],
    columns: Optional[List[str]] = None,
) -> int:
    """Concatenate multiple parquet files into one.

    Args:
        paths: List of input parquet file paths
        output_path: Output parquet file path
        columns: Optional columns to include

    Returns:
        Total number of rows in output
    """
    dfs = []
    for p in paths:
        df = read_parquet_safe(p, columns=columns)
        if not df.empty:
            dfs.append(df)

    if not dfs:
        return 0

    combined = pd.concat(dfs, ignore_index=True)
    combined.to_parquet(output_path, index=False)
    return len(combined)


def filter_parquet_by_column(
    input_path: Union[str, Path],
    output_path: Union[str, Path],
    column: str,
    values: List[Any],
) -> int:
    """Filter a parquet file by column values.

    Args:
        input_path: Input parquet file
        output_path: Output parquet file
        column: Column name to filter on
        values: List of values to keep

    Returns:
        Number of rows in filtered output
    """
    # Use PyArrow filtering for efficiency
    filters = [(column, 'in', values)]
    df = pd.read_parquet(input_path, filters=filters)

    if df.empty:
        return 0

    df.to_parquet(output_path, index=False)
    return len(df)


class ParquetBatchWriter:
    """Write DataFrames to parquet in batches for memory efficiency."""

    def __init__(
        self,
        output_path: Union[str, Path],
        schema: Optional[pa.Schema] = None,
        compression: str = "snappy",
    ):
        """Initialize batch writer.

        Args:
            output_path: Output parquet file path
            schema: Optional PyArrow schema (inferred from first batch if not provided)
            compression: Compression codec
        """
        self.output_path = Path(output_path)
        self.output_path.parent.mkdir(parents=True, exist_ok=True)
        self.schema = schema
        self.compression = compression
        self.writer = None
        self.n_rows_written = 0

    def write_batch(self, df: pd.DataFrame) -> None:
        """Write a batch of data.

        Args:
            df: DataFrame batch to write
        """
        if df.empty:
            return

        table = pa.Table.from_pandas(df, schema=self.schema)

        if self.writer is None:
            # Initialize writer with schema from first batch
            self.schema = table.schema
            self.writer = pq.ParquetWriter(
                self.output_path,
                self.schema,
                compression=self.compression,
            )

        self.writer.write_table(table)
        self.n_rows_written += len(df)

    def close(self) -> int:
        """Close the writer and return total rows written."""
        if self.writer is not None:
            self.writer.close()
        return self.n_rows_written

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
        return False


if __name__ == "__main__":
    import tempfile

    print("Testing parquet utilities:\n")

    # Create test data
    df = pd.DataFrame({
        "precursor_id": [1, 2, 3, 4, 5],
        "sequence": ["PEPTIDE", "ANOTHER", "TESTING", "SAMPLE", "DATA"],
        "charge": [2, 2, 3, 2, 3],
        "mz": [500.25, 600.33, 450.12, 550.88, 480.42],
    })

    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)

        # Test write with metadata
        test_path = tmpdir / "test.parquet"
        write_parquet_with_metadata(
            df, test_path,
            metadata={"accession": "PXD019086", "step": "test", "n_precursors": 5}
        )
        print(f"Wrote parquet with metadata to {test_path}")

        # Test read metadata
        meta = read_parquet_metadata(test_path)
        print(f"Read metadata: {meta}")

        # Test get info
        info = get_parquet_info(test_path)
        print(f"File info: rows={info['num_rows']}, columns={info['columns']}")

        # Test batch writer
        batch_path = tmpdir / "batch_test.parquet"
        with ParquetBatchWriter(batch_path) as writer:
            writer.write_batch(df.iloc[:2])
            writer.write_batch(df.iloc[2:])
        print(f"Batch writer wrote {writer.n_rows_written} rows")

        # Verify batch output
        df_read = pd.read_parquet(batch_path)
        print(f"Read back {len(df_read)} rows from batch output")

    print("\nAll tests passed!")
