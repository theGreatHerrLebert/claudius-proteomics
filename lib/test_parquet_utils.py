#!/usr/bin/env python3
"""Tests for parquet_utils.py - read/write, metadata, batch writer."""

import tempfile
import unittest
from pathlib import Path

import pandas as pd
import numpy as np

import sys
sys.path.insert(0, str(Path(__file__).parent))

from parquet_utils import (
    read_parquet_safe,
    write_parquet_with_metadata,
    read_parquet_metadata,
    get_parquet_info,
    concat_parquet_files,
    filter_parquet_by_column,
    ParquetBatchWriter,
)


class TestReadParquetSafe(unittest.TestCase):
    """Test safe parquet reading."""

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        self.df = pd.DataFrame({
            "id": [1, 2, 3], "value": [10.0, 20.0, 30.0], "label": ["a", "b", "c"]
        })
        self.path = Path(self.tmpdir) / "test.parquet"
        self.df.to_parquet(self.path, index=False)

    def tearDown(self):
        import shutil
        shutil.rmtree(self.tmpdir)

    def test_read_all(self):
        result = read_parquet_safe(self.path)
        self.assertEqual(len(result), 3)
        pd.testing.assert_frame_equal(result, self.df)

    def test_read_columns(self):
        result = read_parquet_safe(self.path, columns=["id", "value"])
        self.assertEqual(list(result.columns), ["id", "value"])

    def test_missing_file_returns_empty(self):
        result = read_parquet_safe(Path(self.tmpdir) / "missing.parquet")
        self.assertTrue(result.empty)

    def test_string_path(self):
        result = read_parquet_safe(str(self.path))
        self.assertEqual(len(result), 3)


class TestWriteParquetWithMetadata(unittest.TestCase):
    """Test parquet writing with metadata."""

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        self.df = pd.DataFrame({"x": [1, 2, 3]})

    def tearDown(self):
        import shutil
        shutil.rmtree(self.tmpdir)

    def test_write_basic(self):
        path = Path(self.tmpdir) / "out.parquet"
        write_parquet_with_metadata(self.df, path)
        self.assertTrue(path.exists())
        result = pd.read_parquet(path)
        self.assertEqual(len(result), 3)

    def test_write_with_metadata(self):
        path = Path(self.tmpdir) / "out.parquet"
        meta = {"accession": "PXD019086", "version": "1.0", "n_rows": 3}
        write_parquet_with_metadata(self.df, path, metadata=meta)
        recovered = read_parquet_metadata(path)
        self.assertEqual(recovered["accession"], "PXD019086")
        self.assertEqual(recovered["version"], "1.0")
        self.assertEqual(recovered["n_rows"], 3)

    def test_creates_parent_dirs(self):
        path = Path(self.tmpdir) / "sub" / "dir" / "out.parquet"
        write_parquet_with_metadata(self.df, path)
        self.assertTrue(path.exists())

    def test_no_metadata(self):
        path = Path(self.tmpdir) / "out.parquet"
        write_parquet_with_metadata(self.df, path, metadata=None)
        meta = read_parquet_metadata(path)
        self.assertIsNone(meta)


class TestReadParquetMetadata(unittest.TestCase):
    """Test metadata reading."""

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()

    def tearDown(self):
        import shutil
        shutil.rmtree(self.tmpdir)

    def test_missing_file_returns_none(self):
        result = read_parquet_metadata(Path(self.tmpdir) / "missing.parquet")
        self.assertIsNone(result)

    def test_file_without_metadata(self):
        path = Path(self.tmpdir) / "plain.parquet"
        pd.DataFrame({"x": [1]}).to_parquet(path)
        result = read_parquet_metadata(path)
        self.assertIsNone(result)

    def test_roundtrip_metadata(self):
        path = Path(self.tmpdir) / "meta.parquet"
        meta = {"key1": "value1", "nested": {"a": 1, "b": [1, 2, 3]}}
        write_parquet_with_metadata(pd.DataFrame({"x": [1]}), path, metadata=meta)
        recovered = read_parquet_metadata(path)
        self.assertEqual(recovered["key1"], "value1")
        self.assertEqual(recovered["nested"]["a"], 1)
        self.assertEqual(recovered["nested"]["b"], [1, 2, 3])


class TestGetParquetInfo(unittest.TestCase):
    """Test parquet info extraction."""

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()

    def tearDown(self):
        import shutil
        shutil.rmtree(self.tmpdir)

    def test_basic_info(self):
        path = Path(self.tmpdir) / "test.parquet"
        df = pd.DataFrame({"a": [1, 2, 3], "b": [4.0, 5.0, 6.0]})
        df.to_parquet(path, index=False)
        info = get_parquet_info(path)
        self.assertTrue(info["exists"])
        self.assertEqual(info["num_rows"], 3)
        self.assertEqual(info["num_columns"], 2)
        self.assertIn("a", info["columns"])
        self.assertIn("b", info["columns"])
        self.assertGreaterEqual(info["size_mb"], 0)

    def test_missing_file(self):
        info = get_parquet_info(Path(self.tmpdir) / "missing.parquet")
        self.assertFalse(info["exists"])


class TestConcatParquetFiles(unittest.TestCase):
    """Test parquet concatenation."""

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()

    def tearDown(self):
        import shutil
        shutil.rmtree(self.tmpdir)

    def test_concat_two_files(self):
        p1 = Path(self.tmpdir) / "a.parquet"
        p2 = Path(self.tmpdir) / "b.parquet"
        out = Path(self.tmpdir) / "combined.parquet"

        pd.DataFrame({"x": [1, 2]}).to_parquet(p1, index=False)
        pd.DataFrame({"x": [3, 4]}).to_parquet(p2, index=False)

        n = concat_parquet_files([p1, p2], out)
        self.assertEqual(n, 4)
        result = pd.read_parquet(out)
        self.assertEqual(list(result["x"]), [1, 2, 3, 4])

    def test_concat_with_missing_file(self):
        p1 = Path(self.tmpdir) / "a.parquet"
        p_missing = Path(self.tmpdir) / "missing.parquet"
        out = Path(self.tmpdir) / "combined.parquet"

        pd.DataFrame({"x": [1, 2]}).to_parquet(p1, index=False)

        n = concat_parquet_files([p1, p_missing], out)
        self.assertEqual(n, 2)

    def test_concat_empty_list(self):
        out = Path(self.tmpdir) / "empty.parquet"
        n = concat_parquet_files([], out)
        self.assertEqual(n, 0)

    def test_concat_with_column_filter(self):
        p1 = Path(self.tmpdir) / "a.parquet"
        out = Path(self.tmpdir) / "combined.parquet"

        pd.DataFrame({"x": [1], "y": [2]}).to_parquet(p1, index=False)
        n = concat_parquet_files([p1], out, columns=["x"])
        result = pd.read_parquet(out)
        self.assertEqual(list(result.columns), ["x"])


class TestFilterParquetByColumn(unittest.TestCase):
    """Test parquet filtering."""

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        self.path = Path(self.tmpdir) / "data.parquet"
        pd.DataFrame({
            "raw_file": ["run1", "run1", "run2", "run3"],
            "value": [10, 20, 30, 40],
        }).to_parquet(self.path, index=False)

    def tearDown(self):
        import shutil
        shutil.rmtree(self.tmpdir)

    def test_filter_single_value(self):
        out = Path(self.tmpdir) / "filtered.parquet"
        n = filter_parquet_by_column(self.path, out, "raw_file", ["run1"])
        self.assertEqual(n, 2)
        result = pd.read_parquet(out)
        self.assertTrue(all(result["raw_file"] == "run1"))

    def test_filter_multiple_values(self):
        out = Path(self.tmpdir) / "filtered.parquet"
        n = filter_parquet_by_column(self.path, out, "raw_file", ["run1", "run3"])
        self.assertEqual(n, 3)

    def test_filter_no_match(self):
        out = Path(self.tmpdir) / "filtered.parquet"
        n = filter_parquet_by_column(self.path, out, "raw_file", ["nonexistent"])
        self.assertEqual(n, 0)


class TestParquetBatchWriter(unittest.TestCase):
    """Test batch parquet writing."""

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()

    def tearDown(self):
        import shutil
        shutil.rmtree(self.tmpdir)

    def test_write_multiple_batches(self):
        path = Path(self.tmpdir) / "batch.parquet"
        with ParquetBatchWriter(path) as writer:
            writer.write_batch(pd.DataFrame({"x": [1, 2]}))
            writer.write_batch(pd.DataFrame({"x": [3, 4]}))
        self.assertEqual(writer.n_rows_written, 4)
        result = pd.read_parquet(path)
        self.assertEqual(len(result), 4)
        self.assertEqual(list(result["x"]), [1, 2, 3, 4])

    def test_skip_empty_batch(self):
        path = Path(self.tmpdir) / "batch.parquet"
        with ParquetBatchWriter(path) as writer:
            writer.write_batch(pd.DataFrame({"x": [1]}))
            writer.write_batch(pd.DataFrame())  # Empty, should be skipped
        self.assertEqual(writer.n_rows_written, 1)

    def test_close_returns_count(self):
        path = Path(self.tmpdir) / "batch.parquet"
        writer = ParquetBatchWriter(path)
        writer.write_batch(pd.DataFrame({"x": [1, 2, 3]}))
        n = writer.close()
        self.assertEqual(n, 3)

    def test_creates_parent_dirs(self):
        path = Path(self.tmpdir) / "sub" / "dir" / "batch.parquet"
        with ParquetBatchWriter(path) as writer:
            writer.write_batch(pd.DataFrame({"x": [1]}))
        self.assertTrue(path.exists())

    def test_no_batches_written(self):
        path = Path(self.tmpdir) / "empty.parquet"
        with ParquetBatchWriter(path) as writer:
            pass
        self.assertEqual(writer.n_rows_written, 0)


if __name__ == "__main__":
    unittest.main(verbosity=2)
