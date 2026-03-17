#!/usr/bin/env python3
"""Tests for runner summary - StepSummary, write/create functions, manifest."""

import json
import tempfile
import time
import unittest
from pathlib import Path

import sys
sys.path.insert(0, str(Path(__file__).parent.parent))

from runner.summary import (
    StepSummary,
    write_step_summary,
    create_step1_summary,
    create_step2_summary,
    create_step3_summary,
    create_step4_summary,
    create_step5_summary,
    create_manifest,
)


class TestStepSummary(unittest.TestCase):
    """Test StepSummary dataclass."""

    def test_defaults(self):
        s = StepSummary(step_name="step1", accession="PXD123")
        self.assertEqual(s.status, "running")
        self.assertIsNone(s.completed_at)
        self.assertIsNone(s.error_message)

    def test_complete_success(self):
        s = StepSummary(step_name="step1", accession="PXD123")
        s.complete(success=True)
        self.assertEqual(s.status, "success")
        self.assertIsNotNone(s.completed_at)
        self.assertIsNotNone(s.duration_seconds)
        self.assertGreaterEqual(s.duration_seconds, 0)

    def test_complete_failure(self):
        s = StepSummary(step_name="step1", accession="PXD123")
        s.complete(success=False, error_message="something broke")
        self.assertEqual(s.status, "error")
        self.assertEqual(s.error_message, "something broke")

    def test_to_dict(self):
        s = StepSummary(step_name="step1", accession="PXD123",
                       data={"key": "value"}, outputs=["file1.parquet"])
        d = s.to_dict()
        self.assertEqual(d["step_name"], "step1")
        self.assertEqual(d["data"]["key"], "value")
        self.assertEqual(d["outputs"], ["file1.parquet"])

    def test_duration_calculated(self):
        s = StepSummary(step_name="step1", accession="PXD123")
        time.sleep(0.01)
        s.complete(success=True)
        self.assertGreater(s.duration_seconds, 0)


class TestWriteStepSummary(unittest.TestCase):
    """Test writing summaries to disk."""

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()

    def tearDown(self):
        import shutil
        shutil.rmtree(self.tmpdir)

    def test_writes_json(self):
        s = StepSummary(step_name="step1", accession="PXD123")
        s.complete(success=True)
        path = write_step_summary(s, Path(self.tmpdir))
        self.assertTrue(path.exists())
        self.assertEqual(path.name, "step1_summary.json")

    def test_json_valid(self):
        s = StepSummary(step_name="step2", accession="PXD123",
                       data={"n_psms": 50000})
        s.complete(success=True)
        path = write_step_summary(s, Path(self.tmpdir))
        with open(path) as f:
            data = json.load(f)
        self.assertEqual(data["step_name"], "step2")
        self.assertEqual(data["data"]["n_psms"], 50000)

    def test_creates_directory(self):
        subdir = Path(self.tmpdir) / "sub" / "dir"
        s = StepSummary(step_name="step1", accession="PXD123")
        s.complete(success=True)
        path = write_step_summary(s, subdir)
        self.assertTrue(path.exists())


class TestCreateSummaryFactories(unittest.TestCase):
    """Test step-specific summary creation functions."""

    def test_step1_summary(self):
        s = create_step1_summary(
            accession="PXD123", n_raw_files=10, total_size_gb=5.2,
            raw_files=["f1.d", "f2.d"], metadata={"organism": "human"},
            output_dir=Path("/tmp/raw"),
        )
        self.assertEqual(s.step_name, "step1")
        self.assertEqual(s.data["n_raw_files"], 10)
        self.assertEqual(s.data["total_size_gb"], 5.2)

    def test_step2_summary(self):
        s = create_step2_summary(
            accession="PXD123",
            fragpipe_stats={"n_psms": 50000},
            diann_stats={"n_precursors": 45000},
            sage_stats={"n_psms": 48000},
            fasta_info={"n_proteins": 20000},
            output_dir=Path("/tmp/processed"),
        )
        self.assertEqual(s.data["fragpipe"]["n_psms"], 50000)

    def test_step3_summary(self):
        s = create_step3_summary(
            accession="PXD123",
            overlap_stats={"n_all_three": 5000, "n_union": 10000},
            stratified_counts={"all_three": 5000},
            match_tiers={"sequence": 3000, "coordinate": 2000},
            output_dir=Path("/tmp/processed"),
        )
        self.assertEqual(s.data["overlap_stats"]["n_all_three"], 5000)

    def test_step4_summary(self):
        s = create_step4_summary(
            accession="PXD123", n_precursors=100000,
            quality_stats={"rt_r2_median": 0.85}, blob_size_gb=2.1,
            output_dir=Path("/tmp/extracted"),
        )
        self.assertEqual(s.data["n_precursors_extracted"], 100000)

    def test_step5_summary(self):
        s = create_step5_summary(
            accession="PXD123", n_total=100000,
            n_per_engine={"fragpipe": 50000, "diann": 40000, "sage": 48000},
            n_unidentified=20000,
            quality_summary={"pct_high_quality": 15.0},
            output_path=Path("/tmp/store.parquet"),
        )
        self.assertEqual(s.data["n_total_precursors"], 100000)
        self.assertEqual(s.data["n_per_engine"]["fragpipe"], 50000)

    def test_step1_truncates_raw_files(self):
        """Raw files list should be truncated to 10 in summary."""
        many_files = [f"file{i}.d" for i in range(50)]
        s = create_step1_summary(
            accession="PXD123", n_raw_files=50, total_size_gb=25.0,
            raw_files=many_files, metadata={}, output_dir=Path("/tmp"),
        )
        self.assertEqual(len(s.data["raw_files"]), 10)
        self.assertEqual(s.data["raw_files_total"], 50)


class TestCreateManifest(unittest.TestCase):
    """Test manifest generation."""

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()

    def tearDown(self):
        import shutil
        shutil.rmtree(self.tmpdir)

    def test_creates_manifest(self):
        s1 = StepSummary(step_name="step1", accession="PXD123")
        s1.complete(success=True)
        s2 = StepSummary(step_name="step2", accession="PXD123")
        s2.complete(success=True)

        path = create_manifest("PXD123", [s1, s2], Path(self.tmpdir))
        self.assertTrue(path.exists())
        self.assertEqual(path.name, "manifest.json")

    def test_manifest_contents(self):
        s1 = StepSummary(step_name="step1", accession="PXD123",
                        outputs=["raw/"])
        s1.complete(success=True)

        path = create_manifest("PXD123", [s1], Path(self.tmpdir))
        with open(path) as f:
            m = json.load(f)
        self.assertEqual(m["accession"], "PXD123")
        self.assertEqual(m["pipeline_version"], "1.0")
        self.assertIn("step1", m["steps"])
        self.assertEqual(m["final_outputs"], ["raw/"])

    def test_total_duration(self):
        s1 = StepSummary(step_name="step1", accession="PXD123")
        s1.duration_seconds = 10.0
        s2 = StepSummary(step_name="step2", accession="PXD123")
        s2.duration_seconds = 20.0

        path = create_manifest("PXD123", [s1, s2], Path(self.tmpdir))
        with open(path) as f:
            m = json.load(f)
        self.assertAlmostEqual(m["total_duration_seconds"], 30.0)

    def test_empty_summaries(self):
        path = create_manifest("PXD123", [], Path(self.tmpdir))
        with open(path) as f:
            m = json.load(f)
        self.assertEqual(m["total_duration_seconds"], 0)
        self.assertEqual(m["final_outputs"], [])


if __name__ == "__main__":
    unittest.main(verbosity=2)
