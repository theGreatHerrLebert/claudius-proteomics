#!/usr/bin/env python3
"""Tests for runner state management - checkpointing, step progression, persistence."""

import json
import tempfile
import unittest
from pathlib import Path

import sys
sys.path.insert(0, str(Path(__file__).parent.parent))

from runner.state import StepStatus, StepState, RunnerState, CheckpointManager


class TestStepStatus(unittest.TestCase):
    """Test StepStatus enum."""

    def test_values(self):
        self.assertEqual(StepStatus.PENDING.value, "pending")
        self.assertEqual(StepStatus.RUNNING.value, "running")
        self.assertEqual(StepStatus.COMPLETED.value, "completed")
        self.assertEqual(StepStatus.FAILED.value, "failed")
        self.assertEqual(StepStatus.SKIPPED.value, "skipped")


class TestStepState(unittest.TestCase):
    """Test StepState dataclass."""

    def test_defaults(self):
        s = StepState()
        self.assertEqual(s.status, StepStatus.PENDING)
        self.assertIsNone(s.started_at)
        self.assertIsNone(s.completed_at)
        self.assertEqual(s.outputs, {})

    def test_to_dict(self):
        s = StepState(status=StepStatus.COMPLETED)
        d = s.to_dict()
        self.assertEqual(d["status"], "completed")

    def test_from_dict(self):
        d = {"status": "completed", "started_at": "2026-01-01T00:00:00",
             "completed_at": "2026-01-01T00:01:00", "duration_seconds": 60.0,
             "error_message": None, "outputs": {}, "summary": {}}
        s = StepState.from_dict(d)
        self.assertEqual(s.status, StepStatus.COMPLETED)
        self.assertEqual(s.duration_seconds, 60.0)

    def test_roundtrip(self):
        s = StepState(status=StepStatus.FAILED, error_message="test error")
        d = s.to_dict()
        s2 = StepState.from_dict(d)
        self.assertEqual(s2.status, StepStatus.FAILED)
        self.assertEqual(s2.error_message, "test error")


class TestRunnerState(unittest.TestCase):
    """Test RunnerState - step progression and persistence."""

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()

    def tearDown(self):
        import shutil
        shutil.rmtree(self.tmpdir)

    def test_init_creates_all_steps(self):
        state = RunnerState(accession="PXD123", config_path="config.yaml")
        expected_steps = ["step1_download", "step2_search", "step3_stratify",
                         "step4_extract", "step5_merge", "step6_package"]
        for step in expected_steps:
            self.assertIn(step, state.steps)
            self.assertEqual(state.steps[step].status, StepStatus.PENDING)

    def test_start_step(self):
        state = RunnerState(accession="PXD123", config_path="config.yaml")
        state.start_step("step1_download")
        self.assertEqual(state.steps["step1_download"].status, StepStatus.RUNNING)
        self.assertIsNotNone(state.steps["step1_download"].started_at)
        self.assertEqual(state.current_step, "step1_download")

    def test_complete_step(self):
        state = RunnerState(accession="PXD123", config_path="config.yaml")
        state.start_step("step1_download")
        state.complete_step("step1_download", outputs={"dir": "/tmp/test"})
        step = state.steps["step1_download"]
        self.assertEqual(step.status, StepStatus.COMPLETED)
        self.assertIsNotNone(step.completed_at)
        self.assertIsNotNone(step.duration_seconds)
        self.assertEqual(step.outputs["dir"], "/tmp/test")

    def test_fail_step(self):
        state = RunnerState(accession="PXD123", config_path="config.yaml")
        state.start_step("step2_search")
        state.fail_step("step2_search", "Engine crashed")
        step = state.steps["step2_search"]
        self.assertEqual(step.status, StepStatus.FAILED)
        self.assertEqual(step.error_message, "Engine crashed")

    def test_skip_step(self):
        state = RunnerState(accession="PXD123", config_path="config.yaml")
        state.skip_step("step6_package", "not requested")
        self.assertEqual(state.steps["step6_package"].status, StepStatus.SKIPPED)

    def test_is_step_done(self):
        state = RunnerState(accession="PXD123", config_path="config.yaml")
        self.assertFalse(state.is_step_done("step1_download"))
        state.start_step("step1_download")
        self.assertFalse(state.is_step_done("step1_download"))  # Running != done
        state.complete_step("step1_download")
        self.assertTrue(state.is_step_done("step1_download"))

    def test_is_step_done_skipped(self):
        state = RunnerState(accession="PXD123", config_path="config.yaml")
        state.skip_step("step1_download")
        self.assertTrue(state.is_step_done("step1_download"))

    def test_get_next_step(self):
        state = RunnerState(accession="PXD123", config_path="config.yaml")
        self.assertEqual(state.get_next_step(), "step1_download")
        state.start_step("step1_download")
        state.complete_step("step1_download")
        self.assertEqual(state.get_next_step(), "step2_search")

    def test_get_next_step_all_done(self):
        state = RunnerState(accession="PXD123", config_path="config.yaml")
        for step in state.steps:
            state.start_step(step)
            state.complete_step(step)
        self.assertIsNone(state.get_next_step())

    def test_save_and_load(self):
        state = RunnerState(accession="PXD123", config_path="config.yaml",
                           test_mode=True, max_files=3)
        state.start_step("step1_download")
        state.complete_step("step1_download", outputs={"dir": "/data"})

        checkpoint = Path(self.tmpdir) / "checkpoint"
        state.save(checkpoint)

        loaded = RunnerState.load(checkpoint)
        self.assertIsNotNone(loaded)
        self.assertEqual(loaded.accession, "PXD123")
        self.assertTrue(loaded.test_mode)
        self.assertEqual(loaded.max_files, 3)
        self.assertEqual(loaded.steps["step1_download"].status, StepStatus.COMPLETED)
        self.assertEqual(loaded.steps["step2_search"].status, StepStatus.PENDING)

    def test_load_nonexistent(self):
        result = RunnerState.load(Path(self.tmpdir) / "missing")
        self.assertIsNone(result)


class TestCheckpointManager(unittest.TestCase):
    """Test CheckpointManager."""

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        self.mgr = CheckpointManager(Path(self.tmpdir))

    def tearDown(self):
        import shutil
        shutil.rmtree(self.tmpdir)

    def test_no_checkpoint_initially(self):
        self.assertFalse(self.mgr.has_checkpoint("PXD123"))

    def test_save_and_load(self):
        state = RunnerState(accession="PXD123", config_path="config.yaml")
        state.start_step("step1_download")
        state.complete_step("step1_download")
        self.mgr.save_state(state)

        self.assertTrue(self.mgr.has_checkpoint("PXD123"))
        loaded = self.mgr.load_state("PXD123")
        self.assertEqual(loaded.accession, "PXD123")

    def test_mark_step_done(self):
        self.assertFalse(self.mgr.is_step_done("PXD123", "step1_download"))
        self.mgr.mark_step_done("PXD123", "step1_download")
        self.assertTrue(self.mgr.is_step_done("PXD123", "step1_download"))

    def test_list_checkpoints(self):
        state = RunnerState(accession="PXD123", config_path="config.yaml")
        state.start_step("step1_download")
        state.complete_step("step1_download")
        self.mgr.save_state(state)

        checkpoints = self.mgr.list_checkpoints()
        self.assertIn("PXD123", checkpoints)
        self.assertEqual(checkpoints["PXD123"]["steps_done"], 1)

    def test_list_empty(self):
        checkpoints = self.mgr.list_checkpoints()
        self.assertEqual(checkpoints, {})


if __name__ == "__main__":
    unittest.main(verbosity=2)
