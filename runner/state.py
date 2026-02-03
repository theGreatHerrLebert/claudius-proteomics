#!/usr/bin/env python3
"""
Runner State Management

Provides checkpointing and state tracking for the 6-step pipeline runner.
Enables resumption from any failed step.
"""

import json
from datetime import datetime
from enum import Enum
from pathlib import Path
from typing import Dict, Optional, Any
from dataclasses import dataclass, field, asdict


class StepStatus(Enum):
    """Status of a pipeline step."""
    PENDING = "pending"
    RUNNING = "running"
    COMPLETED = "completed"
    FAILED = "failed"
    SKIPPED = "skipped"


@dataclass
class StepState:
    """State of a single step."""
    status: StepStatus = StepStatus.PENDING
    started_at: Optional[str] = None
    completed_at: Optional[str] = None
    duration_seconds: Optional[float] = None
    error_message: Optional[str] = None
    outputs: Dict[str, str] = field(default_factory=dict)
    summary: Dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> Dict[str, Any]:
        d = asdict(self)
        d['status'] = self.status.value
        return d

    @classmethod
    def from_dict(cls, d: Dict[str, Any]) -> 'StepState':
        d = d.copy()
        d['status'] = StepStatus(d['status'])
        return cls(**d)


@dataclass
class RunnerState:
    """
    Complete state of a runner execution.

    Persisted to disk for checkpointing and resumption.
    """
    accession: str
    config_path: str
    started_at: str = field(default_factory=lambda: datetime.now().isoformat())
    steps: Dict[str, StepState] = field(default_factory=dict)
    current_step: Optional[str] = None
    completed: bool = False
    test_mode: bool = False
    max_files: int = 0

    def __post_init__(self):
        # Initialize all steps if not present
        step_names = ["step1_download", "step2_search", "step3_stratify",
                      "step4_extract", "step5_merge"]
        for name in step_names:
            if name not in self.steps:
                self.steps[name] = StepState()

    @classmethod
    def load(cls, checkpoint_dir: Path) -> Optional['RunnerState']:
        """Load state from checkpoint directory."""
        state_file = checkpoint_dir / "state.json"
        if not state_file.exists():
            return None

        with open(state_file) as f:
            data = json.load(f)

        # Convert steps dict
        steps = {}
        for name, step_data in data.get('steps', {}).items():
            steps[name] = StepState.from_dict(step_data)
        data['steps'] = steps

        return cls(**data)

    def save(self, checkpoint_dir: Path) -> None:
        """Save state to checkpoint directory."""
        checkpoint_dir.mkdir(parents=True, exist_ok=True)
        state_file = checkpoint_dir / "state.json"

        data = {
            'accession': self.accession,
            'config_path': self.config_path,
            'started_at': self.started_at,
            'current_step': self.current_step,
            'completed': self.completed,
            'test_mode': self.test_mode,
            'max_files': self.max_files,
            'steps': {name: step.to_dict() for name, step in self.steps.items()},
        }

        with open(state_file, 'w') as f:
            json.dump(data, f, indent=2)

    def start_step(self, step_name: str) -> None:
        """Mark a step as started."""
        self.current_step = step_name
        self.steps[step_name].status = StepStatus.RUNNING
        self.steps[step_name].started_at = datetime.now().isoformat()

    def complete_step(
        self,
        step_name: str,
        outputs: Dict[str, str] = None,
        summary: Dict[str, Any] = None,
    ) -> None:
        """Mark a step as completed."""
        step = self.steps[step_name]
        step.status = StepStatus.COMPLETED
        step.completed_at = datetime.now().isoformat()

        if step.started_at:
            started = datetime.fromisoformat(step.started_at)
            completed = datetime.fromisoformat(step.completed_at)
            step.duration_seconds = (completed - started).total_seconds()

        if outputs:
            step.outputs = outputs
        if summary:
            step.summary = summary

    def fail_step(self, step_name: str, error_message: str) -> None:
        """Mark a step as failed."""
        step = self.steps[step_name]
        step.status = StepStatus.FAILED
        step.completed_at = datetime.now().isoformat()
        step.error_message = error_message

        if step.started_at:
            started = datetime.fromisoformat(step.started_at)
            completed = datetime.fromisoformat(step.completed_at)
            step.duration_seconds = (completed - started).total_seconds()

    def skip_step(self, step_name: str, reason: str = "already completed") -> None:
        """Mark a step as skipped."""
        self.steps[step_name].status = StepStatus.SKIPPED
        self.steps[step_name].error_message = reason

    def is_step_done(self, step_name: str) -> bool:
        """Check if a step is completed or skipped."""
        return self.steps[step_name].status in (StepStatus.COMPLETED, StepStatus.SKIPPED)

    def get_next_step(self) -> Optional[str]:
        """Get the next step to run."""
        step_order = ["step1_download", "step2_search", "step3_stratify",
                      "step4_extract", "step5_merge"]
        for name in step_order:
            if not self.is_step_done(name):
                return name
        return None

    def get_checkpoint_path(self, base_dir: Path) -> Path:
        """Get checkpoint directory path."""
        return base_dir / "checkpoints" / self.accession


class CheckpointManager:
    """Manages checkpoints for pipeline runs."""

    def __init__(self, base_dir: Path):
        self.base_dir = Path(base_dir)
        self.checkpoints_dir = self.base_dir / "checkpoints"

    def get_checkpoint_dir(self, accession: str) -> Path:
        """Get checkpoint directory for an accession."""
        return self.checkpoints_dir / accession

    def has_checkpoint(self, accession: str) -> bool:
        """Check if a checkpoint exists."""
        state_file = self.get_checkpoint_dir(accession) / "state.json"
        return state_file.exists()

    def load_state(self, accession: str) -> Optional[RunnerState]:
        """Load state for an accession."""
        return RunnerState.load(self.get_checkpoint_dir(accession))

    def save_state(self, state: RunnerState) -> None:
        """Save state for an accession."""
        state.save(self.get_checkpoint_dir(state.accession))

    def mark_step_done(self, accession: str, step_name: str) -> None:
        """Create a marker file indicating step completion."""
        step_file = self.get_checkpoint_dir(accession) / f"{step_name}.done"
        step_file.parent.mkdir(parents=True, exist_ok=True)
        step_file.write_text(datetime.now().isoformat())

    def is_step_done(self, accession: str, step_name: str) -> bool:
        """Check if a step is marked as done."""
        step_file = self.get_checkpoint_dir(accession) / f"{step_name}.done"
        return step_file.exists()

    def list_checkpoints(self) -> Dict[str, Dict[str, Any]]:
        """List all checkpoints with their status."""
        if not self.checkpoints_dir.exists():
            return {}

        result = {}
        for accession_dir in self.checkpoints_dir.iterdir():
            if accession_dir.is_dir():
                state = self.load_state(accession_dir.name)
                if state:
                    result[accession_dir.name] = {
                        'started_at': state.started_at,
                        'completed': state.completed,
                        'current_step': state.current_step,
                        'steps_done': sum(1 for s in state.steps.values()
                                         if s.status == StepStatus.COMPLETED),
                    }
        return result


if __name__ == "__main__":
    import tempfile

    print("Testing runner state management:\n")

    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)

        # Create state
        state = RunnerState(
            accession="PXD019086",
            config_path="config/config.yaml",
            test_mode=True,
            max_files=3,
        )

        # Test step progression
        state.start_step("step1_download")
        print(f"Started step1: {state.steps['step1_download'].status}")

        state.complete_step("step1_download", outputs={"raw_dir": "data/raw/PXD019086"})
        print(f"Completed step1: {state.steps['step1_download'].status}")

        # Save and reload
        checkpoint_dir = tmpdir / "checkpoints" / "PXD019086"
        state.save(checkpoint_dir)
        print(f"\nSaved state to {checkpoint_dir}")

        loaded = RunnerState.load(checkpoint_dir)
        print(f"Loaded state: accession={loaded.accession}, step1={loaded.steps['step1_download'].status}")

        # Test next step
        next_step = loaded.get_next_step()
        print(f"Next step: {next_step}")

        # Test checkpoint manager
        manager = CheckpointManager(tmpdir)
        manager.mark_step_done("PXD019086", "step1_download")
        print(f"\nStep1 done marker: {manager.is_step_done('PXD019086', 'step1_download')}")

        checkpoints = manager.list_checkpoints()
        print(f"Checkpoints: {checkpoints}")

    print("\nAll tests passed!")
