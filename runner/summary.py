#!/usr/bin/env python3
"""
Summary Generation Utilities

Creates step summary JSON files and manifest for the runner.
"""

import json
from datetime import datetime
from pathlib import Path
from typing import Dict, Any, List, Optional
from dataclasses import dataclass, field, asdict


@dataclass
class StepSummary:
    """Summary of a pipeline step execution."""

    # Identification
    step_name: str
    accession: str

    # Timing
    started_at: str = field(default_factory=lambda: datetime.now().isoformat())
    completed_at: Optional[str] = None
    duration_seconds: Optional[float] = None

    # Status
    status: str = "running"  # running, success, error
    error_message: Optional[str] = None

    # Step-specific data
    data: Dict[str, Any] = field(default_factory=dict)

    # Output files
    outputs: List[str] = field(default_factory=list)

    def complete(self, success: bool = True, error_message: str = None) -> None:
        """Mark summary as completed."""
        self.completed_at = datetime.now().isoformat()
        self.status = "success" if success else "error"
        self.error_message = error_message

        if self.started_at:
            started = datetime.fromisoformat(self.started_at)
            completed = datetime.fromisoformat(self.completed_at)
            self.duration_seconds = (completed - started).total_seconds()

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


def write_step_summary(summary: StepSummary, output_dir: Path) -> Path:
    """Write step summary to JSON file.

    Args:
        summary: StepSummary object
        output_dir: Directory to write summary file

    Returns:
        Path to written summary file
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # File name: step1_summary.json, step2_summary.json, etc.
    summary_file = output_dir / f"{summary.step_name}_summary.json"

    with open(summary_file, 'w') as f:
        json.dump(summary.to_dict(), f, indent=2)

    return summary_file


def create_step1_summary(
    accession: str,
    n_raw_files: int,
    total_size_gb: float,
    raw_files: List[str],
    metadata: Dict[str, Any],
    output_dir: Path,
) -> StepSummary:
    """Create summary for Step 1: Download."""
    summary = StepSummary(
        step_name="step1",
        accession=accession,
        data={
            "n_raw_files": n_raw_files,
            "total_size_gb": total_size_gb,
            "raw_files": raw_files[:10],  # First 10 only for summary
            "raw_files_total": len(raw_files),
            "metadata": metadata,
        },
        outputs=[str(output_dir)],
    )
    return summary


def create_step2_summary(
    accession: str,
    fragpipe_stats: Dict[str, Any],
    diann_stats: Dict[str, Any],
    sage_stats: Dict[str, Any],
    fasta_info: Dict[str, Any],
    output_dir: Path,
) -> StepSummary:
    """Create summary for Step 2: Search."""
    summary = StepSummary(
        step_name="step2",
        accession=accession,
        data={
            "fragpipe": fragpipe_stats,
            "diann": diann_stats,
            "sage": sage_stats,
            "fasta": fasta_info,
        },
        outputs=[str(output_dir)],
    )
    return summary


def create_step3_summary(
    accession: str,
    overlap_stats: Dict[str, Any],
    stratified_counts: Dict[str, int],
    match_tiers: Dict[str, int],
    output_dir: Path,
) -> StepSummary:
    """Create summary for Step 3: Stratify/Merge."""
    summary = StepSummary(
        step_name="step3",
        accession=accession,
        data={
            "overlap_stats": overlap_stats,
            "stratified_counts": stratified_counts,
            "match_tiers": match_tiers,
        },
        outputs=[str(output_dir)],
    )
    return summary


def create_step4_summary(
    accession: str,
    n_precursors: int,
    quality_stats: Dict[str, Any],
    blob_size_gb: float,
    output_dir: Path,
) -> StepSummary:
    """Create summary for Step 4: Extract."""
    summary = StepSummary(
        step_name="step4",
        accession=accession,
        data={
            "n_precursors_extracted": n_precursors,
            "quality_stats": quality_stats,
            "blob_size_gb": blob_size_gb,
        },
        outputs=[str(output_dir)],
    )
    return summary


def create_step5_summary(
    accession: str,
    n_total: int,
    n_per_engine: Dict[str, int],
    n_unidentified: int,
    quality_summary: Dict[str, Any],
    output_path: Path,
) -> StepSummary:
    """Create summary for Step 5: Final Merge."""
    summary = StepSummary(
        step_name="step5",
        accession=accession,
        data={
            "n_total_precursors": n_total,
            "n_per_engine": n_per_engine,
            "n_unidentified": n_unidentified,
            "quality_summary": quality_summary,
        },
        outputs=[str(output_path)],
    )
    return summary


def create_manifest(
    accession: str,
    step_summaries: List[StepSummary],
    output_dir: Path,
) -> Path:
    """Create final manifest.json for a completed run.

    Args:
        accession: PRIDE accession
        step_summaries: List of step summaries
        output_dir: Directory to write manifest

    Returns:
        Path to manifest file
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Aggregate data from all steps
    total_duration = sum(
        s.duration_seconds for s in step_summaries
        if s.duration_seconds is not None
    )

    manifest = {
        "accession": accession,
        "pipeline_version": "1.0",
        "generated_at": datetime.now().isoformat(),
        "total_duration_seconds": total_duration,
        "steps": {s.step_name: s.to_dict() for s in step_summaries},
        "final_outputs": step_summaries[-1].outputs if step_summaries else [],
    }

    manifest_file = output_dir / "manifest.json"
    with open(manifest_file, 'w') as f:
        json.dump(manifest, f, indent=2)

    return manifest_file


if __name__ == "__main__":
    import tempfile

    print("Testing summary utilities:\n")

    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)

        # Create step 1 summary
        s1 = create_step1_summary(
            accession="PXD019086",
            n_raw_files=10,
            total_size_gb=5.2,
            raw_files=["file1.d", "file2.d", "file3.d"],
            metadata={"organism": "Homo sapiens", "instrument": "timsTOF Pro"},
            output_dir=tmpdir / "raw",
        )
        s1.complete(success=True)
        write_step_summary(s1, tmpdir)
        print(f"Step 1 summary: {s1.status}, duration={s1.duration_seconds}")

        # Create step 2 summary
        s2 = create_step2_summary(
            accession="PXD019086",
            fragpipe_stats={"n_psms": 50000, "n_peptides": 10000},
            diann_stats={"n_precursors": 45000},
            sage_stats={"n_psms": 48000},
            fasta_info={"n_proteins": 20000, "path": "db.fasta"},
            output_dir=tmpdir / "processed",
        )
        s2.complete(success=True)
        write_step_summary(s2, tmpdir)
        print(f"Step 2 summary: {s2.status}")

        # Create manifest
        manifest_path = create_manifest("PXD019086", [s1, s2], tmpdir / "merged")
        print(f"\nManifest written to: {manifest_path}")

        # Read back manifest
        with open(manifest_path) as f:
            manifest = json.load(f)
        print(f"Manifest total duration: {manifest['total_duration_seconds']:.1f}s")

    print("\nAll tests passed!")
