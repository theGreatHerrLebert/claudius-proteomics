#!/usr/bin/env python3
"""
San José Runner - Main Entry Point

Processes a PRIDE dataset through the 6-step pipeline with checkpointing.
Designed for SLURM cluster submission.

Steps:
1. Download - Fetch raw data from PRIDE
2. Search - Run FragPipe, DIA-NN, Sage
3. Stratify - Merge and stratify search results
4. Extract - Extract raw 4D signal with quality metrics
5. Merge - Final merge of search + raw data
6. Package - Create distributable archive (optional, requires --package)

Usage:
    # Local execution
    python runner/run_dataset.py PXD019086 --config config/config.yaml

    # Test mode (3 files)
    python runner/run_dataset.py PXD019086 --test-mode --max-files 3

    # Resume from checkpoint
    python runner/run_dataset.py PXD019086 --resume

    # Run specific steps
    python runner/run_dataset.py PXD019086 --steps 1 2 3

    # Create distributable archive
    python runner/run_dataset.py PXD019086 --package --package-version 1.0
"""

import argparse
import sys
import traceback
from pathlib import Path
from typing import Dict, Any, List, Optional

import yaml

# Add parent to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from runner.state import RunnerState, CheckpointManager, StepStatus
from runner.summary import StepSummary
from runner.steps import (
    run_step1_download,
    run_step2_search,
    run_step3_stratify,
    run_step4_extract,
    run_step5_merge,
    run_step6_package,
)


def run_dataset(
    accession: str,
    config: Dict[str, Any],
    output_base_dir: Path,
    test_mode: bool = False,
    max_files: int = 3,
    resume: bool = False,
    steps: Optional[List[int]] = None,
    num_threads: int = 16,
    local_data_path: Optional[Path] = None,
    package: bool = False,
    package_version: str = "1.0",
) -> bool:
    """
    Run the San José pipeline for a single dataset.

    Args:
        accession: PRIDE accession
        config: Pipeline configuration dict
        output_base_dir: Base directory for outputs
        test_mode: If True, process limited files
        max_files: Max files in test mode
        resume: Resume from checkpoint if exists
        steps: List of step numbers to run (default: all)
        num_threads: Number of threads
        local_data_path: Path to local data (skip download)
        package: If True, run step 6 to create distributable archive
        package_version: Version string for archive (default: 1.0)

    Returns:
        True if all steps completed successfully
    """
    print("=" * 70)
    print(f"  San José Runner - {accession}")
    print("=" * 70)

    # Initialize checkpoint manager
    checkpoint_manager = CheckpointManager(output_base_dir)

    # Load or create state
    if resume and checkpoint_manager.has_checkpoint(accession):
        print(f"\nResuming from checkpoint...")
        state = checkpoint_manager.load_state(accession)
        if state is None:
            state = _create_new_state(accession, config, test_mode, max_files)
        else:
            print(f"  Started: {state.started_at}")
            print(f"  Current step: {state.current_step}")
    else:
        state = _create_new_state(accession, config, test_mode, max_files)

    # Determine which steps to run
    all_steps = [1, 2, 3, 4, 5]
    if package:
        all_steps.append(6)
    if steps:
        steps_to_run = [s for s in steps if s in all_steps]
    else:
        steps_to_run = all_steps

    print(f"\nConfiguration:")
    print(f"  Test mode: {test_mode}")
    print(f"  Max files: {max_files if test_mode else 'unlimited'}")
    print(f"  Threads: {num_threads}")
    print(f"  Steps to run: {steps_to_run}")

    # Run steps
    success = True
    step_summaries = []

    step_names = ["download", "search", "stratify", "extract", "merge", "package"]
    for step_num in steps_to_run:
        step_name = f"step{step_num}_{step_names[step_num - 1]}"

        # Skip if already completed
        if state.is_step_done(step_name):
            print(f"\n[Step {step_num}] {step_name} - SKIPPED (already completed)")
            continue

        print(f"\n[Step {step_num}] {step_name.upper()}")
        print("-" * 50)

        state.start_step(step_name)
        checkpoint_manager.save_state(state)

        try:
            summary = _run_step(
                step_num=step_num,
                accession=accession,
                config=config,
                output_base_dir=output_base_dir,
                state=state,
                num_threads=num_threads,
                local_data_path=local_data_path,
                package_version=package_version,
            )
            step_summaries.append(summary)

            state.complete_step(
                step_name,
                outputs=dict(zip(summary.outputs, summary.outputs)),
                summary=summary.data,
            )
            checkpoint_manager.mark_step_done(accession, step_name)

            print(f"\n  ✓ Step {step_num} completed in {summary.duration_seconds:.1f}s")

        except Exception as e:
            state.fail_step(step_name, str(e))
            checkpoint_manager.save_state(state)

            print(f"\n  ✗ Step {step_num} FAILED: {e}")
            traceback.print_exc()
            success = False
            break

        checkpoint_manager.save_state(state)

    # Mark as completed if all steps done
    if success:
        state.completed = True
        checkpoint_manager.save_state(state)
        print("\n" + "=" * 70)
        print(f"  Pipeline completed successfully for {accession}")
        print("=" * 70)
    else:
        print("\n" + "=" * 70)
        print(f"  Pipeline FAILED for {accession}")
        print(f"  Resume with: python runner/run_dataset.py {accession} --resume")
        print("=" * 70)

    return success


def _create_new_state(
    accession: str,
    config: Dict[str, Any],
    test_mode: bool,
    max_files: int,
) -> RunnerState:
    """Create a new runner state."""
    return RunnerState(
        accession=accession,
        config_path=str(config.get("_config_path", "config/config.yaml")),
        test_mode=test_mode,
        max_files=max_files,
    )


def _run_step(
    step_num: int,
    accession: str,
    config: Dict[str, Any],
    output_base_dir: Path,
    state: RunnerState,
    num_threads: int,
    local_data_path: Optional[Path] = None,
    package_version: str = "1.0",
) -> StepSummary:
    """Run a single pipeline step."""

    if step_num == 1:
        return run_step1_download(
            accession=accession,
            config=config,
            output_base_dir=output_base_dir,
            test_mode=state.test_mode,
            max_files=state.max_files if state.test_mode else 0,
            local_data_path=local_data_path,
        )

    elif step_num == 2:
        return run_step2_search(
            accession=accession,
            config=config,
            output_base_dir=output_base_dir,
            num_threads=num_threads,
            max_files=state.max_files if state.test_mode else 0,
        )

    elif step_num == 3:
        return run_step3_stratify(
            accession=accession,
            config=config,
            output_base_dir=output_base_dir,
        )

    elif step_num == 4:
        return run_step4_extract(
            accession=accession,
            config=config,
            output_base_dir=output_base_dir,
            num_threads=num_threads,
            max_files=state.max_files if state.test_mode else 0,
        )

    elif step_num == 5:
        return run_step5_merge(
            accession=accession,
            config=config,
            output_base_dir=output_base_dir,
        )

    elif step_num == 6:
        return run_step6_package(
            accession=accession,
            config=config,
            output_base_dir=output_base_dir,
            version=package_version,
        )

    else:
        raise ValueError(f"Unknown step number: {step_num}")


def main():
    parser = argparse.ArgumentParser(
        description="San José Runner - Process PRIDE datasets through the pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Full run
    python runner/run_dataset.py PXD019086

    # Test mode (3 files)
    python runner/run_dataset.py PXD019086 --test-mode

    # Resume from checkpoint
    python runner/run_dataset.py PXD019086 --resume

    # Run specific steps
    python runner/run_dataset.py PXD019086 --steps 1 2 3

    # Use local data
    python runner/run_dataset.py PXD019086 --local-data /path/to/data
        """,
    )

    parser.add_argument(
        "accession",
        help="PRIDE accession (e.g., PXD019086)",
    )
    parser.add_argument(
        "--config", "-c",
        default="config/config.yaml",
        help="Configuration file (default: config/config.yaml)",
    )
    parser.add_argument(
        "--output-dir", "-o",
        default="data",
        help="Output base directory (default: data)",
    )
    parser.add_argument(
        "--test-mode", "-t",
        action="store_true",
        help="Test mode - process limited files",
    )
    parser.add_argument(
        "--max-files",
        type=int,
        default=3,
        help="Maximum files in test mode (default: 3)",
    )
    parser.add_argument(
        "--resume", "-r",
        action="store_true",
        help="Resume from checkpoint if exists",
    )
    parser.add_argument(
        "--steps",
        type=int,
        nargs="+",
        help="Specific steps to run (1-5)",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=16,
        help="Number of threads (default: 16)",
    )
    parser.add_argument(
        "--local-data",
        type=Path,
        help="Path to local data (skip download)",
    )
    parser.add_argument(
        "--package", "-p",
        action="store_true",
        help="Create distributable archive after processing (step 6)",
    )
    parser.add_argument(
        "--package-version",
        default="1.0",
        help="Version string for archive (default: 1.0)",
    )
    parser.add_argument(
        "--enrich-metadata",
        action="store_true",
        help="Re-run only paper metadata extraction (for manually-placed PDFs)",
    )

    args = parser.parse_args()

    # Load configuration
    config_path = Path(args.config)
    if not config_path.exists():
        print(f"Error: Config file not found: {config_path}")
        sys.exit(1)

    with open(config_path) as f:
        config = yaml.safe_load(f)

    # Store config path for reference
    config["_config_path"] = str(config_path)

    # Handle --enrich-metadata: re-run only paper extraction, then exit
    if args.enrich_metadata:
        from scripts.paper_metadata import run_paper_extraction

        metadata_dir = Path(args.output_dir) / "metadata" / args.accession
        # Force re-run by overriding skip_if_exists
        config.setdefault("paper_extraction", {})["skip_if_exists"] = False

        print(f"Re-running paper metadata extraction for {args.accession}")
        print(f"  Metadata dir: {metadata_dir}")
        result = run_paper_extraction(args.accession, metadata_dir, config)
        if result and result.get("status") == "success":
            n = len(result.get("fields", {}))
            print(f"\nPaper extraction: {n} fields extracted")
            sys.exit(0)
        elif result:
            print(f"\nPaper extraction: {result.get('status', 'unknown')}")
            sys.exit(0)
        else:
            print("\nPaper extraction: skipped")
            sys.exit(0)

    # Run pipeline
    success = run_dataset(
        accession=args.accession,
        config=config,
        output_base_dir=Path(args.output_dir),
        test_mode=args.test_mode,
        max_files=args.max_files,
        resume=args.resume,
        steps=args.steps,
        num_threads=args.threads,
        local_data_path=args.local_data,
        package=args.package,
        package_version=args.package_version,
    )

    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()
