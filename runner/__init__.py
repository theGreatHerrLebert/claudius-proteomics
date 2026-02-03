"""
San José Runner

Self-contained runner for processing PRIDE datasets through the 6-step pipeline.
Designed for SLURM cluster submission with checkpointing and resumability.

Steps:
1. Download - Fetch raw data from PRIDE
2. Search - Run FragPipe, DIA-NN, Sage
3. Stratify - Merge and stratify search results
4. Extract - Extract raw 4D signal with quality metrics
5. Merge - Final merge of search + raw data
6. Dashboard - Visualization (separate)
"""

from runner.state import RunnerState, StepStatus
from runner.summary import StepSummary, write_step_summary

__all__ = [
    "RunnerState",
    "StepStatus",
    "StepSummary",
    "write_step_summary",
]
