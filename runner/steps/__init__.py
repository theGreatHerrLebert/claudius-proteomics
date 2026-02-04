"""
San José Pipeline Steps

Each step is a self-contained module that can be run independently
with checkpointing support.

Steps:
1. step1_download - Download raw data from PRIDE
2. step2_search - Run FragPipe, DIA-NN, Sage
3. step3_stratify - Merge and stratify search results
4. step4_extract - Extract raw 4D signal with quality metrics
5. step5_merge - Final merge of search + raw data
6. step6_package - Create distributable archive
"""

from runner.steps.step1_download import run_step1_download
from runner.steps.step2_search import run_step2_search
from runner.steps.step3_stratify import run_step3_stratify
from runner.steps.step4_extract import run_step4_extract
from runner.steps.step5_merge import run_step5_merge
from runner.steps.step6_package import run_step6_package

__all__ = [
    "run_step1_download",
    "run_step2_search",
    "run_step3_stratify",
    "run_step4_extract",
    "run_step5_merge",
    "run_step6_package",
]
