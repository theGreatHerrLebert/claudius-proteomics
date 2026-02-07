"""
Sage Engine Job

Self-contained job module for running Sage search.
Can be executed standalone: python -m runner.engines.sage_job --help
"""

import sys
from pathlib import Path
from typing import Any, Dict, List, Optional

sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from runner.engines.base import EngineJob
from scripts.engine_parsers.sage_parser import SageParser


class SageJob(EngineJob):
    """Sage search engine job."""

    @property
    def engine_name(self) -> str:
        return "sage"

    def _get_parser(self):
        return SageParser()

    def _build_command(
        self,
        accession: str,
        config: Dict[str, Any],
        raw_dir: Path,
        processed_dir: Path,
        fasta_path: Path,
        num_threads: int,
        max_files: int,
    ) -> Optional[List[str]]:
        sage_config = config.get("sage", {})
        sage_path = sage_config.get("path")

        output_dir = processed_dir / "sage"
        output_dir.mkdir(parents=True, exist_ok=True)

        # Get input files - Sage reads .d files directly
        d_files = sorted(raw_dir.glob("*.d"))
        if max_files > 0:
            d_files = d_files[:max_files]

        if not d_files:
            return None

        # Build Sage command - requires config JSON as first positional arg
        sage_config_json = Path(__file__).parent.parent.parent / "config" / "sage_config.json"
        cmd = [
            str(sage_path),
            str(sage_config_json),
            "--fasta", str(fasta_path),
            "--output_directory", str(output_dir),
            "--batch-size", str(max(1, num_threads // 2)),
            "--parquet",
            "--annotate-matches",
        ]

        # Add input files as positional arguments
        cmd.extend([str(f) for f in d_files])

        return cmd

    def _post_subprocess(
        self,
        accession: str,
        config: Dict[str, Any],
        raw_dir: Path,
        processed_dir: Path,
        result,
    ) -> Dict[str, Any]:
        """Count raw PSMs from Sage output."""
        import pandas as pd

        output_dir = processed_dir / "sage"
        results_file = output_dir / "results.sage.parquet"
        n_psms = 0
        if results_file.exists():
            df = pd.read_parquet(results_file)
            n_psms = len(df[~df.get("is_decoy", False)])

        return {
            "n_psms_raw": n_psms,
            "output_dir": str(output_dir),
        }


if __name__ == "__main__":
    SageJob.cli()
