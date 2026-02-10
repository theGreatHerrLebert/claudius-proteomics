"""
Sage Engine Job

Self-contained job module for running Sage search.
Can be executed standalone: python -m runner.engines.sage_job --help
"""

import json
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
        d_files: Optional[List[Path]] = None,
        enzyme_config: Optional[Dict[str, Any]] = None,
    ) -> Optional[List[str]]:
        sage_config = config.get("sage", {})
        sage_path = sage_config.get("path")

        output_dir = processed_dir / "sage"
        output_dir.mkdir(parents=True, exist_ok=True)

        # Get input files: prefer explicit list, fall back to glob
        if d_files is None:
            d_files = sorted(raw_dir.glob("*.d"))
            if max_files > 0:
                d_files = d_files[:max_files]

        if not d_files:
            return None

        # Generate per-group sage config JSON with enzyme settings
        sage_config_json = self._write_sage_config(output_dir, enzyme_config)

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

    @staticmethod
    def _write_sage_config(
        output_dir: Path,
        enzyme_config: Optional[Dict[str, Any]] = None,
    ) -> Path:
        """Write a Sage config JSON, applying enzyme overrides if provided."""
        base_config_path = Path(__file__).parent.parent.parent / "config" / "sage_config.json"
        with open(base_config_path) as f:
            sage_cfg = json.load(f)

        if enzyme_config:
            enzyme = sage_cfg.get("database", {}).get("enzyme", {})
            if "sage_cleave_at" in enzyme_config:
                enzyme["cleave_at"] = enzyme_config["sage_cleave_at"]
            if "sage_restrict" in enzyme_config:
                enzyme["restrict"] = enzyme_config["sage_restrict"]
            if enzyme_config.get("restrict") is None and "sage_restrict" not in enzyme_config:
                enzyme["restrict"] = None
            if "missed_cleavages" in enzyme_config:
                enzyme["missed_cleavages"] = enzyme_config["missed_cleavages"]
            sage_cfg.setdefault("database", {})["enzyme"] = enzyme

        config_path = output_dir / "sage_config.json"
        with open(config_path, "w") as f:
            json.dump(sage_cfg, f, indent=2)
        return config_path

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
