"""
DIA-NN Engine Job

Self-contained job module for running DIA-NN search.
Can be executed standalone: python -m runner.engines.diann_job --help

Includes spectral library reuse: if report-lib.predicted.speclib exists
from a prior run, DIA-NN reuses it and skips deep learning prediction.
"""

import sys
from pathlib import Path
from typing import Any, Dict, List, Optional

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from runner.engines.base import EngineJob
from scripts.engine_parsers.diann_parser import DiannParser


class DiannJob(EngineJob):
    """DIA-NN search engine job."""

    @property
    def engine_name(self) -> str:
        return "diann"

    def _get_parser(self):
        return DiannParser()

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
        diann_config = config.get("diann", {})
        diann_path = diann_config.get("path")

        output_dir = processed_dir / "diann"
        output_dir.mkdir(parents=True, exist_ok=True)

        # Get .d files
        d_files = sorted(raw_dir.glob("*.d"))
        if max_files > 0:
            d_files = d_files[:max_files]

        if not d_files:
            return None

        # Check for existing predicted speclib to skip library generation
        existing_speclib = output_dir / "report-lib.predicted.speclib"

        cmd = [
            str(diann_path),
            "--out", str(output_dir / "report.tsv"),
            "--qvalue", "1.0",  # No FDR filtering
            "--threads", str(num_threads),
            "--dda",  # DDA mode for DDA datasets
        ]

        if existing_speclib.exists():
            # Reuse existing predicted spectral library, skip prediction
            print(f"  DIA-NN: reusing existing speclib ({existing_speclib})")
            cmd.extend([
                "--lib", str(existing_speclib),
                "--fasta", str(fasta_path),  # For protein annotation
            ])
        else:
            # Full library-free search with prediction
            cmd.extend([
                "--fasta", str(fasta_path),
                "--fasta-search",  # Enable FASTA digest for library-free search
                "--predictor",  # Enable deep learning
            ])

        # Digestion and modification settings (always needed)
        cmd.extend([
            "--cut", "K*,R*,!*P",  # Trypsin: cleave at K/R, not before P
            "--missed-cleavages", "2",
            "--min-pep-len", "7",
            "--max-pep-len", "50",
            "--min-pr-charge", "1",
            "--max-pr-charge", "4",
            "--var-mod", "UniMod:35,15.994915,M",  # Oxidation (M)
            "--var-mod", "UniMod:1,42.010565,*n",  # N-term Acetyl
            "--fixed-mod", "UniMod:4,57.021464,C",  # Carbamidomethyl (C)
            "--max-var-mods", "3",
            "--met-excision",  # N-terminal methionine excision
        ])

        # Add input files
        for d_file in d_files:
            cmd.extend(["--f", str(d_file)])

        return cmd

    def _post_subprocess(
        self,
        accession: str,
        config: Dict[str, Any],
        raw_dir: Path,
        processed_dir: Path,
        result,
    ) -> Dict[str, Any]:
        """Convert report.tsv to report.parquet for the parser."""
        output_dir = processed_dir / "diann"
        report_tsv = output_dir / "report.tsv"
        report_parquet = output_dir / "report.parquet"

        n_precursors = 0
        if report_tsv.exists():
            df = pd.read_csv(report_tsv, sep="\t")
            df.to_parquet(report_parquet, index=False)
            n_precursors = len(df)

        return {
            "n_precursors_raw": n_precursors,
            "output_dir": str(output_dir),
        }


if __name__ == "__main__":
    DiannJob.cli()
