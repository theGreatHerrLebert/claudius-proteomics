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
        d_files: Optional[List[Path]] = None,
        enzyme_config: Optional[Dict[str, Any]] = None,
        mod_config: Optional[Dict[str, Any]] = None,
    ) -> Optional[List[str]]:
        diann_config = config.get("diann", {})
        diann_path = diann_config.get("path")

        output_dir = processed_dir / "diann"
        output_dir.mkdir(parents=True, exist_ok=True)

        # Get .d files: prefer explicit list, fall back to glob
        if d_files is None:
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

        # Digestion settings
        diann_cut = "K*,R*,!*P"  # Default: trypsin
        missed_cleavages = "2"
        if enzyme_config:
            diann_cut = enzyme_config.get("diann_cut", diann_cut)
            missed_cleavages = str(enzyme_config.get("missed_cleavages", 2))

        # Nonspecific enzyme: omit --cut entirely, set --missed-cleavages 0
        if diann_cut == "":
            cmd.extend(["--missed-cleavages", "0"])
        else:
            cmd.extend([
                "--cut", diann_cut,
                "--missed-cleavages", missed_cleavages,
            ])

        # Modification settings
        if mod_config:
            # Dynamic: build from profile
            for fm in mod_config.get("fixed_modifications", []):
                residues = ",".join(fm["residues"])
                cmd.extend(["--fixed-mod", f"UniMod:{fm['unimod_id']},{fm['mass']},{residues}"])
            for vm in mod_config.get("variable_modifications", []):
                if vm.get("site") == "N-term":
                    cmd.extend(["--var-mod", f"UniMod:{vm['unimod_id']},{vm['mass']},*n"])
                else:
                    residues = ",".join(vm["residues"])
                    cmd.extend(["--var-mod", f"UniMod:{vm['unimod_id']},{vm['mass']},{residues}"])
            cmd.extend(["--max-var-mods", str(mod_config.get("max_variable_mods", 3))])
            cmd.extend([
                "--min-pep-len", str(mod_config.get("min_peptide_length", 7)),
                "--max-pep-len", str(mod_config.get("max_peptide_length", 50)),
            ])
        else:
            # Fallback to hardcoded defaults (standard profile)
            cmd.extend([
                "--var-mod", "UniMod:35,15.994915,M",  # Oxidation (M)
                "--var-mod", "UniMod:1,42.010565,*n",  # N-term Acetyl
                "--fixed-mod", "UniMod:4,57.021464,C",  # Carbamidomethyl (C)
                "--max-var-mods", "3",
                "--min-pep-len", "7",
                "--max-pep-len", "50",
            ])

        cmd.extend([
            "--min-pr-charge", "1",
            "--max-pr-charge", "4",
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
