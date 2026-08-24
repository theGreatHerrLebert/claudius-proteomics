"""
FragPipe Engine Job

Self-contained job module for running FragPipe search.
Can be executed standalone: python -m runner.engines.fragpipe_job --help
"""

import json
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional

sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from runner.engines.base import EngineJob
from scripts.engine_parsers.fragpipe_parser import FragPipeParser


class FragPipeJob(EngineJob):
    """FragPipe search engine job."""

    @property
    def engine_name(self) -> str:
        return "fragpipe"

    def _get_parser(self):
        return FragPipeParser()

    def _skip_reason(self, config: Dict[str, Any]) -> Optional[str]:
        fragpipe_config = config.get("fragpipe", {})
        fragpipe_path = fragpipe_config.get("path")
        if not fragpipe_path or not Path(fragpipe_path).exists():
            return "FragPipe not configured"

        runner_script = Path(__file__).parent.parent.parent / "scripts" / "run_fragpipe.py"
        if not runner_script.exists():
            return "run_fragpipe.py not found"

        return None

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
        fragpipe_config = config.get("fragpipe", {})
        fragpipe_path = fragpipe_config.get("path")

        output_dir = processed_dir / "fragpipe_output"
        output_dir.mkdir(parents=True, exist_ok=True)

        runner_script = Path(__file__).parent.parent.parent / "scripts" / "run_fragpipe.py"

        # Cap threads at the configured value (default 64). MSFragger loads
        # Bruker .d files in parallel and the timsdata native library is
        # contention-prone; too many threads triggers JNA native-load failures
        # ("Timsdata could not open the file") partway through a dataset.
        fp_threads = fragpipe_config.get("threads") or num_threads

        cmd = [
            sys.executable, str(runner_script),
            "--fragpipe", str(fragpipe_path),
            "--input", str(raw_dir),
            "--output", str(output_dir),
            "--fasta", str(fasta_path),
            "--threads", str(fp_threads),
            # Stable per-dataset temp dir so JNA extracts/loads the timsdata
            # native lib from one place rather than racing in the node /tmp.
            "--temp-dir", str(processed_dir / "fragpipe_tmp"),
        ]

        # Pass explicit file list if provided
        if d_files:
            cmd.extend(["--files"] + [str(f) for f in d_files])
        elif max_files > 0:
            cmd.extend(["--max-files", str(max_files)])

        # Pass enzyme override if provided
        if enzyme_config and "fragpipe_name" in enzyme_config:
            cmd.extend(["--enzyme", enzyme_config["fragpipe_name"]])

        # Pass the FragPipe base workflow if the mod profile specifies one
        # (e.g. HLA -> "Nonspecific-HLA"). The base workflow is what actually
        # determines MSFragger's specificity (num_enzyme_termini,
        # search_enzyme_cut_1); the --enzyme flag above only renames the
        # enzyme and cannot make a tryptic workflow nonspecific.
        # (run_fragpipe.py resolves the name to workflows/<name>.workflow and
        # errors clearly if the file isn't shipped with FragPipe.)
        #
        # Without this, a run through THIS path falls back to run_fragpipe.py's
        # default LFQ-MBR base workflow. Note that HLA datasets processed before
        # 2026-08-24 are NOT uniformly affected: rendered fragpipe.workflow
        # artifacts show PXD046535/PXD042316 ran nonspecific (termini=0, 7-25)
        # while PXD046755 ran generic tryptic (termini=2, 7-50) — i.e. some runs
        # reached the HLA workflow by another invocation path. Check the rendered
        # workflow per dataset rather than assuming. See TO_ECCB_POSTER.md §4.1.
        if mod_config and mod_config.get("fragpipe_workflow"):
            cmd.extend(["--workflow", mod_config["fragpipe_workflow"]])

        # Pass modification profile as JSON if provided
        if mod_config:
            cmd.extend(["--mod-config-json", json.dumps(mod_config)])

        # Per-dataset MSBooster skip (raw-filename-space datasets, mzML NPE).
        if mod_config and mod_config.get("skip_msbooster"):
            cmd.append("--skip-msbooster")

        return cmd

    def _post_subprocess(
        self,
        accession: str,
        config: Dict[str, Any],
        raw_dir: Path,
        processed_dir: Path,
        result,
    ) -> Dict[str, Any]:
        """Symlink combined_ion.tsv to standard location."""
        output_dir = processed_dir / "fragpipe_output"
        combined_ion = output_dir / "combined_ion.tsv"
        if combined_ion.exists():
            target = processed_dir / "combined_ion.tsv"
            if not target.exists():
                target.symlink_to(combined_ion)

        # Count PSMs from psm.tsv files
        psm_files = list(output_dir.rglob("psm.tsv"))
        n_psms = 0
        for psm_file in psm_files:
            with open(psm_file) as f:
                n_psms += sum(1 for _ in f) - 1  # Subtract header

        return {
            "n_psms_raw": n_psms,
            "n_psm_files": len(psm_files),
            "output_dir": str(output_dir),
        }


if __name__ == "__main__":
    FragPipeJob.cli()
