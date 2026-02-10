"""
Engine Job Base Module

Defines EngineResult dataclass and EngineJob abstract base class.
Each search engine (FragPipe, DIA-NN, Sage) extends EngineJob to become
a self-contained, independently runnable job that produces canonical parquet output.
"""

import argparse
import json
import subprocess
import sys
from abc import ABC, abstractmethod
from dataclasses import dataclass, field, asdict
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional

@dataclass
class EngineResult:
    """Result of an engine job execution."""

    engine_name: str
    status: str  # "success", "error", "skipped"
    canonical_parquet: Optional[Path] = None
    n_psms: int = 0
    duration_seconds: float = 0.0
    error_message: Optional[str] = None
    extra: Dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> Dict[str, Any]:
        d = asdict(self)
        if d["canonical_parquet"] is not None:
            d["canonical_parquet"] = str(d["canonical_parquet"])
        return d


class EngineJob(ABC):
    """Abstract base class for search engine jobs.

    Subclasses implement:
    - engine_name property
    - _build_command() to construct the subprocess command
    - _get_parser() to return the appropriate engine parser

    The concrete run() method is a template that:
    1. Builds and executes the subprocess (with configurable timeout)
    2. Calls the parser to normalize native output into canonical parquet
    3. Writes status JSON sidecar
    4. Returns EngineResult
    """

    @property
    @abstractmethod
    def engine_name(self) -> str:
        """Return the engine name (lowercase): 'fragpipe', 'diann', 'sage'."""
        pass

    @abstractmethod
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
        """Build the subprocess command list.

        Args:
            d_files: Explicit list of .d paths to process. If None, glob raw_dir.
            enzyme_config: Enzyme settings dict from config["enzymes"][name].
                If None, use hardcoded trypsin defaults.

        Returns:
            Command list, or None if the engine should be skipped
            (e.g. binary not configured).
        """
        pass

    @abstractmethod
    def _get_parser(self):
        """Return an instance of the engine's parser (from scripts.engine_parsers)."""
        pass

    def _post_subprocess(
        self,
        accession: str,
        config: Dict[str, Any],
        raw_dir: Path,
        processed_dir: Path,
        result: subprocess.CompletedProcess,
    ) -> Dict[str, Any]:
        """Optional hook for engine-specific post-processing after subprocess completes.

        For example, DIA-NN converts report.tsv to report.parquet here.

        Returns:
            Dict of extra metadata to include in EngineResult.extra
        """
        return {}

    def _skip_reason(
        self,
        config: Dict[str, Any],
    ) -> Optional[str]:
        """Check if the engine should be skipped. Return reason string or None."""
        engine_config = config.get(self.engine_name, {})
        engine_path = engine_config.get("path")
        if not engine_path or not Path(engine_path).exists():
            return f"{self.engine_name} not configured"
        return None

    def run(
        self,
        accession: str,
        config: Dict[str, Any],
        raw_dir: Path,
        processed_dir: Path,
        fasta_path: Path,
        num_threads: int = 16,
        max_files: int = 0,
        d_files: Optional[List[Path]] = None,
        enzyme_config: Optional[Dict[str, Any]] = None,
    ) -> EngineResult:
        """Execute the engine job: subprocess -> parse -> canonical parquet.

        Args:
            d_files: Explicit list of .d paths. If None, engine globs raw_dir.
            enzyme_config: Enzyme settings dict. If None, use hardcoded trypsin.

        This is the template method that orchestrates the full engine run.
        """
        started_at = datetime.now()

        # Check if engine should be skipped
        skip = self._skip_reason(config)
        if skip is not None:
            result = EngineResult(
                engine_name=self.engine_name,
                status="skipped",
                error_message=skip,
            )
            self._write_status(result, processed_dir, started_at)
            return result

        try:
            # Build command
            cmd = self._build_command(
                accession, config, raw_dir, processed_dir, fasta_path,
                num_threads, max_files,
                d_files=d_files, enzyme_config=enzyme_config,
            )

            if cmd is None:
                result = EngineResult(
                    engine_name=self.engine_name,
                    status="skipped",
                    error_message="Command could not be built",
                )
                self._write_status(result, processed_dir, started_at)
                return result

            # Execute subprocess with timeout
            engine_config = config.get(self.engine_name, {})
            timeout_hours = engine_config.get("timeout_hours", 4)
            timeout_sec = int(timeout_hours * 3600)

            proc = subprocess.run(
                cmd, capture_output=True, text=True, timeout=timeout_sec,
            )

            if proc.returncode != 0:
                elapsed = (datetime.now() - started_at).total_seconds()
                result = EngineResult(
                    engine_name=self.engine_name,
                    status="error",
                    duration_seconds=elapsed,
                    error_message=proc.stderr[:500],
                )
                self._write_status(result, processed_dir, started_at)
                return result

            # Post-subprocess hook (e.g. tsv->parquet conversion)
            extra = self._post_subprocess(
                accession, config, raw_dir, processed_dir, proc,
            )

            # Parse native output into canonical DataFrame
            parser = self._get_parser()
            # Parsers expect (base_dir, accession) where base_dir/accession = processed_dir
            base_dir = processed_dir.parent
            canonical_df = parser.parse(base_dir, accession)

            # Write canonical parquet
            canonical_path = processed_dir / f"{self.engine_name}_canonical.parquet"
            if canonical_df is not None and not canonical_df.empty:
                canonical_df.to_parquet(canonical_path, index=False)
                n_psms = len(canonical_df)
            else:
                n_psms = 0

            elapsed = (datetime.now() - started_at).total_seconds()
            result = EngineResult(
                engine_name=self.engine_name,
                status="success",
                canonical_parquet=canonical_path if n_psms > 0 else None,
                n_psms=n_psms,
                duration_seconds=elapsed,
                extra=extra,
            )
            self._write_status(result, processed_dir, started_at)
            return result

        except subprocess.TimeoutExpired:
            elapsed = (datetime.now() - started_at).total_seconds()
            timeout_hours = config.get(self.engine_name, {}).get("timeout_hours", 4)
            result = EngineResult(
                engine_name=self.engine_name,
                status="error",
                duration_seconds=elapsed,
                error_message=f"Timeout after {timeout_hours} hours",
            )
            self._write_status(result, processed_dir, started_at)
            return result
        except Exception as e:
            elapsed = (datetime.now() - started_at).total_seconds()
            result = EngineResult(
                engine_name=self.engine_name,
                status="error",
                duration_seconds=elapsed,
                error_message=str(e),
            )
            self._write_status(result, processed_dir, started_at)
            return result

    def _write_status(
        self,
        result: EngineResult,
        processed_dir: Path,
        started_at: datetime,
    ) -> None:
        """Write {engine}_status.json sidecar file."""
        processed_dir.mkdir(parents=True, exist_ok=True)
        status_path = processed_dir / f"{self.engine_name}_status.json"

        status = {
            "engine_name": result.engine_name,
            "status": result.status,
            "canonical_parquet": str(result.canonical_parquet) if result.canonical_parquet else None,
            "n_psms": result.n_psms,
            "duration_seconds": result.duration_seconds,
            "error_message": result.error_message,
            "started_at": started_at.isoformat(),
            "completed_at": datetime.now().isoformat(),
        }

        with open(status_path, "w") as f:
            json.dump(status, f, indent=2)

    @classmethod
    def cli(cls) -> None:
        """Common CLI entry point for standalone engine execution.

        Usage in each engine module:
            if __name__ == "__main__":
                FragPipeJob.cli()
        """
        import yaml

        parser = argparse.ArgumentParser(
            description=f"Run {cls.__name__} search engine job",
        )
        parser.add_argument("--accession", "-a", required=True, help="PRIDE accession")
        parser.add_argument("--config", "-c", default="config/config.yaml", help="Config file")
        parser.add_argument("--output-dir", "-o", required=True, help="Output base directory")
        parser.add_argument("--raw-dir", type=Path, help="Raw data directory")
        parser.add_argument("--fasta", type=Path, required=True, help="FASTA database path")
        parser.add_argument("--threads", type=int, default=16, help="Number of threads")
        parser.add_argument("--max-files", type=int, default=0, help="Max input files (0=all)")

        args = parser.parse_args()

        with open(args.config) as f:
            config = yaml.safe_load(f)

        output_base_dir = Path(args.output_dir)
        raw_dir = args.raw_dir or (output_base_dir / "raw" / args.accession)
        processed_dir = output_base_dir / "processed" / args.accession

        job = cls()
        result = job.run(
            accession=args.accession,
            config=config,
            raw_dir=raw_dir,
            processed_dir=processed_dir,
            fasta_path=args.fasta,
            num_threads=args.threads,
            max_files=args.max_files,
        )

        print(f"\n{job.engine_name} job completed: {result.status}")
        if result.canonical_parquet:
            print(f"  Canonical parquet: {result.canonical_parquet}")
        print(f"  PSMs: {result.n_psms}")
        print(f"  Duration: {result.duration_seconds:.1f}s")
        if result.error_message:
            print(f"  Error: {result.error_message}")
