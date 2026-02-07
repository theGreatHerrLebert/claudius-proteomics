#!/usr/bin/env python3
"""
Step 2: Execute Third-Party Search Engines

Runs FragPipe, DIA-NN, and Sage on raw data with FDR=1.0 (full results).
Each engine is a self-contained job (runner/engines/) that produces
canonical parquet + status JSON.

Input: data/raw/{accession}/*.d, FASTA database

Outputs:
- data/processed/{accession}/{engine}_canonical.parquet (per engine)
- data/processed/{accession}/{engine}_status.json (per engine)
- step2_summary.json
"""

import sys
from pathlib import Path
from typing import Dict, Any, List, Optional

# Add parent to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from runner.engines import ENGINES
from runner.summary import StepSummary, write_step_summary


def run_step2_search(
    accession: str,
    config: Dict[str, Any],
    output_base_dir: Path,
    raw_dir: Optional[Path] = None,
    fasta_path: Optional[Path] = None,
    engines: Optional[List[str]] = None,
    num_threads: int = 16,
    max_files: int = 0,
) -> StepSummary:
    """
    Execute Step 2: Run search engines.

    Args:
        accession: PRIDE accession
        config: Pipeline configuration dict
        output_base_dir: Base directory for outputs
        raw_dir: Path to raw data directory (default: data/raw/{accession})
        fasta_path: Path to FASTA database (default from config)
        engines: List of engines to run (default: all three)
        num_threads: Number of threads
        max_files: Maximum number of raw files to process (0=all)

    Returns:
        StepSummary with results
    """
    summary = StepSummary(
        step_name="step2",
        accession=accession,
    )

    if engines is None:
        engines = list(ENGINES.keys())

    if raw_dir is None:
        raw_dir = output_base_dir / "raw" / accession

    processed_dir = output_base_dir / "processed" / accession

    try:
        # Get FASTA path
        if fasta_path is None:
            fasta_path = _get_fasta_path(accession, config, output_base_dir)

        # Verify inputs
        d_files = list(raw_dir.glob("*.d"))
        if not d_files:
            raise FileNotFoundError(f"No .d files found in {raw_dir}")

        print(f"  Running search on {len(d_files)} .d files")
        print(f"  FASTA: {fasta_path}")

        results = {}
        for engine_name in engines:
            if engine_name not in ENGINES:
                results[engine_name] = {"status": "skipped", "reason": f"Unknown engine: {engine_name}"}
                continue

            job = ENGINES[engine_name]()
            result = job.run(
                accession=accession,
                config=config,
                raw_dir=raw_dir,
                processed_dir=processed_dir,
                fasta_path=fasta_path,
                num_threads=num_threads,
                max_files=max_files,
            )
            results[engine_name] = result.to_dict()

        # Ensure all default engines have an entry (for summary consistency)
        for name in ["fragpipe", "diann", "sage"]:
            if name not in results:
                results[name] = {"status": "skipped"}

        # FASTA info
        fasta_info = {
            "path": str(fasta_path),
            "n_proteins": _count_proteins(fasta_path),
        }

        # Update summary
        summary.data = {
            "fragpipe": results.get("fragpipe", {"status": "skipped"}),
            "diann": results.get("diann", {"status": "skipped"}),
            "sage": results.get("sage", {"status": "skipped"}),
            "fasta": fasta_info,
        }
        summary.outputs = [str(processed_dir)]
        summary.complete(success=True)

    except Exception as e:
        summary.complete(success=False, error_message=str(e))
        raise

    # Write summary file
    write_step_summary(summary, processed_dir)

    return summary


def _get_fasta_path(accession: str, config: Dict[str, Any], output_base_dir: Path) -> Path:
    """Get FASTA database path for accession."""
    # Check for accession-specific FASTA in resources/fasta/search_db
    resources_dir = Path(__file__).parent.parent.parent / "resources" / "fasta" / "search_db"
    accession_fasta = resources_dir / f"{accession}_decoys.fasta"
    if accession_fasta.exists():
        return accession_fasta

    # Also check output dir
    fasta_dir = output_base_dir / "resources" / "fasta"
    accession_fasta_alt = fasta_dir / f"{accession}_decoys.fasta"
    if accession_fasta_alt.exists():
        return accession_fasta_alt

    # Check dataset metadata for organism
    organism = config.get("dataset_metadata", {}).get(accession, {}).get("organism")
    if organism:
        organisms = config.get("organisms", {})
        if organism in organisms:
            org_config = organisms[organism]
            # Handle nested config format: {local_fasta: path, ...}
            if isinstance(org_config, dict):
                local_fasta = org_config.get("local_fasta")
                if local_fasta:
                    org_fasta = Path(local_fasta)
                    if org_fasta.exists():
                        return org_fasta
            # Handle simple string path format
            elif isinstance(org_config, str):
                org_fasta = Path(org_config)
                if org_fasta.exists():
                    return org_fasta

    # Fallback to default FASTA
    default_fasta = config.get("fasta", {}).get("default")
    if default_fasta:
        default_path = Path(default_fasta)
        if default_path.exists():
            return default_path

    raise FileNotFoundError(f"No FASTA database found for {accession}")


def _count_proteins(fasta_path: Path) -> int:
    """Count proteins in FASTA file."""
    count = 0
    with open(fasta_path) as f:
        for line in f:
            if line.startswith(">"):
                count += 1
    return count


if __name__ == "__main__":
    import argparse
    import yaml

    parser = argparse.ArgumentParser(description="Step 2: Run search engines")
    parser.add_argument("--accession", "-a", required=True, help="PRIDE accession")
    parser.add_argument("--config", "-c", default="config/config.yaml", help="Config file")
    parser.add_argument("--output-dir", "-o", default="data", help="Output base directory")
    parser.add_argument("--raw-dir", type=Path, help="Raw data directory")
    parser.add_argument("--fasta", type=Path, help="FASTA database path")
    parser.add_argument("--engines", nargs="+", default=["fragpipe", "diann", "sage"])
    parser.add_argument("--threads", type=int, default=16, help="Number of threads")
    parser.add_argument("--max-files", type=int, default=0, help="Max input files (0=all)")

    args = parser.parse_args()

    # Load config
    with open(args.config) as f:
        config = yaml.safe_load(f)

    # Run step
    summary = run_step2_search(
        accession=args.accession,
        config=config,
        output_base_dir=Path(args.output_dir),
        raw_dir=args.raw_dir,
        fasta_path=args.fasta,
        engines=args.engines,
        num_threads=args.threads,
        max_files=args.max_files,
    )

    print(f"\nStep 2 completed: {summary.status}")
    for engine in ["fragpipe", "diann", "sage"]:
        eng_data = summary.data.get(engine, {})
        print(f"  {engine}: {eng_data.get('status')}")
