#!/usr/bin/env python3
"""
Step 2: Execute Third-Party Search Engines

Runs FragPipe, DIA-NN, and Sage on raw data with FDR=1.0 (full results).

Input: data/raw/{accession}/*.d, FASTA database

Outputs:
- data/processed/{accession}/fragpipe/combined_ion.tsv
- data/processed/{accession}/diann/report.parquet
- data/processed/{accession}/sage/results.sage.parquet
- step2_summary.json
"""

import subprocess
import sys
from pathlib import Path
from typing import Dict, Any, List, Optional

# Add parent to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

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

    Returns:
        StepSummary with results
    """
    summary = StepSummary(
        step_name="step2",
        accession=accession,
    )

    if engines is None:
        engines = ["fragpipe", "diann", "sage"]

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

        results = {
            "fragpipe": {"status": "skipped"},
            "diann": {"status": "skipped"},
            "sage": {"status": "skipped"},
        }

        # Run FragPipe
        if "fragpipe" in engines:
            fp_result = _run_fragpipe(
                accession, config, raw_dir, processed_dir, fasta_path, num_threads, max_files
            )
            results["fragpipe"] = fp_result

        # Run DIA-NN
        if "diann" in engines:
            dn_result = _run_diann(
                accession, config, raw_dir, processed_dir, fasta_path, num_threads, max_files
            )
            results["diann"] = dn_result

        # Run Sage
        if "sage" in engines:
            sg_result = _run_sage(
                accession, config, raw_dir, processed_dir, fasta_path, num_threads, max_files
            )
            results["sage"] = sg_result

        # FASTA info
        fasta_info = {
            "path": str(fasta_path),
            "n_proteins": _count_proteins(fasta_path),
        }

        # Update summary
        summary.data = {
            "fragpipe": results["fragpipe"],
            "diann": results["diann"],
            "sage": results["sage"],
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


def _run_fragpipe(
    accession: str,
    config: Dict[str, Any],
    raw_dir: Path,
    processed_dir: Path,
    fasta_path: Path,
    num_threads: int,
    max_files: int = 0,
) -> Dict[str, Any]:
    """Run FragPipe search."""
    fragpipe_config = config.get("fragpipe", {})
    fragpipe_path = fragpipe_config.get("path")

    if not fragpipe_path or not Path(fragpipe_path).exists():
        return {"status": "skipped", "reason": "FragPipe not configured"}

    output_dir = processed_dir / "fragpipe_output"
    output_dir.mkdir(parents=True, exist_ok=True)

    # Use Snakemake rule or direct execution
    try:
        # Check if run_fragpipe.py exists
        runner_script = Path(__file__).parent.parent.parent / "scripts" / "run_fragpipe.py"

        if runner_script.exists():
            cmd = [
                sys.executable, str(runner_script),
                "--fragpipe", str(fragpipe_path),
                "--input", str(raw_dir),
                "--output", str(output_dir),
                "--fasta", str(fasta_path),
                "--threads", str(num_threads),
            ]
            if max_files > 0:
                cmd.extend(["--max-files", str(max_files)])
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=3600 * 4)

            if result.returncode != 0:
                return {
                    "status": "error",
                    "error": result.stderr[:500],
                }
        else:
            # Direct FragPipe execution would go here
            return {"status": "skipped", "reason": "run_fragpipe.py not found"}

        # Parse results
        psm_files = list(output_dir.rglob("psm.tsv"))
        n_psms = 0
        for psm_file in psm_files:
            with open(psm_file) as f:
                n_psms += sum(1 for _ in f) - 1  # Subtract header

        # Copy to standard location
        combined_ion = output_dir / "combined_ion.tsv"
        if combined_ion.exists():
            target = processed_dir / "combined_ion.tsv"
            if not target.exists():
                target.symlink_to(combined_ion)

        return {
            "status": "success",
            "n_psms": n_psms,
            "n_files": len(psm_files),
            "output_dir": str(output_dir),
        }

    except subprocess.TimeoutExpired:
        return {"status": "error", "error": "Timeout after 4 hours"}
    except Exception as e:
        return {"status": "error", "error": str(e)}


def _run_diann(
    accession: str,
    config: Dict[str, Any],
    raw_dir: Path,
    processed_dir: Path,
    fasta_path: Path,
    num_threads: int,
    max_files: int = 0,
) -> Dict[str, Any]:
    """Run DIA-NN search."""
    diann_config = config.get("diann", {})
    diann_path = diann_config.get("path")

    if not diann_path or not Path(diann_path).exists():
        return {"status": "skipped", "reason": "DIA-NN not configured"}

    output_dir = processed_dir / "diann"
    output_dir.mkdir(parents=True, exist_ok=True)

    try:
        # Get .d files (DIA-NN reads Bruker .d directly)
        d_files = sorted(raw_dir.glob("*.d"))
        if max_files > 0:
            d_files = d_files[:max_files]

        # Build DIA-NN command
        cmd = [
            str(diann_path),
            "--fasta", str(fasta_path),
            "--fasta-search",  # Enable FASTA digest for library-free search
            "--out", str(output_dir / "report.tsv"),
            "--qvalue", "1.0",  # No FDR filtering
            "--threads", str(num_threads),
            "--predictor",  # Enable deep learning
            "--dda",  # DDA mode for DDA datasets
            # Modifications to match FragPipe/Sage
            "--var-mod", "UniMod:35,15.994915,M",  # Oxidation (M)
            "--fixed-mod", "UniMod:4,57.021464,C",  # Carbamidomethyl (C)
        ]

        # Add input files
        for d_file in d_files:
            cmd.extend(["--f", str(d_file)])

        result = subprocess.run(cmd, capture_output=True, text=True, timeout=3600 * 4)

        if result.returncode != 0:
            return {
                "status": "error",
                "error": result.stderr[:500],
            }

        # Convert to parquet
        report_tsv = output_dir / "report.tsv"
        report_parquet = output_dir / "report.parquet"

        if report_tsv.exists():
            import pandas as pd
            df = pd.read_csv(report_tsv, sep="\t")
            df.to_parquet(report_parquet, index=False)
            n_precursors = len(df)
        else:
            n_precursors = 0

        return {
            "status": "success",
            "n_precursors": n_precursors,
            "output_dir": str(output_dir),
        }

    except subprocess.TimeoutExpired:
        return {"status": "error", "error": "Timeout after 4 hours"}
    except Exception as e:
        return {"status": "error", "error": str(e)}


def _run_sage(
    accession: str,
    config: Dict[str, Any],
    raw_dir: Path,
    processed_dir: Path,
    fasta_path: Path,
    num_threads: int,
    max_files: int = 0,
) -> Dict[str, Any]:
    """Run Sage search."""
    sage_config = config.get("sage", {})
    sage_path = sage_config.get("path")

    if not sage_path or not Path(sage_path).exists():
        return {"status": "skipped", "reason": "Sage not configured"}

    output_dir = processed_dir / "sage"
    output_dir.mkdir(parents=True, exist_ok=True)

    try:
        # Get input files - Sage can read .d files directly, prefer them over mzML
        d_files = sorted(raw_dir.glob("*.d"))
        if max_files > 0:
            d_files = d_files[:max_files]

        input_files = d_files
        if not input_files:
            return {"status": "skipped", "reason": "No .d files found"}

        # Build Sage command - requires config JSON as first positional arg
        sage_config = Path(__file__).parent.parent.parent / "config" / "sage_config.json"
        cmd = [
            str(sage_path),
            str(sage_config),
            "--fasta", str(fasta_path),
            "--output_directory", str(output_dir),
            "--batch-size", str(max(1, num_threads // 2)),
            "--parquet",
            "--annotate-matches",
        ]

        # Add input files as positional arguments
        cmd.extend([str(f) for f in input_files])

        result = subprocess.run(cmd, capture_output=True, text=True, timeout=3600 * 2)

        if result.returncode != 0:
            return {
                "status": "error",
                "error": result.stderr[:500],
            }

        # Find output parquet
        results_file = output_dir / "results.sage.parquet"
        if results_file.exists():
            import pandas as pd
            df = pd.read_parquet(results_file)
            n_psms = len(df[~df.get("is_decoy", False)])
        else:
            n_psms = 0

        return {
            "status": "success",
            "n_psms": n_psms,
            "output_dir": str(output_dir),
        }

    except subprocess.TimeoutExpired:
        return {"status": "error", "error": "Timeout after 2 hours"}
    except Exception as e:
        return {"status": "error", "error": str(e)}


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
    )

    print(f"\nStep 2 completed: {summary.status}")
    for engine in ["fragpipe", "diann", "sage"]:
        eng_data = summary.data.get(engine, {})
        print(f"  {engine}: {eng_data.get('status')}")
