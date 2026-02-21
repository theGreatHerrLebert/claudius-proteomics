#!/usr/bin/env python3
"""
Step 2: Execute Third-Party Search Engines

Runs FragPipe, DIA-NN, and Sage on raw data with FDR=1.0 (full results).
Each engine is a self-contained job (runner/engines/) that produces
canonical parquet + status JSON.

Supports multi-organism datasets: loads sample_groups.yaml from step 1 and
runs each engine per sample group with the correct FASTA and enzyme settings.
Falls back to single-group (legacy) behavior when no sample_groups.yaml exists.

Input: data/raw/{accession}/*.d, FASTA database(s)

Outputs (per-group):
- data/processed/{accession}/{group_id}/{engine}_canonical.parquet
- data/processed/{accession}/{group_id}/{engine}_status.json
- data/processed/{accession}/step2_summary.json
"""

import sys
from pathlib import Path
from typing import Dict, Any, List, Optional

# Add parent to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from runner.engines import ENGINES
from runner.summary import StepSummary, write_step_summary
from scripts.sample_group_resolver import SampleGroupManifest


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

    If sample_groups.yaml exists (written by step 1), runs each engine per sample
    group with the group-specific FASTA and enzyme. Otherwise falls back to
    single-group legacy behavior.

    Args:
        accession: PRIDE accession
        config: Pipeline configuration dict
        output_base_dir: Base directory for outputs
        raw_dir: Path to raw data directory (default: data/raw/{accession})
        fasta_path: Explicit FASTA path (overrides per-group resolution)
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
    metadata_dir = output_base_dir / "metadata" / accession

    try:
        # Try to load sample groups from step 1
        sg_path = metadata_dir / "sample_groups.yaml"
        if sg_path.exists():
            manifest = SampleGroupManifest.from_yaml(sg_path)
            print(f"  Loaded {len(manifest.groups)} sample groups from {sg_path}")
        else:
            manifest = None
            print(f"  No sample_groups.yaml found — using single-group mode")

        if manifest and len(manifest.groups) > 0:
            results = _run_per_group(
                manifest=manifest,
                accession=accession,
                config=config,
                raw_dir=raw_dir,
                processed_dir=processed_dir,
                output_base_dir=output_base_dir,
                engines=engines,
                num_threads=num_threads,
                max_files=max_files,
                fasta_override=fasta_path,
            )
        else:
            results = _run_single_group(
                accession=accession,
                config=config,
                raw_dir=raw_dir,
                processed_dir=processed_dir,
                output_base_dir=output_base_dir,
                engines=engines,
                num_threads=num_threads,
                max_files=max_files,
                fasta_path=fasta_path,
            )

        # Check for engine errors (raises only when ALL engines failed)
        engine_check = _check_engine_errors(results)
        results["_engine_check"] = engine_check

        # Update summary
        summary.data = results
        summary.outputs = [str(processed_dir)]
        summary.complete(success=True)

    except Exception as e:
        summary.complete(success=False, error_message=str(e))
        raise

    # Write summary file
    write_step_summary(summary, processed_dir)

    return summary


def _check_engine_errors(results: Dict[str, Any]) -> Dict[str, Any]:
    """Check engine statuses and return a summary.

    Tolerates partial failures: only raises RuntimeError when ALL engines
    failed.  When at least one engine succeeded, prints warnings for the
    failed engines and returns normally so that Steps 3-6 can proceed.

    Returns a dict with keys: succeeded, failed, skipped (each a list of
    engine labels), plus a boolean ``all_failed``.
    """
    succeeded: List[str] = []
    failed: List[str] = []
    skipped: List[str] = []
    error_details: List[str] = []

    if results.get("mode") == "per_group":
        for group_id, group_data in results.get("groups", {}).items():
            engines = group_data.get("engines", {})
            for engine_name, engine_data in engines.items():
                label = f"{group_id}/{engine_name}"
                status = engine_data.get("status")
                if status == "error":
                    msg = engine_data.get("error_message", "unknown error")
                    failed.append(label)
                    error_details.append(f"{label}: {msg}")
                elif status == "skipped":
                    skipped.append(label)
                else:
                    succeeded.append(label)
    else:
        for engine_name in ["fragpipe", "diann", "sage"]:
            engine_data = results.get(engine_name, {})
            status = engine_data.get("status")
            if status == "error":
                msg = engine_data.get("error_message", "unknown error")
                failed.append(engine_name)
                error_details.append(f"{engine_name}: {msg}")
            elif status == "skipped":
                skipped.append(engine_name)
            else:
                succeeded.append(engine_name)

    check = {
        "succeeded": succeeded,
        "failed": failed,
        "skipped": skipped,
        "all_failed": len(succeeded) == 0 and len(failed) > 0,
    }

    if failed:
        for detail in error_details:
            print(f"  WARNING: engine error — {detail}")

    if check["all_failed"]:
        raise RuntimeError(
            f"All {len(failed)} engine(s) failed — cannot continue: "
            + "; ".join(error_details)
        )

    if failed:
        print(
            f"  Step 2 partial success: {len(succeeded)} engine(s) OK, "
            f"{len(failed)} failed, {len(skipped)} skipped"
        )

    return check


def _run_per_group(
    manifest: SampleGroupManifest,
    accession: str,
    config: Dict[str, Any],
    raw_dir: Path,
    processed_dir: Path,
    output_base_dir: Path,
    engines: List[str],
    num_threads: int,
    max_files: int,
    fasta_override: Optional[Path] = None,
) -> Dict[str, Any]:
    """Run engines for each sample group."""
    group_results = {}

    for group in manifest.groups:
        print(f"\n  === Sample group: {group.group_id} ===")
        print(f"      Organism:    {group.organism_name} ({group.organism_key})")
        print(f"      Enzyme:      {group.enzyme}")
        print(f"      Instrument:  {group.instrument_model or 'unknown'}")
        print(f"      Runs:        {group.n_runs}")

        if group.n_runs == 0:
            print(f"      Skipping (no runs)")
            group_results[group.group_id] = {"status": "skipped", "reason": "no runs"}
            continue

        # Resolve FASTA for this group's organism
        if fasta_override:
            group_fasta = fasta_override
        else:
            group_fasta = _get_fasta_for_organism(
                group.organism_key, config, output_base_dir
            )
        print(f"      FASTA:    {group_fasta}")

        # Resolve modification profile
        mod_config = config.get("mod_profiles", {}).get(group.mod_profile)
        if mod_config:
            print(f"      Mod profile:   {group.mod_profile}")
        else:
            print(f"      Mod profile:   {group.mod_profile} (not in config, using defaults)")

        # Resolve enzyme config (may be overridden by mod profile)
        enzyme_name = group.enzyme
        if mod_config and mod_config.get("enzyme_override"):
            enzyme_name = mod_config["enzyme_override"]
            print(f"      Enzyme override: {enzyme_name} (from mod profile)")
        enzyme_config = config.get("enzymes", {}).get(enzyme_name)
        if enzyme_config:
            print(f"      Enzyme config: {enzyme_name}")
        else:
            print(f"      Warning: no enzyme config for '{enzyme_name}', using defaults")

        # Resolve .d file paths
        group_d_files = [raw_dir / run for run in group.runs]
        if max_files > 0:
            group_d_files = group_d_files[:max_files]

        # Per-group output directory
        group_processed_dir = processed_dir / group.group_id

        # Run each engine
        engine_results = {}
        for engine_name in engines:
            if engine_name not in ENGINES:
                engine_results[engine_name] = {
                    "status": "skipped",
                    "reason": f"Unknown engine: {engine_name}",
                }
                continue

            job = ENGINES[engine_name]()
            result = job.run(
                accession=accession,
                config=config,
                raw_dir=raw_dir,
                processed_dir=group_processed_dir,
                fasta_path=group_fasta,
                num_threads=num_threads,
                max_files=0,  # Already sliced via d_files
                d_files=group_d_files,
                enzyme_config=enzyme_config,
                mod_config=mod_config,
            )
            engine_results[engine_name] = result.to_dict()

        # Ensure all default engines have an entry
        for name in ["fragpipe", "diann", "sage"]:
            if name not in engine_results:
                engine_results[name] = {"status": "skipped"}

        group_results[group.group_id] = {
            "organism": group.organism_key,
            "enzyme": group.enzyme,
            "mod_profile": group.mod_profile,
            "n_runs": len(group_d_files),
            "fasta": str(group_fasta),
            "engines": engine_results,
        }

    return {
        "mode": "per_group",
        "n_groups": len(manifest.groups),
        "groups": group_results,
    }


def _run_single_group(
    accession: str,
    config: Dict[str, Any],
    raw_dir: Path,
    processed_dir: Path,
    output_base_dir: Path,
    engines: List[str],
    num_threads: int,
    max_files: int,
    fasta_path: Optional[Path] = None,
) -> Dict[str, Any]:
    """Legacy single-group execution (backward compatible)."""
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

    return {
        "mode": "single",
        "fragpipe": results.get("fragpipe", {"status": "skipped"}),
        "diann": results.get("diann", {"status": "skipped"}),
        "sage": results.get("sage", {"status": "skipped"}),
        "fasta": fasta_info,
    }


def _get_fasta_for_organism(
    organism_key: str,
    config: Dict[str, Any],
    output_base_dir: Path,
) -> Path:
    """
    Get FASTA database path for a specific organism.

    Prefers decoy FASTAs (required by FragPipe/Sage). If only a target FASTA
    exists (local_fasta), auto-generates a decoy version with rev_ prefix.
    """
    resources_dir = Path(__file__).parent.parent.parent / "resources" / "fasta" / "search_db"

    # 1. Check for pre-built decoy FASTA in search_db
    organism_decoy = resources_dir / f"{organism_key}_decoys.fasta"
    if organism_decoy.exists():
        return organism_decoy

    # 2. Check output dir for decoy FASTA
    fasta_dir = output_base_dir / "resources" / "fasta"
    organism_decoy_alt = fasta_dir / f"{organism_key}_decoys.fasta"
    if organism_decoy_alt.exists():
        return organism_decoy_alt

    # 3. No decoy FASTA exists — try to generate from local_fasta
    source_fasta = None
    organisms = config.get("organisms", {})
    if organism_key in organisms:
        org_config = organisms[organism_key]
        if isinstance(org_config, dict):
            local_fasta = org_config.get("local_fasta")
            if local_fasta:
                candidate = Path(local_fasta)
                if candidate.exists():
                    source_fasta = candidate
        elif isinstance(org_config, str):
            candidate = Path(org_config)
            if candidate.exists():
                source_fasta = candidate

    # Also check symlinked FASTA in resources/fasta/
    if source_fasta is None:
        symlink_fasta = Path(__file__).parent.parent.parent / "resources" / "fasta" / f"{organism_key}.fasta"
        if symlink_fasta.exists():
            source_fasta = symlink_fasta

    if source_fasta is not None:
        print(f"      Generating decoy FASTA from {source_fasta.name}...")
        _generate_decoy_fasta(source_fasta, organism_decoy)
        return organism_decoy

    raise FileNotFoundError(
        f"No FASTA database found for organism '{organism_key}'. "
        f"Configure it in config.yaml under organisms.{organism_key}.local_fasta "
        f"or place {organism_key}_decoys.fasta in resources/fasta/search_db/"
    )


def _generate_decoy_fasta(source_fasta: Path, output_fasta: Path) -> None:
    """
    Generate a target+decoy FASTA from a target-only FASTA.

    Appends reversed sequences with 'rev_' header prefix, matching the
    Snakemake add_decoys rule format.
    """
    output_fasta.parent.mkdir(parents=True, exist_ok=True)

    with open(source_fasta) as f_in, open(output_fasta, "w") as f_out:
        sequences = []
        current_header = None
        current_seq = []

        for line in f_in:
            f_out.write(line)  # Write original
            if line.startswith(">"):
                if current_header and current_seq:
                    sequences.append((current_header, "".join(current_seq)))
                current_header = line.strip()
                current_seq = []
            else:
                current_seq.append(line.strip())

        if current_header and current_seq:
            sequences.append((current_header, "".join(current_seq)))

        # Write reversed decoys
        f_out.write("\n")
        for header, seq in sequences:
            decoy_header = f">rev_{header[1:]}"
            reversed_seq = seq[::-1]
            f_out.write(f"{decoy_header}\n")
            for i in range(0, len(reversed_seq), 80):
                f_out.write(f"{reversed_seq[i:i+80]}\n")

    n_target = len(sequences)
    print(f"      Generated {output_fasta.name}: {n_target} targets + {n_target} decoys")


def _get_fasta_path(accession: str, config: Dict[str, Any], output_base_dir: Path) -> Path:
    """Get FASTA database path for accession (legacy single-organism lookup)."""
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
        try:
            return _get_fasta_for_organism(organism, config, output_base_dir)
        except FileNotFoundError:
            pass

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
    mode = summary.data.get("mode", "single")
    if mode == "per_group":
        print(f"  Mode: per-group ({summary.data.get('n_groups', 0)} groups)")
        for gid, gdata in summary.data.get("groups", {}).items():
            eng = gdata.get("engines", {})
            statuses = {e: eng.get(e, {}).get("status", "?") for e in ["fragpipe", "diann", "sage"]}
            print(f"    {gid}: {statuses}")
    else:
        for engine in ["fragpipe", "diann", "sage"]:
            eng_data = summary.data.get(engine, {})
            print(f"  {engine}: {eng_data.get('status')}")
