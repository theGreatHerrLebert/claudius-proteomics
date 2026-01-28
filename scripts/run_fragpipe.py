#!/usr/bin/env python3
"""
Run FragPipe on timsTOF DDA data in headless mode.

Based on the fragpipe_executor pattern from timsim-validate.
"""

import argparse
import os
import subprocess
import sys
from pathlib import Path


def create_manifest(raw_files: list[Path], output_path: Path, data_type: str = "DDA") -> None:
    """
    Create FragPipe manifest file.

    Format: path\texperiment\tbiorep\ttechrep\tdata_type
    """
    with open(output_path, 'w') as f:
        for raw_file in raw_files:
            exp_name = raw_file.stem  # filename without .d
            f.write(f"{raw_file}\t{exp_name}\t1\t1\t{data_type}\n")
    print(f"Created manifest with {len(raw_files)} files: {output_path}")


def update_workflow(
    base_workflow: Path,
    output_workflow: Path,
    fasta_path: Path,
    is_timstof: bool = True,
) -> None:
    """
    Update workflow file with FASTA path and timsTOF settings.
    """
    with open(base_workflow) as f:
        lines = f.readlines()

    updated_lines = []
    found_db_path = False

    for line in lines:
        line = line.rstrip('\n')

        # Update FASTA path
        if line.startswith('database.db-path='):
            updated_lines.append(f'database.db-path={fasta_path}')
            found_db_path = True
        # Enable ion mobility for timsTOF
        elif is_timstof and line.startswith('workflow.input.data-type.im-ms='):
            updated_lines.append('workflow.input.data-type.im-ms=true')
        elif is_timstof and line.startswith('workflow.input.data-type.regular-ms='):
            updated_lines.append('workflow.input.data-type.regular-ms=false')
        else:
            updated_lines.append(line)

    # Add database path if not found
    if not found_db_path:
        # Insert near top, after any header comments
        insert_idx = 0
        for i, line in enumerate(updated_lines):
            if line.startswith('#'):
                insert_idx = i + 1
            else:
                break
        updated_lines.insert(insert_idx, f'database.db-path={fasta_path}')

    with open(output_workflow, 'w') as f:
        f.write('\n'.join(updated_lines))

    print(f"Created workflow: {output_workflow}")
    print(f"  FASTA: {fasta_path}")
    print(f"  timsTOF IM-MS: {is_timstof}")


def find_raw_files(input_dir: Path, max_files: int = 0) -> list[Path]:
    """Find all .d folders in input directory."""
    raw_files = sorted(input_dir.glob("*.d"))

    if not raw_files:
        # Check one level deeper
        raw_files = sorted(input_dir.glob("*/*.d"))

    if max_files > 0:
        raw_files = raw_files[:max_files]
        print(f"Limited to {max_files} files (test mode)")

    print(f"Found {len(raw_files)} raw files")
    return raw_files


def run_fragpipe(
    fragpipe_path: Path,
    workflow_path: Path,
    manifest_path: Path,
    output_dir: Path,
    threads: int = 16,
    ram: int = 0,
    temp_dir: Path | None = None,
) -> int:
    """Run FragPipe in headless mode."""
    fragpipe_bin = fragpipe_path / "bin" / "fragpipe"

    if not fragpipe_bin.exists():
        print(f"ERROR: FragPipe binary not found: {fragpipe_bin}")
        sys.exit(1)

    cmd = [
        str(fragpipe_bin),
        "--headless",
        "--workflow", str(workflow_path),
        "--manifest", str(manifest_path),
        "--workdir", str(output_dir),
    ]

    if threads > 0:
        cmd.extend(["--threads", str(threads)])

    if ram > 0:
        cmd.extend(["--ram", str(ram)])

    print(f"\nRunning FragPipe:")
    print(f"  Command: {' '.join(cmd)}")

    # Set up environment with temp directory for MSFragger cache files
    env = os.environ.copy()
    if temp_dir:
        temp_dir.mkdir(parents=True, exist_ok=True)
        # Set Java temp directory to avoid writing mzBIN cache next to raw files
        java_opts = env.get("_JAVA_OPTIONS", "")
        java_opts = f"{java_opts} -Djava.io.tmpdir={temp_dir}".strip()
        env["_JAVA_OPTIONS"] = java_opts
        env["TMPDIR"] = str(temp_dir)
        env["TEMP"] = str(temp_dir)
        env["TMP"] = str(temp_dir)
        print(f"  Temp dir: {temp_dir}")

    print()

    result = subprocess.run(cmd, cwd=output_dir, env=env)
    return result.returncode


def main():
    parser = argparse.ArgumentParser(
        description="Run FragPipe on timsTOF DDA data"
    )
    parser.add_argument(
        "--fragpipe", required=True,
        help="Path to FragPipe installation directory"
    )
    parser.add_argument(
        "--input", required=True,
        help="Input directory containing .d folders"
    )
    parser.add_argument(
        "--output", required=True,
        help="Output directory for results"
    )
    parser.add_argument(
        "--fasta", required=True,
        help="Path to FASTA database"
    )
    parser.add_argument(
        "--workflow", default="LFQ-MBR",
        help="Workflow name or path to .workflow file (default: LFQ-MBR)"
    )
    parser.add_argument(
        "--threads", type=int, default=16,
        help="Number of threads (default: 16)"
    )
    parser.add_argument(
        "--ram", type=int, default=0,
        help="RAM limit in GB (0 = auto, default: 0)"
    )
    parser.add_argument(
        "--max-files", type=int, default=0,
        help="Maximum number of files to process (0 = all, default: 0)"
    )
    parser.add_argument(
        "--no-timstof", action="store_true",
        help="Disable timsTOF ion mobility settings"
    )
    parser.add_argument(
        "--temp-dir", type=str, default=None,
        help="Temp directory for MSFragger cache files (avoids writing next to raw data)"
    )

    args = parser.parse_args()

    # Resolve paths
    fragpipe_path = Path(args.fragpipe).resolve()
    input_dir = Path(args.input).resolve()
    output_dir = Path(args.output).resolve()
    fasta_path = Path(args.fasta).resolve()

    # Validate
    if not fragpipe_path.exists():
        print(f"ERROR: FragPipe not found: {fragpipe_path}")
        sys.exit(1)

    if not input_dir.exists():
        print(f"ERROR: Input directory not found: {input_dir}")
        sys.exit(1)

    if not fasta_path.exists():
        print(f"ERROR: FASTA not found: {fasta_path}")
        sys.exit(1)

    # Find raw files
    raw_files = find_raw_files(input_dir, args.max_files)
    if not raw_files:
        print(f"ERROR: No .d folders found in {input_dir}")
        sys.exit(1)

    # Create output directory
    output_dir.mkdir(parents=True, exist_ok=True)

    # Resolve workflow path
    if Path(args.workflow).exists():
        base_workflow = Path(args.workflow)
    else:
        # Look in FragPipe workflows directory
        base_workflow = fragpipe_path / "workflows" / f"{args.workflow}.workflow"
        if not base_workflow.exists():
            print(f"ERROR: Workflow not found: {args.workflow}")
            print(f"  Checked: {base_workflow}")
            print(f"  Available workflows:")
            for wf in (fragpipe_path / "workflows").glob("*.workflow"):
                print(f"    {wf.stem}")
            sys.exit(1)

    print(f"Using base workflow: {base_workflow}")

    # Create manifest
    manifest_path = output_dir / "manifest.fp-manifest"
    create_manifest(raw_files, manifest_path, "DDA")

    # Create updated workflow
    workflow_path = output_dir / "workflow.workflow"
    update_workflow(
        base_workflow,
        workflow_path,
        fasta_path,
        is_timstof=not args.no_timstof,
    )

    # Resolve temp directory
    temp_dir = Path(args.temp_dir).resolve() if args.temp_dir else None

    # Run FragPipe
    return_code = run_fragpipe(
        fragpipe_path,
        workflow_path,
        manifest_path,
        output_dir,
        args.threads,
        args.ram,
        temp_dir,
    )

    if return_code == 0:
        print("\nFragPipe completed successfully!")

        # Check for output files
        psm_file = output_dir / "psm.tsv"
        if psm_file.exists():
            print(f"  PSM file: {psm_file}")
        else:
            # Check in experiment subdirectory
            for subdir in output_dir.iterdir():
                if subdir.is_dir() and (subdir / "psm.tsv").exists():
                    print(f"  PSM file: {subdir / 'psm.tsv'}")
                    break
    else:
        print(f"\nFragPipe failed with exit code {return_code}")

    sys.exit(return_code)


if __name__ == "__main__":
    main()
