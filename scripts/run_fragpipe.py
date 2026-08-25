#!/usr/bin/env python3
"""
Run FragPipe on timsTOF DDA data in headless mode.

Based on the fragpipe_executor pattern from timsim-validate.
"""

import argparse
import json
import os
import subprocess
import sys
from pathlib import Path
from typing import Any, Dict, Optional


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
    disable_fdr_filter: bool = True,
    enzyme_name: str | None = None,
    mod_config: Optional[Dict[str, Any]] = None,
    skip_msbooster: bool = False,
) -> None:
    """
    Update workflow file with FASTA path and timsTOF settings.

    If disable_fdr_filter=True (San José mode), sets all FDR thresholds to 1.0
    so that ALL PSMs are reported regardless of confidence.

    If enzyme_name is provided, override the enzyme name in the workflow
    (msfragger.search_enzyme_name_1, or the legacy unsuffixed key). Raises if
    the key is absent or if the base workflow's specificity
    (num_enzyme_termini / search_enzyme_cut_1) contradicts the requested
    enzyme — a name-only override cannot change MSFragger's digestion.

    If mod_config is provided, override MSFragger modification settings
    from the canonical UNIMOD-based profile.
    """
    with open(base_workflow) as f:
        lines = f.readlines()

    updated_lines = []
    found_db_path = False

    # FDR-related settings to override for San José (report ALL results).
    #
    # phi-report.filter is the literal `philosopher filter` command line.
    # Setting --psm/--pep/--ion to 1 makes Philosopher report ALL PSMs.
    #
    # The --prot value, however, is NOT honored: FragPipe's CmdPhilosopherFilter
    # hardcodes `--prot 0.01` (regex-replacing whatever we pass) whenever
    # IonQuant runs — see CmdPhilosopherFilter.configure(), gated on
    # isRunIonQuant. With --sequential, that 1% protein FDR cascades down and
    # silently FDR-filters the PSM table (~500k -> ~40k observed on PXD046777).
    # So we also disable label-free quant / IonQuant below: that flips
    # isRunIonQuant to false, the --prot override is skipped, and our --prot 1
    # is used verbatim. San José does its own raw 4D extraction downstream and
    # never consumes IonQuant's MS1 quant, so this loses nothing.
    #
    # NOTE: keys must exist in the base workflow to take effect — FragPipe
    # silently ignores unknown keys. Every key below is present in the stock
    # LFQ-MBR workflow.
    fdr_overrides = {
        'msfragger.report-fdr': '1.0',
        'percolator.min-prob': '0.0',
        'ptmprophet.fdr': '1.0',
        'phi-report.filter': '--sequential --picked --psm 1 --pep 1 --ion 1 --prot 1 --pepProb 0 --protProb 0',
        # Disable IonQuant/LFQ so FragPipe does not force `--prot 0.01`.
        'quantitation.run-label-free-quant': 'false',
        'ionquant.run-ionquant': 'false',
    }

    applied_overrides = set()

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
        # Override FDR settings if San José mode
        elif disable_fdr_filter:
            key = line.split('=')[0] if '=' in line else None
            if key and key in fdr_overrides:
                updated_lines.append(f'{key}={fdr_overrides[key]}')
                applied_overrides.add(key)
            else:
                updated_lines.append(line)
        else:
            updated_lines.append(line)

    # Append any FDR overrides not already in the base workflow
    if disable_fdr_filter:
        for key, value in fdr_overrides.items():
            if key not in applied_overrides:
                updated_lines.append(f'{key}={value}')

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

    # Override enzyme if requested.
    #
    # This used to match `msfragger.search_enzyme_name=`, but FragPipe >= 20
    # keys the enzyme per slot (`msfragger.search_enzyme_name_1`), so the
    # unsuffixed key matched nothing in any shipped workflow: the override was
    # a silent no-op and nothing checked that it applied.
    #
    # Matching the right key is necessary but NOT sufficient. MSFragger derives
    # specificity from `search_enzyme_cut_1` / `num_enzyme_termini`, not from
    # the name, so writing the name alone would leave the workflow internally
    # inconsistent (name=nonspecific, cut=KR, termini=2) — worse than the
    # no-op. Rewriting digestion semantics here would silently change results
    # for every profile, so instead: require the base workflow to already agree
    # with the requested enzyme, and fail loudly when it does not. The
    # supported way to get a non-tryptic search is to select the matching base
    # workflow (`mod_profile.fragpipe_workflow`, e.g. `Nonspecific-HLA`).
    if enzyme_name:
        enzyme_keys = ('msfragger.search_enzyme_name_1=',   # FragPipe >= 20
                       'msfragger.search_enzyme_name=')     # legacy
        current = {}
        for line in updated_lines:
            for key in ('msfragger.num_enzyme_termini=', 'msfragger.search_enzyme_cut_1='):
                if line.startswith(key):
                    current[key] = line.split('=', 1)[1].strip()

        wants_nonspecific = enzyme_name.strip().lower() == 'nonspecific'
        termini = current.get('msfragger.num_enzyme_termini=')
        cut = current.get('msfragger.search_enzyme_cut_1=')
        workflow_is_nonspecific = termini == '0' or cut == '@'

        if wants_nonspecific != workflow_is_nonspecific:
            raise ValueError(
                f"enzyme override '{enzyme_name}' does not match base workflow "
                f"{Path(base_workflow).name} (num_enzyme_termini={termini}, "
                f"search_enzyme_cut_1={cut}). A name-only override cannot change "
                f"MSFragger's specificity — select a base workflow that matches "
                f"(mod_profile.fragpipe_workflow)."
            )

        replaced = False
        final_lines = []
        for line in updated_lines:
            matched = next((k for k in enzyme_keys if line.startswith(k)), None)
            if matched:
                final_lines.append(f'{matched}{enzyme_name}')
                replaced = True
            else:
                final_lines.append(line)
        if not replaced:
            raise ValueError(
                f"enzyme override '{enzyme_name}' requested but no enzyme key "
                f"({' or '.join(k.rstrip('=') for k in enzyme_keys)}) exists in "
                f"{base_workflow} — refusing to report a no-op as success."
            )
        updated_lines = final_lines

    # Override modifications if mod_config provided
    if mod_config:
        mod_overrides = _build_msfragger_mod_overrides(mod_config)
        final_lines = []
        for line in updated_lines:
            key = line.split('=')[0] if '=' in line else None
            if key and key in mod_overrides:
                final_lines.append(f'{key}={mod_overrides.pop(key)}')
            else:
                final_lines.append(line)
        # Append any overrides not found in existing lines
        for key, value in mod_overrides.items():
            final_lines.append(f'{key}={value}')
        updated_lines = final_lines

    # Disable MSBooster for datasets with raw-filename spaces (NoSuchFileException
    # in PinMzmlMatcher) or mzML scan-null NPEs in MzmlReader.getScanNumObject.
    # FragPipe still runs MSFragger + Philosopher; we lose the booster's rescoring.
    if skip_msbooster:
        found = False
        final_lines = []
        for line in updated_lines:
            if line.startswith('msbooster.run-msbooster='):
                final_lines.append('msbooster.run-msbooster=false')
                found = True
            else:
                final_lines.append(line)
        if not found:
            final_lines.append('msbooster.run-msbooster=false')
        updated_lines = final_lines

    with open(output_workflow, 'w') as f:
        f.write('\n'.join(updated_lines))

    print(f"Created workflow: {output_workflow}")
    print(f"  FASTA: {fasta_path}")
    print(f"  timsTOF IM-MS: {is_timstof}")
    if enzyme_name:
        print(f"  Enzyme: {enzyme_name}")
    if mod_config:
        print(f"  Mod profile: {mod_config.get('description', 'custom')}")
    if disable_fdr_filter:
        print(f"  FDR filtering: DISABLED (San José mode - report ALL PSMs)")
    if skip_msbooster:
        print(f"  MSBooster: DISABLED (workaround for filename-space / mzML NPE bugs)")


def _build_msfragger_mod_overrides(mod_config: Dict[str, Any]) -> Dict[str, str]:
    """Build MSFragger workflow key overrides from a mod profile.

    The authoritative variable-mod set is the FragPipe UI table
    ``msfragger.table.var-mods`` — FragPipe regenerates MSFragger's
    ``variable_mod_*`` params from it, so writing ``variable_mod_0X`` directly has
    no effect. Table entry format: ``mass,residues,enabled,max_per_peptide``,
    with ``[^`` marking the N-terminus.

    Returns:
        Dict of workflow keys -> values (``table.var-mods``, digest lengths, and
        any profile ``msfragger_overrides``).
    """
    overrides = {}

    var_mods = mod_config.get("variable_modifications", [])
    max_var = mod_config.get("max_variable_mods", 3)
    if len(var_mods) > 16:
        raise ValueError(
            f"mod profile has {len(var_mods)} variable modifications; FragPipe's "
            f"msfragger.table.var-mods supports at most 16 slots"
        )

    # Peptide length
    min_len = mod_config.get("min_peptide_length", 7)
    max_len = mod_config.get("max_peptide_length", 50)
    overrides["msfragger.digest_min_length"] = str(min_len)
    overrides["msfragger.digest_max_length"] = str(max_len)

    # --- A4: build the authoritative variable-mod table FragPipe actually uses ---
    # Profile mods ENABLED; padded to 16 slots with disabled placeholders to match
    # FragPipe's fixed table size. N-term mods use the [^ marker and max 1 (a
    # peptide has a single N-terminus); side-chain mods use the profile's
    # max_variable_mods as the per-mod occurrence cap.
    table = []
    for vm in var_mods:
        if vm.get("site") == "N-term":
            table.append(f"{vm['mass']:.6f},[^,true,1")
        else:
            residues = "".join(vm["residues"])
            table.append(f"{vm['mass']:.6f},{residues},true,{max_var}")
    for i in range(len(table) + 1, 17):
        table.append(f"0.0,site_{i:02d},false,1")
    overrides["msfragger.table.var-mods"] = "; ".join(table)

    # Generic MSFragger key overrides from mod_profile (e.g. num_database_chunks
    # for HLA nonspecific searches that would otherwise OOM). Values are
    # written verbatim to the workflow as `msfragger.<key>=<value>`.
    extra = mod_config.get("msfragger_overrides", {}) or {}
    for key, value in extra.items():
        overrides[f"msfragger.{key}"] = str(value)

    return overrides


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
    parser.add_argument(
        "--enable-fdr-filter", action="store_true",
        help="Enable FDR filtering (default: disabled for San José - report ALL PSMs)"
    )
    parser.add_argument(
        "--files", nargs="+", type=str, default=None,
        help="Explicit list of .d file paths (bypasses find_raw_files)"
    )
    parser.add_argument(
        "--enzyme", type=str, default=None,
        help="Enzyme name to override in workflow (e.g. stricttrypsin, lysc, lysn)"
    )
    parser.add_argument(
        "--mod-config-json", type=str, default=None,
        help="Modification profile as JSON string (from mod_profiles config)"
    )
    parser.add_argument(
        "--skip-msbooster", action="store_true",
        help="Disable MSBooster in the workflow. Use for datasets where MSBooster "
             "fails (raw filenames with spaces → PinMzmlMatcher NoSuchFileException, "
             "or mzML scan NPE in MzmlReader)."
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

    # Find raw files: prefer explicit --files list, fall back to directory scan
    if args.files:
        raw_files = [Path(f).resolve() for f in args.files]
        print(f"Using {len(raw_files)} explicitly provided files")
    else:
        raw_files = find_raw_files(input_dir, args.max_files)
    if not raw_files:
        print(f"ERROR: No .d folders found in {input_dir}")
        sys.exit(1)

    # Create output directory
    output_dir.mkdir(parents=True, exist_ok=True)

    # Sanitize raw filenames: FragPipe's MSBooster (PinMzmlMatcher) and
    # Percolator (PercolatorOutputToPepXML.getSpectrum) both parse paths on
    # whitespace and crash on raw filenames containing spaces. Symlink each
    # offending .d under a space-free name and feed FragPipe the symlinks.
    if any(' ' in p.name for p in raw_files):
        nosp_dir = output_dir / "_raw_nosp"
        nosp_dir.mkdir(parents=True, exist_ok=True)
        sanitized = []
        for p in raw_files:
            new_name = p.name.replace(' ', '_')
            sym = nosp_dir / new_name
            if sym.is_symlink() or sym.exists():
                sym.unlink()
            sym.symlink_to(p.resolve())
            sanitized.append(sym)
        print(f"Sanitized {len(sanitized)} raw paths (spaces → '_') in {nosp_dir}")
        raw_files = sanitized

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

    # Parse mod config if provided
    mod_config_dict = None
    if args.mod_config_json:
        try:
            mod_config_dict = json.loads(args.mod_config_json)
        except json.JSONDecodeError as e:
            print(f"WARNING: Could not parse --mod-config-json: {e}")

    # Create updated workflow
    workflow_path = output_dir / "workflow.workflow"
    update_workflow(
        base_workflow,
        workflow_path,
        fasta_path,
        is_timstof=not args.no_timstof,
        disable_fdr_filter=not args.enable_fdr_filter,  # Default: disable FDR (San José mode)
        enzyme_name=args.enzyme,
        mod_config=mod_config_dict,
        skip_msbooster=args.skip_msbooster,
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
