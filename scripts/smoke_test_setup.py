#!/usr/bin/env python3
"""Generate 3 simulated DDA timsTOF datasets for San Jose smoke testing.

Creates SMOKE_REP_01 (full simulation), SMOKE_REP_02 and SMOKE_REP_03
(variation replicates from REP_01 with RT/IM/intensity noise).

Usage:
    .venv/bin/python scripts/smoke_test_setup.py          # Generate all 3
    .venv/bin/python scripts/smoke_test_setup.py --check   # Check if ready
"""

import argparse
import subprocess
import sys
import tempfile
from pathlib import Path

REFERENCE_D = Path.home() / "validate-sim/sim/blanks/ref-dda/G230913_002_Slot2-2_1_11392.d"
FASTA_PATH = Path.home() / "validate-sim/sim/hela.fasta"
OUTPUT_DIR = Path(__file__).resolve().parent.parent / "data" / "smoke_test" / "simulations"

REPLICATE_NAMES = ["SMOKE_REP_01", "SMOKE_REP_02", "SMOKE_REP_03"]


def generate_base_config(name: str, save_path: str) -> str:
    """Config for SMOKE_REP_01: full simulation, no variation."""
    return f"""[paths]
save_path = "{save_path}"
reference_path = "{REFERENCE_D}"
fasta_path = "{FASTA_PATH}"

[experiment]
experiment_name = "{name}"
acquisition_type = "DDA"
gradient_length = 3600.0
use_reference_layout = true
reference_in_memory = false
use_bruker_sdk = true
apply_fragmentation = true
from_existing = false
silent_mode = true

[digestion]
n_proteins = 20000
num_peptides_total = 250000
num_sample_peptides = 150000
sample_peptides = true
sample_seed = 42
cleave_at = "KR"
restrict = "P"
missed_cleavages = 2
min_len = 7
max_len = 30
decoys = false
upscale_factor = 100000

[retention_time]
sigma_alpha_rt = 1
sigma_beta_rt = 1
k_lower_rt = 0.1
k_upper_rt = 2
k_alpha_rt = 1
k_beta_rt = 1
target_p = 0.999
sampling_step_size = 0.0001

[ion_mobility]
use_inverse_mobility_std_mean = false
inverse_mobility_std_mean = 0.0075

[charge_states]
p_charge = 0.5
min_charge_contrib = 0.25

[isotopic_pattern]
isotope_k = 8
isotope_min_intensity = 1
isotope_centroid = true

[fragment_intensity]
down_sample_factor = 0.5

[noise]
mz_noise_precursor = true
precursor_noise_ppm = 6.5
mz_noise_fragment = true
fragment_noise_ppm = 6.5
mz_noise_uniform = false
add_real_data_noise = true
reference_noise_intensity_max = 9999999

[dda]
precursors_every = 7
max_precursors = 8
exclusion_width = 25
precursor_intensity_threshold = 1000
selection_mode = "topN"

[performance]
num_threads = -1
batch_size = 256
use_gpu = false
"""


def generate_variation_config(name: str, save_path: str, existing_path: str) -> str:
    """Config for SMOKE_REP_02/03: reuse REP_01, apply variation."""
    return f"""[paths]
save_path = "{save_path}"
reference_path = "{REFERENCE_D}"
fasta_path = "{FASTA_PATH}"
existing_path = "{existing_path}"

[experiment]
experiment_name = "{name}"
acquisition_type = "DDA"
gradient_length = 3600.0
use_reference_layout = true
reference_in_memory = false
use_bruker_sdk = true
apply_fragmentation = true
from_existing = true
silent_mode = true

[variation]
re_scale_rt = true
rt_variation_std = 7
ion_mobility_variation_std = 0.01
intensity_variation_std = 0.02

[isotopic_pattern]
isotope_k = 8
isotope_min_intensity = 1
isotope_centroid = true

[fragment_intensity]
down_sample_factor = 0.5

[noise]
mz_noise_precursor = true
precursor_noise_ppm = 6.5
mz_noise_fragment = true
fragment_noise_ppm = 6.5
mz_noise_uniform = false
add_real_data_noise = true
reference_noise_intensity_max = 9999999

[dda]
precursors_every = 7
max_precursors = 8
exclusion_width = 25
precursor_intensity_threshold = 1000
selection_mode = "topN"

[performance]
num_threads = -1
batch_size = 256
use_gpu = false
"""


def d_folder_exists(name: str) -> bool:
    """Check if a simulation's .d folder exists."""
    return (OUTPUT_DIR / name / name / f"{name}.d").exists()


def db_exists(name: str) -> bool:
    """Check if a simulation's ground truth database exists."""
    return (OUTPUT_DIR / name / name / "synthetic_data.db").exists()


def run_simulation(config_text: str, name: str) -> None:
    """Write config to temp file and run timsim."""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".toml", delete=False) as f:
        f.write(config_text)
        config_path = f.name

    print(f"  Running simulation for {name}...")
    print(f"  Config: {config_path}")

    result = subprocess.run(
        [sys.executable, "-m", "imspy_simulation.timsim.simulator", config_path],
        capture_output=False,
        text=True,
    )

    Path(config_path).unlink(missing_ok=True)

    if result.returncode != 0:
        raise RuntimeError(f"Simulation failed for {name} (exit code {result.returncode})")


def check_simulations() -> bool:
    """Check if all 3 simulations are ready."""
    all_ready = True
    for name in REPLICATE_NAMES:
        has_d = d_folder_exists(name)
        has_db = db_exists(name)
        status = "ready" if (has_d and has_db) else "MISSING"
        print(f"  {name}: {status}")
        if not has_d:
            print(f"    .d folder: {'found' if has_d else 'NOT FOUND'}")
        if not has_db:
            print(f"    synthetic_data.db: {'found' if has_db else 'NOT FOUND'}")
        if not (has_d and has_db):
            all_ready = False
    return all_ready


def main():
    parser = argparse.ArgumentParser(description="Generate simulated datasets for smoke testing")
    parser.add_argument("--check", action="store_true", help="Check if simulations exist (don't generate)")
    args = parser.parse_args()

    if args.check:
        print("Checking smoke test simulations...")
        ready = check_simulations()
        print(f"\nAll simulations ready: {ready}")
        sys.exit(0 if ready else 1)

    # Validate prerequisites
    if not REFERENCE_D.exists():
        print(f"Error: Reference .d not found: {REFERENCE_D}")
        sys.exit(1)
    if not FASTA_PATH.exists():
        print(f"Error: FASTA not found: {FASTA_PATH}")
        sys.exit(1)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # Step 1: Generate base simulation (SMOKE_REP_01)
    rep1_name = "SMOKE_REP_01"
    if d_folder_exists(rep1_name) and db_exists(rep1_name):
        print(f"{rep1_name}: already exists, skipping")
    else:
        print(f"\nGenerating {rep1_name} (full simulation, ~10-15 min)...")
        rep1_dir = OUTPUT_DIR / rep1_name
        rep1_dir.mkdir(parents=True, exist_ok=True)
        config = generate_base_config(rep1_name, str(rep1_dir))
        run_simulation(config, rep1_name)

        if not d_folder_exists(rep1_name):
            print(f"Error: {rep1_name} simulation did not produce .d folder")
            sys.exit(1)
        print(f"  {rep1_name}: done")

    # Step 2: Generate variation replicates (SMOKE_REP_02, 03)
    existing_path = str(OUTPUT_DIR / rep1_name / rep1_name)
    for name in ["SMOKE_REP_02", "SMOKE_REP_03"]:
        if d_folder_exists(name) and db_exists(name):
            print(f"{name}: already exists, skipping")
            continue

        print(f"\nGenerating {name} (variation replicate, ~5 min)...")
        rep_dir = OUTPUT_DIR / name
        rep_dir.mkdir(parents=True, exist_ok=True)
        config = generate_variation_config(name, str(rep_dir), existing_path)
        run_simulation(config, name)

        if not d_folder_exists(name):
            print(f"Error: {name} simulation did not produce .d folder")
            sys.exit(1)
        print(f"  {name}: done")

    print("\nAll simulations generated successfully.")
    check_simulations()


if __name__ == "__main__":
    main()
