#!/bin/bash
#SBATCH --account=ki-proanagi
#SBATCH --partition=ki-smallcpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=04:00:00
#SBATCH --job-name=cp-build-imspy
#SBATCH --output=/lustre/project/ki-proanagi/dateschn/logs/build-imspy-%j.log
#
# Build the Python environment for claudius-proteomics on Mogon NHR:
#   - .venv with the scientific Python stack
#   - imspy / imspy_connector / imspy-predictors built from rustims source
#
# Submit with:  sbatch scripts/cluster/build_imspy_env.sbatch

set -euo pipefail

ROOT=/lustre/project/ki-proanagi/dateschn
export CARGO_HOME="$ROOT/.cargo"
export PIP_CACHE_DIR="$ROOT/.pip-cache"
# GCC 13 emits DWARF5 but the system assembler (binutils 2.30) rejects it.
# Suppress C/asm debug info so cc-rs crate builds succeed.
export CFLAGS="-g0"
export CXXFLAGS="-g0"
mkdir -p "$ROOT/logs" "$PIP_CACHE_DIR"

module purge
module load lang/Python/3.12.3-GCCcore-13.3.0
module load lang/Rust/1.91.1-GCCcore-13.3.0

cd "$ROOT/claudius-proteomics"

echo "############ Creating venv ############"
[ -d .venv ] || python3 -m venv .venv
# shellcheck disable=SC1091
source .venv/bin/activate
pip install --upgrade pip wheel setuptools

echo "############ Installing scientific Python stack ############"
pip install numpy pandas pyarrow scipy matplotlib seaborn requests pyyaml \
    zarr tqdm h5py scikit-learn snakemake

echo "############ Building imspy from rustims source ############"
export INSTALL_DIR="$ROOT/rustims"
export RUSTIMS_BRANCH="main"
# Clear the clone only if it is missing or partial (keeps Cargo build cache
# across resubmits for a fast incremental rebuild).
if [ ! -d "$INSTALL_DIR/.git" ]; then
    rm -rf "$INSTALL_DIR"
fi
bash scripts/setup_imspy.sh

echo "############ Verifying ############"
python3 -c "import imspy_connector, pyarrow, pandas, numpy, scipy, yaml; print('core imports OK')"
python3 -c "import imspy_predictors; print('imspy_predictors OK')"
python3 -c "from imspy_core.timstof import TimsDatasetDDA; print('TimsDatasetDDA OK')"
echo "############ DONE ############"
