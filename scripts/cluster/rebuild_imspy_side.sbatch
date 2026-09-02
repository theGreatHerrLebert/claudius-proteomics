#!/bin/bash
#SBATCH --account=ki-proanagi
#SBATCH --partition=ki-smallcpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=04:00:00
#SBATCH --job-name=cp-imspy-side
#SBATCH --output=/lustre/project/ki-proanagi/dateschn/logs/build-imspy-side-%j.log
#
# Build a SIDE rustims/imspy environment, leaving the production .venv alone.
#
# The deployed reader is a 2026-05-20 build and predates the SDK-free Bruker
# calibration work, so we want to measure the new one against it before trusting
# either. That means both have to exist at once: this writes a fresh clone to
# rustims-new/ and a fresh venv to .venv-readernew/, and touches nothing the
# running pipeline imports.
#
# Submit with:  sbatch scripts/cluster/rebuild_imspy_side.sbatch

set -euo pipefail

ROOT=/lustre/project/ki-proanagi/dateschn
P="$ROOT/claudius-proteomics"
SIDE_VENV="$P/.venv-readernew"
export INSTALL_DIR="$ROOT/rustims-new"
export RUSTIMS_BRANCH="main"

export CARGO_HOME="$ROOT/.cargo"
export PIP_CACHE_DIR="$ROOT/.pip-cache"
# GCC 13 emits DWARF5 but binutils 2.30 rejects it; kill debug info so the
# cc-rs crate builds go through (same reason as build_imspy_env.sbatch).
export CFLAGS="-g0"
export CXXFLAGS="-g0"
mkdir -p "$PIP_CACHE_DIR"

module purge
module load lang/Python/3.12.3-GCCcore-13.3.0
module load lang/Rust/1.91.1-GCCcore-13.3.0

cd "$P"

echo "############ Side venv: $SIDE_VENV ############"
[ -d "$SIDE_VENV" ] || python3 -m venv "$SIDE_VENV"
# shellcheck disable=SC1091
source "$SIDE_VENV/bin/activate"
pip install --upgrade pip wheel setuptools maturin

# Only what the probe needs — this env is for reading .d files, not for running
# the pipeline.
pip install numpy pandas pyarrow scipy pyyaml

echo "############ Building rustims @ $RUSTIMS_BRANCH into $INSTALL_DIR ############"
bash scripts/setup_imspy.sh

echo "############ Fingerprint ############"
git -C "$INSTALL_DIR" log --format='rustims %h %ad %s' --date=short -1
python3 -c "import importlib.metadata as m; print('imspy_connector', m.version('imspy_connector')); print('imspy_core', m.version('imspy_core'))"
python3 -c "from imspy_core.timstof import TimsDatasetDDA; print('TimsDatasetDDA OK')"

echo "############ Production venv untouched ############"
git -C "$ROOT/rustims" log --format='deployed rustims %h %ad %s' --date=short -1 || true
echo "############ DONE ############"
