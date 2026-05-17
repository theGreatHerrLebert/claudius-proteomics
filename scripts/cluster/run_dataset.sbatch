#!/bin/bash
#SBATCH --account=ki-proanagi
#SBATCH --partition=ki-parallel
#SBATCH --cpus-per-task=32
#SBATCH --mem=120G
#SBATCH --time=12:00:00
#SBATCH --job-name=cp-runner
#SBATCH --output=/lustre/project/ki-proanagi/dateschn/logs/run-%x-%j.log
#
# Run the claudius-proteomics 6-step runner on Mogon NHR.
#
# Usage:
#   sbatch scripts/cluster/run_dataset.sbatch PXD019086 --test-mode --max-files 1
#
# All arguments after the script name are forwarded to runner/run_dataset.py.

set -uo pipefail

ROOT=/lustre/project/ki-proanagi/dateschn
PROJ="$ROOT/claudius-proteomics"

export APPTAINER_CACHEDIR="$ROOT/.apptainer/cache"
export APPTAINER_TMPDIR="$ROOT/.apptainer/tmp"
export PIP_CACHE_DIR="$ROOT/.pip-cache"

module purge
module load lang/Python/3.12.3-GCCcore-13.3.0
module load tools/Apptainer/1.3.4-GCCcore-13.3.0

cd "$PROJ"
# shellcheck disable=SC1091
source .venv/bin/activate

echo "Python: $(which python) $(python --version 2>&1)"
echo "Args:   $*"
echo "========================================================================"

python runner/run_dataset.py "$@" \
    --config config/config.mogon.yaml \
    --output-dir "$ROOT/data" \
    --threads "${SLURM_CPUS_PER_TASK:-16}"
