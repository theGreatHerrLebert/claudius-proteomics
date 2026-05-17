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

# --- Provenance: record exactly how this run was executed (see docs/RUN_LEDGER.md) ---
ACC="${1:-unknown}"
PROV_DIR="$ROOT/data/provenance/$ACC"
mkdir -p "$PROV_DIR"
PROV_FILE="$PROV_DIR/run-${SLURM_JOB_ID:-local}.json"
GIT_COMMIT=$(git -C "$PROJ" rev-parse HEAD 2>/dev/null || echo unknown)
GIT_DIRTY=$(git -C "$PROJ" status --porcelain 2>/dev/null | grep -q . && echo true || echo false)
CONFIG_SHA=$(sha256sum "$PROJ/config/config.mogon.yaml" 2>/dev/null | cut -d' ' -f1 || echo unknown)
SAGE_VERSION=$("$ROOT/engines/sage/sage" --version 2>/dev/null || echo unknown)
cat > "$PROV_FILE" <<JSON
{
  "kind": "runner",
  "accession": "$ACC",
  "runner_args": "$*",
  "slurm": {
    "job_id": "${SLURM_JOB_ID:-}",
    "job_name": "${SLURM_JOB_NAME:-}",
    "partition": "${SLURM_JOB_PARTITION:-}",
    "account": "${SLURM_JOB_ACCOUNT:-}",
    "cpus_per_task": "${SLURM_CPUS_PER_TASK:-}",
    "mem_per_node_mb": "${SLURM_MEM_PER_NODE:-}",
    "nodelist": "${SLURM_JOB_NODELIST:-}"
  },
  "code": { "git_commit": "$GIT_COMMIT", "working_tree_dirty": $GIT_DIRTY },
  "config": { "path": "config/config.mogon.yaml", "sha256": "$CONFIG_SHA" },
  "engines": { "sage": "$SAGE_VERSION" },
  "host": "$(hostname)",
  "started_at_utc": "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
}
JSON
echo "Provenance: $PROV_FILE"

python runner/run_dataset.py "$@" \
    --config config/config.mogon.yaml \
    --output-dir "$ROOT/data" \
    --threads "${SLURM_CPUS_PER_TASK:-16}"
