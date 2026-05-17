#!/bin/bash
#SBATCH --account=ki-proanagi
#SBATCH --partition=ki-smallcpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=2-00:00:00
#SBATCH --job-name=cp-dl
#SBATCH --output=/lustre/project/ki-proanagi/dateschn/logs/dl-%x-%A_%a.log
#
# Download a full PRIDE dataset to data/raw/{accession} on Mogon NHR.
#
# Driven as a SLURM job array: each task picks its accession (one per line)
# from scripts/cluster/download_accessions.txt by array index.
#
# Usage:
#   sbatch --array=0-22%6 scripts/cluster/download_dataset.sbatch
#
# The %6 caps concurrency at 6 simultaneous downloads — polite to the EBI
# HTTPS mirror and keeps the ki-smallcpu footprint small.

set -uo pipefail

ROOT=/lustre/project/ki-proanagi/dateschn
PROJ="$ROOT/claudius-proteomics"
LIST="$PROJ/scripts/cluster/download_accessions.txt"

export PIP_CACHE_DIR="$ROOT/.pip-cache"

module purge
module load lang/Python/3.12.3-GCCcore-13.3.0

cd "$PROJ"
# shellcheck disable=SC1091
source .venv/bin/activate

ACC=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$LIST" | awk '{print $1}')
if [ -z "${ACC:-}" ]; then
    echo "ERROR: no accession at array index $SLURM_ARRAY_TASK_ID in $LIST" >&2
    exit 1
fi

echo "Array task $SLURM_ARRAY_TASK_ID -> $ACC"
echo "Target:  $ROOT/data/raw/$ACC"
echo "Python:  $(which python) $(python --version 2>&1)"
echo "========================================================================"

# --- Provenance (see docs/RUN_LEDGER.md) ---
PROV_DIR="$ROOT/data/provenance/$ACC"
mkdir -p "$PROV_DIR"
GIT_COMMIT=$(git -C "$PROJ" rev-parse HEAD 2>/dev/null || echo unknown)
GIT_DIRTY=$(git -C "$PROJ" status --porcelain 2>/dev/null | grep -q . && echo true || echo false)
PROV_FILE="$PROV_DIR/download-${SLURM_ARRAY_JOB_ID:-local}_${SLURM_ARRAY_TASK_ID:-0}.json"
cat > "$PROV_FILE" <<JSON
{
  "kind": "download",
  "accession": "$ACC",
  "slurm": {
    "array_job_id": "${SLURM_ARRAY_JOB_ID:-}",
    "array_task_id": "${SLURM_ARRAY_TASK_ID:-}",
    "job_id": "${SLURM_JOB_ID:-}",
    "partition": "${SLURM_JOB_PARTITION:-}",
    "nodelist": "${SLURM_JOB_NODELIST:-}"
  },
  "code": { "git_commit": "$GIT_COMMIT", "working_tree_dirty": $GIT_DIRTY },
  "tool": "scripts/download_pride.py",
  "host": "$(hostname)",
  "started_at_utc": "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
}
JSON
echo "Provenance: $PROV_FILE"

python scripts/download_pride.py \
    --accession "$ACC" \
    --output "$ROOT/data/raw/$ACC" \
    --metadata-dir "$ROOT/data/metadata/$ACC" \
    --retry 3
