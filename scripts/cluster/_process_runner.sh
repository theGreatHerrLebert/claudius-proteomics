#!/bin/bash
# Shared processing logic for the per-dataset SLURM scripts in process/.
# Invoked by a thin #SBATCH wrapper:  _process_runner.sh <ACCESSION>
#
# Runs the San José 6-step runner with Sage + FragPipe (DIA-NN deferred — see
# the search-engine-phasing decision). Writes a provenance record per job.
#
# QC gate: when run as a SLURM array task, after each dataset it checks the
# most recent completions; if they show a systemic failure it halts the rest
# of the array (see scripts/cluster/qc_gate.py).
set -uo pipefail

ACC="${1:?accession required}"
ROOT=/lustre/project/ki-proanagi/dateschn
PROJ="$ROOT/claudius-proteomics"

# QC-gate halt sentinel, scoped to this SLURM array job so a new batch starts
# clean. If a prior task in this array tripped the gate, skip immediately.
HALT_SENTINEL="$ROOT/logs/BATCH_QC_HALT.${SLURM_ARRAY_JOB_ID:-none}"
if [ -n "${SLURM_ARRAY_JOB_ID:-}" ] && [ -f "$HALT_SENTINEL" ]; then
    echo "QC gate: batch ${SLURM_ARRAY_JOB_ID} halted (see $HALT_SENTINEL) — skipping $ACC"
    exit 0
fi

export APPTAINER_CACHEDIR="$ROOT/.apptainer/cache"
export APPTAINER_TMPDIR="$ROOT/.apptainer/tmp"
export PIP_CACHE_DIR="$ROOT/.pip-cache"

module purge
module load lang/Python/3.12.3-GCCcore-13.3.0
module load tools/Apptainer/1.3.4-GCCcore-13.3.0

cd "$PROJ"
# shellcheck disable=SC1091
source .venv/bin/activate

echo "Processing $ACC | Sage + FragPipe | $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "Python: $(which python) $(python --version 2>&1)"
echo "========================================================================"

# --- Provenance: record exactly how this run was executed ---
PROV_DIR="$ROOT/data/provenance/$ACC"
mkdir -p "$PROV_DIR"
PROV_FILE="$PROV_DIR/run-${SLURM_JOB_ID:-local}.json"
GIT_COMMIT=$(git -C "$PROJ" rev-parse HEAD 2>/dev/null || echo unknown)
GIT_DIRTY=$(git -C "$PROJ" status --porcelain 2>/dev/null | grep -q . && echo true || echo false)
CONFIG_SHA=$(sha256sum "$PROJ/config/config.mogon.yaml" 2>/dev/null | cut -d' ' -f1 || echo unknown)
SAGE_VERSION=$("$ROOT/engines/sage/sage" --version 2>/dev/null || echo unknown)
# sha256 of the compiled imspy_connector extension — the version string alone
# does not change between rebuilds, so the .so hash is what pins the backend.
IMSPY_SO_SHA=$(python -c "import imspy_connector, glob, os, hashlib; d=os.path.dirname(imspy_connector.__file__); fs=sorted(glob.glob(os.path.join(d,'*.so'))); print(hashlib.sha256(open(fs[0],'rb').read()).hexdigest() if fs else 'unknown')" 2>/dev/null || echo unknown)
cat > "$PROV_FILE" <<JSON
{
  "kind": "runner",
  "accession": "$ACC",
  "engines": ["sage", "fragpipe"],
  "slurm": {
    "job_id": "${SLURM_JOB_ID:-}",
    "job_name": "${SLURM_JOB_NAME:-}",
    "partition": "${SLURM_JOB_PARTITION:-}",
    "cpus_per_task": "${SLURM_CPUS_PER_TASK:-}",
    "mem_per_node_mb": "${SLURM_MEM_PER_NODE:-}",
    "nodelist": "${SLURM_JOB_NODELIST:-}"
  },
  "code": { "git_commit": "$GIT_COMMIT", "working_tree_dirty": $GIT_DIRTY },
  "config": { "path": "config/config.mogon.yaml", "sha256": "$CONFIG_SHA" },
  "engines_versions": { "sage": "$SAGE_VERSION", "imspy_connector_so_sha256": "$IMSPY_SO_SHA" },
  "host": "$(hostname)",
  "started_at_utc": "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
}
JSON
echo "Provenance: $PROV_FILE"

python runner/run_dataset.py "$ACC" \
    --config config/config.mogon.yaml \
    --engines sage fragpipe \
    --output-dir "$ROOT/data" \
    --threads "${SLURM_CPUS_PER_TASK:-16}" \
    --local-data "$ROOT/data/raw/$ACC" \
    --resume
RUN_RC=$?

# --- QC gate: halt the array if recent completions show a systemic failure ---
if [ -n "${SLURM_ARRAY_JOB_ID:-}" ]; then
    if ! python "$PROJ/scripts/cluster/qc_gate.py" --batch-check --data-dir "$ROOT/data"; then
        echo "QC gate: halting batch ${SLURM_ARRAY_JOB_ID} — recent datasets failed QC"
        touch "$HALT_SENTINEL"
        scancel --state=PENDING "$SLURM_ARRAY_JOB_ID" 2>/dev/null || true
    fi
fi

exit "$RUN_RC"
