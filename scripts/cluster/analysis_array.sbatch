#!/bin/bash
#SBATCH --account=ki-proanagi
#SBATCH --partition=ki-smallcpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=08:00:00
#SBATCH --job-name=cp-ianalysis
#SBATCH --output=/lustre/project/ki-proanagi/dateschn/logs/ianalysis-%x-%A_%a.log
#
# Cluster-side independent intensity analysis -- one dataset per array task.
# Reads blobs in place on /lustre; emits a compact JSON per dataset under
# data/analysis/{acc}/intensity_independent.json. Runs only on datasets that
# have completed step 5 (merge); others self-report status "missing_inputs".
#
# Usage:
#   sbatch --array=0-49%8 \
#          --export=ALL,ACCESSION_LIST=$PWD/scripts/cluster/top50_accessions.txt \
#          scripts/cluster/analysis_array.sbatch
set -uo pipefail

ROOT=/lustre/project/ki-proanagi/dateschn
PROJ="$ROOT/claudius-proteomics"
LIST="${ACCESSION_LIST:-$PROJ/scripts/cluster/top50_accessions.txt}"

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

echo "Analysis array task $SLURM_ARRAY_TASK_ID -> $ACC"
echo "Python: $(which python) $(python --version 2>&1)"
echo "========================================================================"

python scripts/cluster/cluster_intensity_analysis.py "$ACC" --data-root "$ROOT/data"
