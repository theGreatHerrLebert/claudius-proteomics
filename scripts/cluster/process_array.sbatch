#!/bin/bash
#SBATCH --account=ki-proanagi
#SBATCH --partition=ki-parallel
#SBATCH --cpus-per-task=128
#SBATCH --mem=240G
#SBATCH --time=2-00:00:00
#SBATCH --job-name=cp-proc
#SBATCH --output=/lustre/project/ki-proanagi/dateschn/logs/proc-%x-%A_%a.log
#
# Process a batch of datasets through the San Jose runner (Sage + FragPipe,
# DIA-NN deferred). SLURM array: each task picks its accession from a list
# by array index, then reuses the shared _process_runner.sh logic.
#
# Usage (cap concurrency with %N so the cluster is not flooded):
#   sbatch --array=0-49%4 \
#          --export=ALL,ACCESSION_LIST=$PWD/scripts/cluster/top50_accessions.txt \
#          scripts/cluster/process_array.sbatch
set -uo pipefail

ROOT=/lustre/project/ki-proanagi/dateschn
PROJ="$ROOT/claudius-proteomics"
LIST="${ACCESSION_LIST:-$PROJ/scripts/cluster/top50_accessions.txt}"

ACC=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$LIST" | awk '{print $1}')
if [ -z "${ACC:-}" ]; then
    echo "ERROR: no accession at array index $SLURM_ARRAY_TASK_ID in $LIST" >&2
    exit 1
fi

echo "Process array task $SLURM_ARRAY_TASK_ID -> $ACC"
exec "$PROJ/scripts/cluster/_process_runner.sh" "$ACC"
