#!/bin/bash
#SBATCH --account=ki-proanagi
#SBATCH --partition=ki-parallel
#SBATCH --cpus-per-task=64
#SBATCH --mem=200G
#SBATCH --time=06:00:00
#SBATCH --job-name=fp-retest
#SBATCH --output=/lustre/project/ki-proanagi/dateschn/logs/fp-retest-%j.log
#
# Standalone FragPipe re-test on PXD046777 with the concurrency fix:
# capped threads (64) + stable per-dataset temp dir. Validates the fix
# without touching the runner state of telemetry job 259116.
set -uo pipefail
ROOT=/lustre/project/ki-proanagi/dateschn
PROJ="$ROOT/claudius-proteomics"
OUT="$ROOT/data/fragpipe_retest_PXD046777"

module purge
module load lang/Python/3.12.3-GCCcore-13.3.0
module load tools/Apptainer/1.3.4-GCCcore-13.3.0
module load lang/Java/17.0.6

cd "$PROJ"
source .venv/bin/activate
mkdir -p "$OUT/tmp"

echo "FragPipe re-test PXD046777 — $(date -u +%Y-%m-%dT%H:%M:%SZ)"
python scripts/run_fragpipe.py \
    --fragpipe "$ROOT/engines/fragpipe-run" \
    --input "$ROOT/data/raw/PXD046777" \
    --output "$OUT" \
    --fasta "$PROJ/resources/fasta/search_db/human_decoys.fasta" \
    --threads 64 \
    --temp-dir "$OUT/tmp"
echo "exit=$? — $(date -u +%Y-%m-%dT%H:%M:%SZ)"
