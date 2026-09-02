#!/bin/bash
#SBATCH --job-name=hf_aggregate
#SBATCH --partition=ki-smallcpu
#SBATCH --time=06:00:00
#SBATCH --mem=96G
#SBATCH --cpus-per-task=2
#SBATCH --output=/lustre/project/ki-proanagi/dateschn/logs/hf_aggregate_%j.log
set -uo pipefail
module load lang/Python/3.12.3-GCCcore-13.3.0
BASE=/lustre/project/ki-proanagi/dateschn
echo "=== aggregate START $(date) ==="
"$BASE/claudius-proteomics/.venv/bin/python" "$BASE/aggregate_hf_corpus.py"
rc=$?
echo "=== aggregate done rc=$rc $(date) ==="
exit $rc
