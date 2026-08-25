#!/bin/bash
#SBATCH --job-name=hf_agg_v2
#SBATCH --partition=ki-smallcpu
#SBATCH --time=06:00:00
#SBATCH --mem=96G
#SBATCH --cpus-per-task=2
#SBATCH --output=/lustre/project/ki-proanagi/dateschn/logs/hf_agg_v2_%j.log
set -uo pipefail
module load lang/Python/3.12.3-GCCcore-13.3.0
BASE=/lustre/project/ki-proanagi/dateschn
echo "=== aggregate v0.2 START $(date) ==="
"$BASE/claudius-proteomics/.venv/bin/python" "$BASE/aggregate_hf_corpus.py" \
    --version v0.2 --tier1-dir hf_tier1_v2 --tier3-dir hf_tier3_v2
echo "=== done rc=$? $(date) ==="
