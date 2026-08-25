#!/bin/bash
#SBATCH --job-name=hf_tier1
#SBATCH --partition=ki-smallcpu
#SBATCH --time=12:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=2
#SBATCH --array=0-57%10
#SBATCH --output=/lustre/project/ki-proanagi/dateschn/logs/hf_tier1_%A_%a.log
set -uo pipefail
module load lang/Python/3.12.3-GCCcore-13.3.0
BASE=/lustre/project/ki-proanagi/dateschn
ACC=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$BASE/hf_tier1_cohort.txt")
echo "=== task $SLURM_ARRAY_TASK_ID : $ACC : $(date) ==="
"$BASE/claudius-proteomics/.venv/bin/python" "$BASE/build_hf_tier1.py" \
    --accession "$ACC" \
    --data-root "$BASE/data" \
    --out "$BASE/hf_tier1" \
    --repo-root "$BASE/claudius-proteomics"
rc=$?
echo "=== done $ACC rc=$rc : $(date) ==="
exit $rc
