#!/bin/bash
#SBATCH --account=ki-proanagi
#SBATCH --partition=ki-smallcpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=1-00:00:00
#SBATCH --job-name=cp-normalize
#SBATCH --output=/lustre/project/ki-proanagi/dateschn/logs/normalize-%j.log
#
# Normalise downloaded raw-data layout: extract leftover archives, flatten
# nested .d wrappers, delete redundant archives. See scripts/normalize_raw.py.

set -uo pipefail
ROOT=/lustre/project/ki-proanagi/dateschn
PROJ="$ROOT/claudius-proteomics"

module purge
module load lang/Python/3.12.3-GCCcore-13.3.0

cd "$PROJ"
# shellcheck disable=SC1091
source .venv/bin/activate

echo "Normalising $ROOT/data/raw — $(date -u +%Y-%m-%dT%H:%M:%SZ)"
python scripts/normalize_raw.py --raw-root "$ROOT/data/raw"
