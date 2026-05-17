#!/bin/bash
#SBATCH --account=ki-proanagi
#SBATCH --partition=ki-smallcpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH --job-name=cp-fetch-pride
#SBATCH --output=/lustre/project/ki-proanagi/dateschn/logs/fetch-pride-%j.log
#
# Fetch one .d folder from a PRIDE zip bundle via HTTP range requests.
#
# Usage:
#   sbatch scripts/cluster/fetch_pride.sbatch <zip-url> <out-dir> [max-d] [contains]

set -uo pipefail

URL="${1:?zip URL required}"
OUT="${2:?output dir required}"
MAXD="${3:-1}"
CONTAINS="${4:-}"

ROOT=/lustre/project/ki-proanagi/dateschn
PROJ="$ROOT/claudius-proteomics"
export PIP_CACHE_DIR="$ROOT/.pip-cache"

module purge
module load lang/Python/3.12.3-GCCcore-13.3.0

cd "$PROJ"
# shellcheck disable=SC1091
source .venv/bin/activate
python -c "import remotezip" 2>/dev/null || pip install --quiet remotezip

python scripts/cluster/fetch_pride_member.py "$URL" --out "$OUT" --max-d "$MAXD" --contains "$CONTAINS"
