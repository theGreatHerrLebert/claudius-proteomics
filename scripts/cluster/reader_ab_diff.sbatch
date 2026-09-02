#!/bin/bash
#SBATCH --account=ki-proanagi
#SBATCH --partition=ki-smallcpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --job-name=cp-reader-ab
#SBATCH --output=/lustre/project/ki-proanagi/dateschn/logs/reader-ab-%j.log
#
# A/B the deployed reader against the rebuilt one on the two calibrated axes.
#
# Controls are chosen from the ModelType survey: one MzCalibration ModelType 1
# file (51 of our 54 raw datasets) and one ModelType 2 (3 of 54 — the case the
# C3/C4-duplicate fix touches). The TimsCalibration is ModelType 2 on all of
# them, so the scan->1/K0 axis is the same question everywhere.
#
# Submit after the side build:
#   sbatch --dependency=afterok:<build_jobid> scripts/cluster/reader_ab_diff.sbatch

set -uo pipefail

ROOT=/lustre/project/ki-proanagi/dateschn
P="$ROOT/claudius-proteomics"
OUT="$ROOT/data/analysis/reader_ab"
OLD_VENV="$P/.venv"
NEW_VENV="$P/.venv-readernew"
PROBE="$P/scripts/cluster/probe_reader_axes.py"
DIFF="$P/scripts/cluster/diff_reader_axes.py"

# accession:label
CONTROLS="PXD050530:modeltype1 PXD051790:modeltype2"

module purge
module load lang/Python/3.12.3-GCCcore-13.3.0
mkdir -p "$OUT"

rc_all=0
for entry in $CONTROLS; do
  ACC="${entry%%:*}"; LABEL="${entry##*:}"
  echo ""
  echo "======================================================================"
  echo "  $ACC ($LABEL)"
  echo "======================================================================"

  D=$(find "$ROOT/data/raw/$ACC" -maxdepth 2 -name '*.d' -type d 2>/dev/null | sort | head -1)
  if [ -z "$D" ]; then
    echo "  SKIP: no .d folder under data/raw/$ACC"; rc_all=1; continue
  fi
  echo "  file: $D"

  for side in old new; do
    if [ "$side" = old ]; then V="$OLD_VENV"; else V="$NEW_VENV"; fi
    if [ ! -x "$V/bin/python" ]; then
      echo "  SKIP $side: no interpreter at $V"; rc_all=1; continue 2
    fi
    echo "  --- probing with $side reader ($V)"
    "$V/bin/python" "$PROBE" "$D" -o "$OUT/${ACC}_${side}.npz" || {
      echo "  FAILED: $side probe on $ACC"; rc_all=1; continue 2; }
  done

  echo "  --- diff"
  "$NEW_VENV/bin/python" "$DIFF" "$OUT/${ACC}_old.npz" "$OUT/${ACC}_new.npz" \
    | tee "$OUT/${ACC}_diff.txt"
  [ "${PIPESTATUS[0]}" -ne 0 ] && rc_all=1
done

echo ""
echo "=== reader A/B complete (rc=$rc_all), artifacts in $OUT ==="
exit $rc_all
