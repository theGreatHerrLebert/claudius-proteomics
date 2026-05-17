#!/bin/bash
#SBATCH --account=ki-proanagi
#SBATCH --partition=ki-smallcpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=03:00:00
#SBATCH --job-name=cp-build-engines
#SBATCH --output=/lustre/project/ki-proanagi/dateschn/logs/build-engines-%j.log
#
# Build the three search engines on Mogon NHR:
#   1. Sage      - compiled from source (cluster glibc 2.28)
#   2. DIA-NN    - Apptainer container (debian:12, glibc 2.36)
#   3. FragPipe  - Apptainer container (ubuntu:22.04 + JDK 17)
#
# Each build runs independently; a failure in one does not abort the others.
# Submit with:  sbatch scripts/cluster/build_engines.sbatch

set -uo pipefail

ROOT=/lustre/project/ki-proanagi/dateschn
export APPTAINER_CACHEDIR="$ROOT/.apptainer/cache"
export APPTAINER_TMPDIR="$ROOT/.apptainer/tmp"
export CARGO_HOME="$ROOT/.cargo"
# GCC 13 on this cluster emits DWARF5, but the system assembler (binutils 2.30)
# rejects --gdwarf-5. Suppress C/asm debug info so cc-rs builds (e.g. `ring`) work.
export CFLAGS="-g0"
export CXXFLAGS="-g0"
mkdir -p "$APPTAINER_CACHEDIR" "$APPTAINER_TMPDIR" "$ROOT/logs"

module purge
module load lang/Rust/1.91.1-GCCcore-13.3.0
module load tools/Apptainer/1.3.4-GCCcore-13.3.0

sage_status=SKIPPED
diann_status=SKIPPED
fragpipe_status=SKIPPED

echo "############ [1/3] Building Sage from source ############"
if ( cd "$ROOT/engines/sage-src" && cargo build --release -p sage-cli ); then
    cp -f "$ROOT/engines/sage-src/target/release/sage" "$ROOT/engines/sage/sage"
    "$ROOT/engines/sage/sage" --version && sage_status=OK || sage_status=FAILED
else
    sage_status=FAILED
fi
echo "Sage: $sage_status"
echo

echo "############ [2/3] Building DIA-NN container ############"
if apptainer build --fakeroot --force \
        "$ROOT/engines/diann.sif" \
        "$ROOT/claudius-proteomics/containers/diann/Singularity.def"; then
    diann_status=OK
else
    diann_status=FAILED
fi
echo "DIA-NN container: $diann_status"
echo

echo "############ [3/3] Building FragPipe container ############"
if apptainer build --fakeroot --force \
        "$ROOT/engines/fragpipe-base.sif" \
        "$ROOT/claudius-proteomics/containers/fragpipe-base/Singularity.def"; then
    fragpipe_status=OK
else
    fragpipe_status=FAILED
fi
echo "FragPipe container: $fragpipe_status"
echo

echo "############ Verifying ############"
if [ "$diann_status" = OK ]; then
    echo "--- DIA-NN ---"
    apptainer exec --bind /lustre "$ROOT/engines/diann.sif" \
        "$ROOT/engines/diann/diann-2.3.2/diann-linux" 2>&1 | head -5 || true
fi
if [ "$fragpipe_status" = OK ]; then
    echo "--- FragPipe (java) ---"
    apptainer exec --bind /lustre "$ROOT/engines/fragpipe-base.sif" java -version 2>&1 | head -2 || true
fi

echo
echo "############ SUMMARY ############"
echo "Sage:              $sage_status"
echo "DIA-NN container:  $diann_status"
echo "FragPipe container:$fragpipe_status"
[ "$sage_status" = OK ] && [ "$diann_status" = OK ] && [ "$fragpipe_status" = OK ]
