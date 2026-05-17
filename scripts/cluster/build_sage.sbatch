#!/bin/bash
#SBATCH --account=ki-proanagi
#SBATCH --partition=ki-smallcpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --job-name=cp-build-sage
#SBATCH --output=/lustre/project/ki-proanagi/dateschn/logs/build-sage-%j.log
#
# Build Sage from the latest upstream lazear/sage (includes the fragment-matching
# fix). Links against the cluster glibc 2.28.
#
# Submit with:  sbatch scripts/cluster/build_sage.sbatch

set -uo pipefail

ROOT=/lustre/project/ki-proanagi/dateschn
export CARGO_HOME="$ROOT/.cargo"
# GCC 13 emits DWARF5; the system assembler (binutils 2.30) rejects it.
export CFLAGS="-g0"
export CXXFLAGS="-g0"

module purge
module load lang/Rust/1.91.1-GCCcore-13.3.0

SRC="$ROOT/engines/sage-src"
rm -rf "$SRC"
git clone --depth 1 https://github.com/lazear/sage.git "$SRC"
cd "$SRC"
echo "Sage commit: $(git rev-parse HEAD)"

cargo build --release -p sage-cli
cp -f target/release/sage "$ROOT/engines/sage/sage"
echo "############ Sage version ############"
"$ROOT/engines/sage/sage" --version
