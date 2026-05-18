#!/bin/bash
#SBATCH --account=ki-proanagi
#SBATCH --partition=ki-parallel
#SBATCH --cpus-per-task=128
#SBATCH --mem=240G
#SBATCH --time=2-00:00:00
#SBATCH --job-name=cp-PXD073076
#SBATCH --output=/lustre/project/ki-proanagi/dateschn/logs/proc-%x-%j.log
#
# Process PXD073076 — San José runner, Sage + FragPipe (DIA-NN deferred).
exec /lustre/project/ki-proanagi/dateschn/claudius-proteomics/scripts/cluster/_process_runner.sh PXD073076
