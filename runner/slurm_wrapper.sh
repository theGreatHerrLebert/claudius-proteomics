#!/bin/bash
#SBATCH --job-name=sanjose_%j
#SBATCH --partition=parallel
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=logs/sanjose_%j_%x.out
#SBATCH --error=logs/sanjose_%j_%x.err
#SBATCH --mail-type=END,FAIL

# San José Runner - SLURM Wrapper Script
#
# Usage:
#   sbatch runner/slurm_wrapper.sh PXD019086
#   sbatch runner/slurm_wrapper.sh PXD019086 --test-mode
#   sbatch runner/slurm_wrapper.sh PXD019086 --steps 1 2 3
#   sbatch runner/slurm_wrapper.sh PXD019086 --resume
#
# Resource overrides (examples):
#   sbatch --cpus-per-task=32 --mem=128G runner/slurm_wrapper.sh PXD019086
#   sbatch --partition=gpu --gres=gpu:1 runner/slurm_wrapper.sh PXD019086

set -euo pipefail

# ============================================================================
# Configuration
# ============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"

# Paths (adjust for your HPC environment)
CONDA_ENV="${CONDA_ENV:-sanjose}"
CONFIG_FILE="${CONFIG_FILE:-${PROJECT_DIR}/config/config.yaml}"
OUTPUT_DIR="${OUTPUT_DIR:-${PROJECT_DIR}/data}"
LOG_DIR="${PROJECT_DIR}/logs"

# ============================================================================
# Parse arguments
# ============================================================================

if [[ $# -lt 1 ]]; then
    echo "Usage: sbatch runner/slurm_wrapper.sh ACCESSION [OPTIONS]"
    echo ""
    echo "Arguments:"
    echo "  ACCESSION     PRIDE accession (e.g., PXD019086)"
    echo ""
    echo "Options (passed to run_dataset.py):"
    echo "  --test-mode   Process limited files for testing"
    echo "  --max-files N Maximum files in test mode (default: 3)"
    echo "  --resume      Resume from checkpoint if exists"
    echo "  --steps N...  Specific steps to run (1-5)"
    echo "  --threads N   Number of threads (default: \$SLURM_CPUS_PER_TASK)"
    echo ""
    echo "Environment variables:"
    echo "  CONDA_ENV     Conda environment name (default: sanjose)"
    echo "  CONFIG_FILE   Config file path (default: config/config.yaml)"
    echo "  OUTPUT_DIR    Output base directory (default: data)"
    exit 1
fi

ACCESSION="$1"
shift
EXTRA_ARGS=("$@")

# ============================================================================
# Environment setup
# ============================================================================

echo "========================================================================"
echo "  San José Runner - SLURM Job"
echo "========================================================================"
echo ""
echo "Job Information:"
echo "  Job ID:       ${SLURM_JOB_ID:-local}"
echo "  Job Name:     ${SLURM_JOB_NAME:-local}"
echo "  Node:         ${SLURM_NODELIST:-$(hostname)}"
echo "  CPUs:         ${SLURM_CPUS_PER_TASK:-$(nproc)}"
echo "  Memory:       ${SLURM_MEM_PER_NODE:-unknown} MB"
echo "  Partition:    ${SLURM_JOB_PARTITION:-local}"
echo ""
echo "Run Parameters:"
echo "  Accession:    ${ACCESSION}"
echo "  Config:       ${CONFIG_FILE}"
echo "  Output:       ${OUTPUT_DIR}"
echo "  Extra args:   ${EXTRA_ARGS[*]:-none}"
echo ""
echo "Start time:     $(date)"
echo "========================================================================"

# Create log directory
mkdir -p "$LOG_DIR"

# Change to project directory
cd "$PROJECT_DIR"

# Load modules (adjust for your HPC environment)
if command -v module &> /dev/null; then
    echo "Loading modules..."
    module purge 2>/dev/null || true
    module load lang/Python/3.12 2>/dev/null || true
    module load bio/FragPipe 2>/dev/null || true
    module load bio/DIA-NN 2>/dev/null || true
fi

# Activate conda environment if available
if command -v conda &> /dev/null; then
    echo "Activating conda environment: ${CONDA_ENV}"
    # shellcheck disable=SC1091
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate "$CONDA_ENV" 2>/dev/null || {
        echo "Warning: Could not activate conda env '${CONDA_ENV}', using current environment"
    }
elif [[ -f "${PROJECT_DIR}/.venv/bin/activate" ]]; then
    echo "Activating virtual environment..."
    # shellcheck disable=SC1091
    source "${PROJECT_DIR}/.venv/bin/activate"
fi

# Print Python info
echo ""
echo "Python environment:"
which python
python --version
echo ""

# ============================================================================
# Run the pipeline
# ============================================================================

# Set number of threads from SLURM allocation
THREADS="${SLURM_CPUS_PER_TASK:-16}"

# Build command
CMD=(
    python
    "${PROJECT_DIR}/runner/run_dataset.py"
    "$ACCESSION"
    --config "$CONFIG_FILE"
    --output-dir "$OUTPUT_DIR"
    --threads "$THREADS"
)

# Add extra arguments
if [[ ${#EXTRA_ARGS[@]} -gt 0 ]]; then
    CMD+=("${EXTRA_ARGS[@]}")
fi

echo "Executing command:"
echo "  ${CMD[*]}"
echo ""
echo "========================================================================"

# Run with timing
START_TIME=$(date +%s)

"${CMD[@]}"
EXIT_CODE=$?

END_TIME=$(date +%s)
DURATION=$((END_TIME - START_TIME))
HOURS=$((DURATION / 3600))
MINUTES=$(((DURATION % 3600) / 60))
SECONDS=$((DURATION % 60))

echo ""
echo "========================================================================"
echo "  Job Complete"
echo "========================================================================"
echo "Exit code:      ${EXIT_CODE}"
echo "Duration:       ${HOURS}h ${MINUTES}m ${SECONDS}s"
echo "End time:       $(date)"
echo "========================================================================"

exit $EXIT_CODE
