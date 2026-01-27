#!/usr/bin/env bash
#
# Setup script for imspy development environment
#
# This script clones rustims from the feature branch, builds the Rust backend,
# and installs the Python packages.
#
# Prerequisites:
#   - Python 3.12
#   - Rust toolchain (rustup)
#   - maturin (pip install maturin)
#
# Usage:
#   ./scripts/setup_imspy.sh [--dev]
#
# Options:
#   --dev    Install in editable mode for development
#

set -euo pipefail

# Configuration
RUSTIMS_REPO="https://github.com/theGreatHerrLebert/rustims.git"
RUSTIMS_BRANCH="feature/koina-online"
INSTALL_DIR="${INSTALL_DIR:-$HOME/.local/share/claudius-proteomics/rustims}"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

log_info() {
    echo -e "${GREEN}[INFO]${NC} $1"
}

log_warn() {
    echo -e "${YELLOW}[WARN]${NC} $1"
}

log_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

# Parse arguments
DEV_MODE=false
while [[ $# -gt 0 ]]; do
    case $1 in
        --dev)
            DEV_MODE=true
            shift
            ;;
        *)
            log_error "Unknown option: $1"
            exit 1
            ;;
    esac
done

# Check prerequisites
log_info "Checking prerequisites..."

# Check Python version
PYTHON_VERSION=$(python3 --version 2>&1 | cut -d' ' -f2 | cut -d'.' -f1,2)
if [[ "$PYTHON_VERSION" != "3.12" ]]; then
    log_warn "Python 3.12 recommended, found $PYTHON_VERSION"
fi

# Check Rust
if ! command -v cargo &> /dev/null; then
    log_error "Rust not found. Install from https://rustup.rs/"
    exit 1
fi
log_info "Rust version: $(rustc --version)"

# Check maturin
if ! command -v maturin &> /dev/null; then
    log_info "Installing maturin..."
    pip install maturin
fi
log_info "maturin version: $(maturin --version)"

# Clone or update rustims
if [[ -d "$INSTALL_DIR" ]]; then
    log_info "Updating existing rustims installation..."
    cd "$INSTALL_DIR"
    git fetch origin
    git checkout "$RUSTIMS_BRANCH"
    git pull origin "$RUSTIMS_BRANCH"
else
    log_info "Cloning rustims from $RUSTIMS_BRANCH branch..."
    mkdir -p "$(dirname "$INSTALL_DIR")"
    git clone --branch "$RUSTIMS_BRANCH" "$RUSTIMS_REPO" "$INSTALL_DIR"
    cd "$INSTALL_DIR"
fi

log_info "Building and installing rustims packages..."

# Build imspy_connector (Rust -> Python bindings)
log_info "Building imspy_connector..."
cd "$INSTALL_DIR/imspy_connector"
maturin build --release

# Install the wheel
WHEEL=$(ls target/wheels/*.whl | head -1)
log_info "Installing $WHEEL..."
pip install --force-reinstall "$WHEEL"

# Install imspy
log_info "Installing imspy..."
cd "$INSTALL_DIR/imspy"
if $DEV_MODE; then
    pip install -e .
else
    pip install .
fi

# Install imspy-predictors
log_info "Installing imspy-predictors..."
cd "$INSTALL_DIR/packages/imspy-predictors"
if $DEV_MODE; then
    pip install -e .
else
    pip install .
fi

# Verify installation
log_info "Verifying installation..."
python3 -c "
import imspy
print(f'imspy version: {imspy.__version__}')

from imspy.timstof import TimsDatasetDDA
print('TimsDatasetDDA imported successfully')

import imspy_predictors
print('imspy_predictors imported successfully')

from imspy_predictors.utilities.tokenizer import PeptideTokenizer
print('PeptideTokenizer imported successfully')
"

log_info "============================================"
log_info "imspy setup completed successfully!"
log_info "Installation directory: $INSTALL_DIR"
log_info "============================================"

# Print usage example
echo ""
echo "Example usage:"
echo "  from imspy.timstof import TimsDatasetDDA"
echo "  dataset = TimsDatasetDDA('/path/to/data.d')"
echo ""
echo "  from imspy_predictors.ccs import load_deep_ccs_predictor"
echo "  predictor = load_deep_ccs_predictor()"
