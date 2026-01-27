#!/usr/bin/env bash
#
# Run FragPipe on timsTOF data
#
# FragPipe must be provided by user due to licensing requirements.
# This script wraps the FragPipe headless execution.
#

set -euo pipefail

# Parse arguments
FRAGPIPE_PATH=""
INPUT_DIR=""
OUTPUT_DIR=""
FASTA=""
WORKFLOW="Default"
THREADS=16

while [[ $# -gt 0 ]]; do
    case $1 in
        --fragpipe)
            FRAGPIPE_PATH="$2"
            shift 2
            ;;
        --input)
            INPUT_DIR="$2"
            shift 2
            ;;
        --output)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        --fasta)
            FASTA="$2"
            shift 2
            ;;
        --workflow)
            WORKFLOW="$2"
            shift 2
            ;;
        --threads)
            THREADS="$2"
            shift 2
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

# Validate required arguments
if [[ -z "$FRAGPIPE_PATH" ]]; then
    echo "ERROR: --fragpipe path is required"
    exit 1
fi

if [[ -z "$INPUT_DIR" ]]; then
    echo "ERROR: --input directory is required"
    exit 1
fi

if [[ -z "$OUTPUT_DIR" ]]; then
    echo "ERROR: --output directory is required"
    exit 1
fi

if [[ -z "$FASTA" ]]; then
    echo "ERROR: --fasta database is required"
    exit 1
fi

# Check FragPipe exists
if [[ ! -d "$FRAGPIPE_PATH" ]]; then
    echo "ERROR: FragPipe not found at $FRAGPIPE_PATH"
    echo "Please download FragPipe from https://fragpipe.nesvilab.org/"
    echo "and update the path in config/config.yaml"
    exit 1
fi

FRAGPIPE_BIN="$FRAGPIPE_PATH/bin/fragpipe"
if [[ ! -f "$FRAGPIPE_BIN" ]]; then
    echo "ERROR: FragPipe binary not found at $FRAGPIPE_BIN"
    exit 1
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Find all .d folders (timsTOF raw data)
echo "Finding raw files in $INPUT_DIR..."
RAW_FILES=$(find "$INPUT_DIR" -maxdepth 2 -name "*.d" -type d | tr '\n' ' ')

if [[ -z "$RAW_FILES" ]]; then
    echo "ERROR: No .d folders found in $INPUT_DIR"
    exit 1
fi

echo "Found raw files: $RAW_FILES"

# Generate manifest file for FragPipe
MANIFEST="$OUTPUT_DIR/fragpipe_manifest.fp-manifest"
echo "Generating manifest: $MANIFEST"

> "$MANIFEST"
for raw in $RAW_FILES; do
    # Format: filepath \t experiment \t bioreplicate \t techreplicate \t datatype
    basename=$(basename "$raw" .d)
    echo -e "$raw\t$basename\t1\t1\tDDA" >> "$MANIFEST"
done

# Run FragPipe headless
echo "Running FragPipe..."
echo "  Workflow: $WORKFLOW"
echo "  Threads: $THREADS"
echo "  Output: $OUTPUT_DIR"

"$FRAGPIPE_BIN" \
    --headless \
    --workflow "$WORKFLOW" \
    --manifest "$MANIFEST" \
    --workdir "$OUTPUT_DIR" \
    --db-path "$FASTA" \
    --threads "$THREADS"

echo "FragPipe completed successfully"
