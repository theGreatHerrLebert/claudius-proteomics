#!/bin/bash
# San José Worker Entrypoint
#
# Modes:
#   ./entrypoint.sh PXD019086           # Process single dataset
#   ./entrypoint.sh --daemon            # Poll job queue continuously
#   ./entrypoint.sh --job job.yaml      # Process from local job manifest

set -euo pipefail

SCRIPT_DIR="/app/scripts"

# Color output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m'

log_info() { echo -e "${GREEN}[INFO]${NC} $1"; }
log_warn() { echo -e "${YELLOW}[WARN]${NC} $1"; }
log_error() { echo -e "${RED}[ERROR]${NC} $1"; }

# Validate environment
validate_env() {
    local missing=0

    if [[ ! -d "$FRAGPIPE_PATH" ]]; then
        log_error "FragPipe not found at $FRAGPIPE_PATH"
        log_error "Mount your FragPipe installation: -v /path/to/fragpipe:/opt/fragpipe"
        missing=1
    fi

    if [[ ! -f "$DIANN_PATH" ]]; then
        log_error "DIA-NN not found at $DIANN_PATH"
        log_error "Mount your DIA-NN installation: -v /path/to/diann:/opt/diann"
        missing=1
    fi

    if [[ $missing -eq 1 ]]; then
        exit 1
    fi

    log_info "FragPipe: $FRAGPIPE_PATH"
    log_info "DIA-NN: $DIANN_PATH"
}

# Process single dataset by accession
process_accession() {
    local accession=$1
    log_info "Processing dataset: $accession"

    # Generate job manifest from accession
    python3 "$SCRIPT_DIR/generate_job.py" \
        --accession "$accession" \
        --output /tmp/job.yaml

    # Run the job
    python3 "$SCRIPT_DIR/run_job.py" \
        --job /tmp/job.yaml \
        --fragpipe "$FRAGPIPE_PATH" \
        --diann "$DIANN_PATH" \
        --output "$OUTPUT_DIR"

    # Submit results if upload URL configured
    if [[ -n "${RESULTS_UPLOAD_URL:-}" ]]; then
        python3 "$SCRIPT_DIR/submit_results.py" \
            --job /tmp/job.yaml \
            --results "$OUTPUT_DIR/$accession" \
            --url "$RESULTS_UPLOAD_URL"
    fi

    log_info "Completed: $accession"
}

# Process from local job manifest
process_job_file() {
    local job_file=$1
    log_info "Processing job from: $job_file"

    python3 "$SCRIPT_DIR/run_job.py" \
        --job "$job_file" \
        --fragpipe "$FRAGPIPE_PATH" \
        --diann "$DIANN_PATH" \
        --output "$OUTPUT_DIR"

    log_info "Job completed"
}

# Daemon mode: poll job queue
run_daemon() {
    if [[ -z "${JOB_QUEUE_URL:-}" ]]; then
        log_error "JOB_QUEUE_URL not set for daemon mode"
        exit 1
    fi

    log_info "Starting daemon mode, polling: $JOB_QUEUE_URL"
    log_info "Worker ID: ${WORKER_ID:-$(hostname)}"

    while true; do
        # Fetch next job
        job_file=$(python3 "$SCRIPT_DIR/fetch_job.py" \
            --url "$JOB_QUEUE_URL" \
            --worker "${WORKER_ID:-$(hostname)}" \
            --output /tmp/current_job.yaml 2>/dev/null || echo "")

        if [[ -n "$job_file" && -f "$job_file" ]]; then
            log_info "Received job"
            process_job_file "$job_file"

            # Submit results
            if [[ -n "${RESULTS_UPLOAD_URL:-}" ]]; then
                python3 "$SCRIPT_DIR/submit_results.py" \
                    --job "$job_file" \
                    --results "$OUTPUT_DIR" \
                    --url "$RESULTS_UPLOAD_URL"
            fi
        else
            log_info "No jobs available, sleeping 60s..."
            sleep 60
        fi
    done
}

# Main
main() {
    echo "========================================"
    echo "  San José Worker v1.0"
    echo "========================================"

    validate_env

    case "${1:-}" in
        --daemon)
            run_daemon
            ;;
        --job)
            if [[ -z "${2:-}" ]]; then
                log_error "Usage: --job <path/to/job.yaml>"
                exit 1
            fi
            process_job_file "$2"
            ;;
        --help|-h)
            echo "Usage:"
            echo "  $0 PXD019086        Process single dataset"
            echo "  $0 --job job.yaml   Process from job manifest"
            echo "  $0 --daemon         Poll job queue continuously"
            exit 0
            ;;
        PXD*)
            process_accession "$1"
            ;;
        *)
            log_error "Unknown argument: ${1:-}"
            log_error "Use --help for usage"
            exit 1
            ;;
    esac
}

main "$@"
