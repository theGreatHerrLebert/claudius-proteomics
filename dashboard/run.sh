#!/bin/bash
# Run both backend and frontend for development

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

# Default store path
STORE_PATH="${1:-$PROJECT_ROOT/data/processed/PXD019086/frac01_full.parquet}"

echo "Starting Precursor Browser Dashboard"
echo "Store: $STORE_PATH"
echo ""

# Start backend
echo "Starting backend on http://localhost:8000..."
cd "$SCRIPT_DIR/backend"
$PROJECT_ROOT/.venv/bin/python main.py --store "$STORE_PATH" --port 8000 &
BACKEND_PID=$!

# Wait for backend to start
sleep 2

# Start frontend
echo "Starting frontend on http://localhost:5173..."
cd "$SCRIPT_DIR/frontend"
npm run dev &
FRONTEND_PID=$!

echo ""
echo "Dashboard running!"
echo "  Frontend: http://localhost:5173"
echo "  Backend:  http://localhost:8000"
echo ""
echo "Press Ctrl+C to stop"

# Handle cleanup
cleanup() {
    echo "Stopping..."
    kill $BACKEND_PID 2>/dev/null
    kill $FRONTEND_PID 2>/dev/null
    exit 0
}
trap cleanup SIGINT SIGTERM

# Wait
wait
