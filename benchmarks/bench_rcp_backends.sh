#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
RESULTS_DIR="$PROJECT_DIR/benchmark-results"
TIMESTAMP="$(date +%Y%m%d_%H%M%S)"
RESULTS_FILE="$RESULTS_DIR/rcp_backends_${TIMESTAMP}.txt"

BENCHMARKS="expand1 add1 add_steal symbench lwbench rcp_throughput"
BUILD_JOBS="$(nproc 2>/dev/null || echo 4)"
CMAKE_EXTRA_ARGS=""

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m'

log_info()  { echo -e "${GREEN}[INFO]${NC} $*"; }
log_warn()  { echo -e "${YELLOW}[WARN]${NC} $*"; }
log_error() { echo -e "${RED}[ERROR]${NC} $*"; }

cleanup_build_dir() {
    local build_dir="$1"
    if [[ -d "$build_dir" ]]; then
        log_info "Cleaning up build directory: $build_dir"
        rm -rf "$build_dir"
    fi
}

build_backend() {
    local backend="$1"
    local build_dir="$2"
    local build_type="$3"
    local thread_safe="${4:-OFF}"
    local extra_cmake_args="${5:-}"

    log_info "Building backend=$backend build_type=$build_type thread_safe=$thread_safe"
    log_info "Build directory: $build_dir"

    cleanup_build_dir "$build_dir"
    mkdir -p "$build_dir"

    local cmake_args=(
        -S "$PROJECT_DIR"
        -B "$build_dir"
        -DCMAKE_BUILD_TYPE="$build_type"
        -DSYMENGINE_RCP_BACKEND="$backend"
        -DBUILD_BENCHMARKS=ON
        -DBUILD_BENCHMARKS_GOOGLE=OFF
        -DWITH_SYMENGINE_THREAD_SAFE="$thread_safe"
        -DBUILD_TESTS=OFF
    )

    if [[ -n "$extra_cmake_args" ]]; then
        cmake_args+=($extra_cmake_args)
    fi

    cmake "${cmake_args[@]}" || {
        log_error "CMake configuration failed for backend=$backend"
        return 1
    }

    cmake --build "$build_dir" --parallel "$BUILD_JOBS" $CMAKE_EXTRA_ARGS 2>&1 || {
        log_error "Build failed for backend=$backend"
        return 1
    }

    log_info "Build succeeded for backend=$backend"
}

run_benchmark() {
    local build_dir="$1"
    local bench_name="$2"
    local backend="$3"
    local label="$4"
    local bench_path="$build_dir/benchmarks/$bench_name"

    if [[ ! -x "$bench_path" ]]; then
        log_warn "Benchmark not found: $bench_path"
        echo "=== $label | $bench_name === NOT FOUND" >> "$RESULTS_FILE"
        return 0
    fi

    log_info "Running: $bench_name ($label)"
    echo "" >> "$RESULTS_FILE"
    echo "=== $label | $bench_name ===" >> "$RESULTS_FILE"
    echo "backend: $backend" >> "$RESULTS_FILE"
    echo "date: $(date -Iseconds)" >> "$RESULTS_FILE"
    echo "---" >> "$RESULTS_FILE"

    local start_time end_time elapsed
    start_time=$(date +%s%N)

    "$bench_path" >> "$RESULTS_FILE" 2>&1 || {
        log_error "Benchmark $bench_name failed for $label"
        echo "ERROR: benchmark failed with exit code $?" >> "$RESULTS_FILE"
    }

    end_time=$(date +%s%N)
    elapsed=$(echo "scale=3; ($end_time - $start_time) / 1000000000" | bc)
    echo "elapsed: ${elapsed}s" >> "$RESULTS_FILE"
    echo "---" >> "$RESULTS_FILE"
}

run_all_benchmarks() {
    local build_dir="$1"
    local backend="$2"
    local label="$3"

    log_info "Running all benchmarks for $label"
    for bench in $BENCHMARKS; do
        run_benchmark "$build_dir" "$bench" "$backend" "$label"
    done
}

main() {
    log_info "RCP Backend Benchmark Suite"
    log_info "==========================="
    log_info "Project: $PROJECT_DIR"
    log_info "Results: $RESULTS_FILE"
    log_info "Benchmarks: $BENCHMARKS"
    echo ""

    mkdir -p "$RESULTS_DIR"
    cat > "$RESULTS_FILE" <<EOF
RCP Backend Benchmark Results
=============================
Generated: $(date)
Hostname: $(hostname)
Kernel: $(uname -r)
CPU: $(grep "model name" /proc/cpuinfo 2>/dev/null | head -1 | cut -d: -f2 | xargs || echo "unknown")
Cores: $(nproc 2>/dev/null || echo "unknown")

EOF

    local build_symengine="$PROJECT_DIR/build-bench-symengine"
    local build_teuchos="$PROJECT_DIR/build-bench-teuchos"
    local build_nb="$PROJECT_DIR/build-bench-cooperative_intrusive"
    local build_nb_ts="$PROJECT_DIR/build-bench-cooperative_intrusive-threadsafe"

    trap 'log_warn "Interrupted. Partial results in $RESULTS_FILE"' INT TERM

    log_info "=== Phase 1: symengine backend (default Release) ==="
    build_backend "symengine" "$build_symengine" "Release"
    run_all_benchmarks "$build_symengine" "symengine" "symengine/Release"

    log_info "=== Phase 2: teuchos backend ==="
    build_backend "teuchos" "$build_teuchos" "Release"
    run_all_benchmarks "$build_teuchos" "teuchos" "teuchos/Release"

    log_info "=== Phase 3: cooperative_intrusive backend ==="
    build_backend "cooperative_intrusive" "$build_nb" "Release"
    run_all_benchmarks "$build_nb" "cooperative_intrusive" "cooperative_intrusive/Release"

    log_info "=== Phase 4: cooperative_intrusive backend (thread-safe) ==="
    build_backend "cooperative_intrusive" "$build_nb_ts" "Release" "ON"
    run_all_benchmarks "$build_nb_ts" "cooperative_intrusive" "cooperative_intrusive/Release/ThreadSafe"

    log_info "=== Results Summary ==="
    echo "" >> "$RESULTS_FILE"
    echo "============================" >> "$RESULTS_FILE"
    echo "Run complete: $(date)" >> "$RESULTS_FILE"

    log_info "All benchmarks complete."
    log_info "Results saved to: $RESULTS_FILE"
    echo ""
    log_info "Quick summary (grep for elapsed times):"
    grep -E "^(===|elapsed:)" "$RESULTS_FILE" | head -40

    local do_cleanup="${RCP_BENCH_KEEP_BUILD:-0}"
    if [[ "$do_cleanup" != "1" ]]; then
        cleanup_build_dir "$build_symengine"
        cleanup_build_dir "$build_teuchos"
        cleanup_build_dir "$build_nb"
        cleanup_build_dir "$build_nb_ts"
    else
        log_info "Build directories kept (RCP_BENCH_KEEP_BUILD=1)"
    fi

    log_info "Done."
}

main "$@"
