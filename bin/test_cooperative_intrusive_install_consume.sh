#!/usr/bin/env bash
# Verify that a SymEngine built with the cooperative_intrusive RCP backend:
#   1. installs correctly and is consumable via find_package(SymEngine), and
#   2. cooperates with a foreign runtime's reference counting -- here CPython,
#      through a minimal, hand-written nanobind extension (no litgen).
#
# The cooperative_intrusive backend itself has no nanobind dependency; nanobind
# is pulled in here (pip install) purely as a representative foreign runtime to
# exercise the generic incref/decref hooks end-to-end.  The fully featured
# bindings live in the sister `nbsymengine` project.
set -euo pipefail
set -x

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SOURCE_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
CONSUME_DIR="${SCRIPT_DIR}/cooperative_intrusive_consume"

BUILD_ROOT="${1:-/tmp/symengine-cooperative-intrusive-consume}"
INSTALL_PREFIX="${BUILD_ROOT}/install"
EXT_BUILD_DIR="${BUILD_ROOT}/consumer"

echo "=== cooperative_intrusive install-consume test"
echo "SOURCE_DIR:     ${SOURCE_DIR}"
echo "BUILD_ROOT:     ${BUILD_ROOT}"
echo "INSTALL_PREFIX: ${INSTALL_PREFIX}"

# ── 1. nanobind (representative foreign runtime) ────────────────────
echo "=== Installing nanobind"
pip install nanobind

# ── 2. Configure, build & install SymEngine ─────────────────────────
echo "=== Configuring SymEngine with cooperative_intrusive backend"
rm -rf "${BUILD_ROOT}"
cmake -S "${SOURCE_DIR}" -B "${BUILD_ROOT}/build" \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX="${INSTALL_PREFIX}" \
    -DSYMENGINE_RCP_BACKEND=cooperative_intrusive \
    -DBUILD_TESTS=OFF \
    -DBUILD_BENCHMARKS=OFF

echo "=== Building SymEngine"
cmake --build "${BUILD_ROOT}/build" -j"$(nproc)"

echo "=== Installing SymEngine"
cmake --install "${BUILD_ROOT}/build"

# ── 3. Build the minimal nanobind extension against the install ──────
echo "=== Configuring consumer (minimal nanobind extension)"
cmake -S "${CONSUME_DIR}" -B "${EXT_BUILD_DIR}" \
    -DCMAKE_PREFIX_PATH="${INSTALL_PREFIX}" \
    -DCMAKE_BUILD_TYPE=Release

echo "=== Building consumer"
cmake --build "${EXT_BUILD_DIR}" -j"$(nproc)"

# ── 4. Run the Python smoke test ────────────────────────────────────
# nanobind reports instances still alive at interpreter shutdown on stderr and
# leaves the exit status alone, so the stream is part of the result here.  A
# binding that never gives the immortal constants' wrappers back leaks exactly
# that way, silently: the objects are never deleted, so nothing crashes and
# nothing fails -- there is just a Python object attached to a C++ static for
# the life of the process.
echo "=== Running consumer smoke test"
SMOKE_ERR="${BUILD_ROOT}/smoke.stderr"
if ! PYTHONPATH="${EXT_BUILD_DIR}" python "${CONSUME_DIR}/smoke.py" \
        2>"${SMOKE_ERR}"; then
    cat "${SMOKE_ERR}" >&2
    exit 1
fi
cat "${SMOKE_ERR}" >&2
if grep -q "nanobind: leaked" "${SMOKE_ERR}"; then
    echo "=== FAIL: instances were still attached at interpreter shutdown" >&2
    exit 1
fi

echo "=== SUCCESS: cooperative_intrusive install-consume test passed"
rm -rf "${BUILD_ROOT}"
