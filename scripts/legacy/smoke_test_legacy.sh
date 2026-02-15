#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "${ROOT}"

# Stashed legacy smoke entrypoint (2DN-era input).
export SMOKE_REQUIRE_CURRENT_INPUT=0
export SMOKE_INPUT="${SMOKE_INPUT:-inputs/legacy/2DN/RZ_smoke.i}"
export SMOKE_JAC_INPUT="${SMOKE_JAC_INPUT:-${SMOKE_INPUT}}"
export SMOKE_LOG_ROOT="${SMOKE_LOG_ROOT:-outputs/smoke_debug_legacy}"

exec ./scripts/smoke_debug_mpi8.sh "$@"
