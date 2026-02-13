#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${ROOT}"

export SMOKE_INPUT="${SMOKE_INPUT:-inputs/current/RZ3_RD_AD_patch_bath_smoke.i}"
export SMOKE_JAC_INPUT="${SMOKE_JAC_INPUT:-${SMOKE_INPUT}}"
export SMOKE_END_TIME="${SMOKE_END_TIME:-0.2}"
export SMOKE_JAC_END_TIME="${SMOKE_JAC_END_TIME:-0.01}"
export SMOKE_LOG_ROOT="${SMOKE_LOG_ROOT:-outputs/smoke_debug_bath}"

exec ./scripts/smoke_debug_mpi8.sh "$@"
