#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${ROOT}"

# Canonical quick MPI smoke/debug command for this repo.
# Override ranks with SMOKE_NP if needed.
SMOKE_NP="${SMOKE_NP:-8}"

mpiexec -n "${SMOKE_NP}" ./collie-opt \
  -i inputs/2DN/RZ3_RD_AD_patch.i \
  Executioner/end_time=0.03 \
  Outputs/exodus=false \
  Outputs/perf_graph=false
