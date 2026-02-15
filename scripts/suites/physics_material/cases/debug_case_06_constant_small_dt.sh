#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../../.." ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)" pwd)"
cd "${ROOT}"

CASE_NAME="06_constant_small_dt"
STAMP="${DEBUG_STAMP:-$(date +%Y%m%d_%H%M%S)}"
OUT_ROOT="${DEBUG_OUT_ROOT:-outputs/debug_suite}"
CASE_DIR="${OUT_ROOT}/${STAMP}_${CASE_NAME}"
NP="${DEBUG_NP:-8}"
END_TIME="${DEBUG_END_TIME:-0.03}"

mkdir -p "${CASE_DIR}"

echo "Running ${CASE_NAME}"
echo "  output: ${CASE_DIR}"

mpiexec -n "${NP}" ./collie-opt \
  -i inputs/current/RZ3_RD_AD_patch.i \
  --error \
  --error-unused \
  Executioner/end_time="${END_TIME}" \
  Executioner/TimeStepper/dt=0.005 \
  Outputs/exodus=false \
  Outputs/perf_graph=false \
  Outputs/mesh_watch/file_base="${CASE_DIR}/mesh_watch" \
  Outputs/solver_watch/file_base="${CASE_DIR}/solver_watch" \
  Debug/show_var_residual_norms=true \
  "Executioner/petsc_options=-snes_monitor -snes_converged_reason -ksp_monitor_short -ksp_converged_reason" \
  2>&1 | tee "${CASE_DIR}/run.log"
