#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${ROOT}"

NP="${NP:-8}"
MPI_LAUNCHER="${MPI_LAUNCHER:-}"
if [[ -z "${MPI_LAUNCHER}" ]]; then
  MPI_LAUNCHER="$(command -v mpiexec || command -v mpirun || true)"
fi
if [[ -z "${MPI_LAUNCHER}" ]]; then
  echo "ERROR: no MPI launcher found (mpiexec/mpirun)." >&2
  exit 1
fi

STAMP="$(date +%Y%m%d_%H%M%S)"
OUT_DIR="${OUT_DIR:-outputs/full24_bath/${STAMP}}"
mkdir -p "${OUT_DIR}"

LOG_FILE="${OUT_DIR}/run.log"

"${MPI_LAUNCHER}" -n "${NP}" ./collie-opt \
  -i inputs/current/RZ3_RD_AD_patch_bath_24h.i \
  --error \
  --error-unused \
  Outputs/file_base="${OUT_DIR}/bath_24h" \
  Outputs/mesh_watch/file_base="${OUT_DIR}/mesh_watch" \
  Outputs/solver_watch/file_base="${OUT_DIR}/solver_watch" \
  Debug/show_var_residual_norms=true \
  Executioner/petsc_options='-snes_monitor -snes_converged_reason -ksp_monitor_short -ksp_converged_reason' \
  2>&1 | tee "${LOG_FILE}"

cat <<MSG
Full 24h run complete.
  Output dir: ${OUT_DIR}
  Log file:   ${LOG_FILE}
MSG
