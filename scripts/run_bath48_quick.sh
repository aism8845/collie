#!/usr/bin/env bash
set -euo pipefail

NP="${NP:-8}"
MPIEXEC="${MPIEXEC:-mpirun}"
EXE="${EXE:-./collie-opt}"
INPUT="${INPUT:-inputs/current/bath_calibration_48h.i}"
SOLVER_OVR="${SOLVER_OVR:-inputs/solver_overrides/bath48_best.i}"
QUIET_OVR="${QUIET_OVR:-inputs/solver_overrides/quiet_console.i}"

mkdir -p outputs outputs/chk

OUT_BASE="outputs/bath48_quick"
LOG="${OUT_BASE}.log"

"${MPIEXEC}" -np "${NP}" "${EXE}" -i "${INPUT}" "${SOLVER_OVR}" "${QUIET_OVR}" \
  Executioner/end_time=2 \
  Outputs/exodus=false \
  Outputs/solver_watch/file_base=${OUT_BASE}_solver_watch \
  Outputs/mesh_watch/file_base=${OUT_BASE}_mesh_watch \
  Outputs/chk/file_base=outputs/chk/bath48_quick \
  | tee "${LOG}"

python3 scripts/audit_one_solver_watch.py \
  "${OUT_BASE}_solver_watch.csv" \
  --mesh-watch "${OUT_BASE}_mesh_watch.csv" \
  --json-out "${OUT_BASE}_audit.json"
