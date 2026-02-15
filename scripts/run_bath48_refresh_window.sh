#!/usr/bin/env bash
set -euo pipefail

NP="${NP:-8}"
MPIEXEC="${MPIEXEC:-mpirun}"
EXE="${EXE:-./collie-opt}"
INPUT="${INPUT:-inputs/current/bath_calibration_48h.i}"
SOLVER_OVR="${SOLVER_OVR:-inputs/solver_overrides/bath48_best.i}"
QUIET_OVR="${QUIET_OVR:-inputs/solver_overrides/quiet_console.i}"

mkdir -p outputs outputs/chk

CHK_BASE="outputs/chk/bath48_refresh"
STAGE1_BASE="outputs/bath48_refresh_stage1"
STAGE2_BASE="outputs/bath48_refresh"

rm -rf "${CHK_BASE}_cp"

latest_recover_prefix() {
  local base="$1"
  local cp_dir="${base}_cp"
  if [[ ! -d "${cp_dir}" ]]; then
    return 1
  fi

  local newest=""
  newest="$(ls -1 "${cp_dir}"/*-restart-*.rd 2>/dev/null | sort | tail -n 1 || true)"
  if [[ -z "${newest}" ]]; then
    return 1
  fi
  printf '%s\n' "${newest%-restart-*}"
}

echo "[refresh] stage 1: 0 -> 23.0"
"${MPIEXEC}" -np "${NP}" "${EXE}" -i "${INPUT}" "${SOLVER_OVR}" "${QUIET_OVR}" \
  Executioner/end_time=23.0 \
  Outputs/exodus=false \
  Outputs/solver_watch/file_base=${STAGE1_BASE}_solver_watch \
  Outputs/mesh_watch/file_base=${STAGE1_BASE}_mesh_watch \
  Outputs/chk/file_base="${CHK_BASE}" \
  | tee "${STAGE1_BASE}.log"

RECOVER_PREFIX="$(latest_recover_prefix "${CHK_BASE}" || true)"
if [[ -z "${RECOVER_PREFIX}" ]]; then
  echo "[refresh] ERROR: no checkpoint found under ${CHK_BASE}_cp" >&2
  exit 1
fi

echo "[refresh] stage 2: recover -> 25.5 from ${RECOVER_PREFIX}"
"${MPIEXEC}" -np "${NP}" "${EXE}" --recover="${RECOVER_PREFIX}" -i "${INPUT}" "${SOLVER_OVR}" "${QUIET_OVR}" \
  Executioner/end_time=25.5 \
  Outputs/exodus=false \
  Outputs/solver_watch/file_base=${STAGE2_BASE}_solver_watch \
  Outputs/mesh_watch/file_base=${STAGE2_BASE}_mesh_watch \
  Outputs/chk/file_base="${CHK_BASE}" \
  | tee "${STAGE2_BASE}.log"

python3 scripts/audit_one_solver_watch.py \
  "${STAGE2_BASE}_solver_watch.csv" \
  --mesh-watch "${STAGE2_BASE}_mesh_watch.csv" \
  --json-out "${STAGE2_BASE}_audit.json"
