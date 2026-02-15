#!/usr/bin/env bash
set -euo pipefail

NP="${NP:-8}"
EXE="${EXE:-./collie-opt}"
MPIEXEC="${MPIEXEC:-mpirun}"

BASE_INPUT="${BASE_INPUT:-inputs/current/bath_calibration_48h.i}"
SOFT_OVR="${SOFT_OVR:-inputs/solver_overrides/bath48_rescue_window.i}"
STRICT_OVR="${STRICT_OVR:-inputs/solver_overrides/bath48_best.i}"
QUIET_OVR="${QUIET_OVR:-inputs/solver_overrides/quiet_console.i}"

START_RECOVER="${START_RECOVER:-outputs/chk/bath48_softtol_test_cp/0655}"

BRIDGE_END="${BRIDGE_END:-32.0}"
FINAL_END="${FINAL_END:-48.0}"

BRIDGE_BASE="${BRIDGE_BASE:-outputs/bath48_bridge32_soft}"
STRICT_BASE="${STRICT_BASE:-outputs/bath48_after32_strict}"

BRIDGE_CHK="${BRIDGE_CHK:-outputs/chk/bath48_bridge32_soft}"
STRICT_CHK="${STRICT_CHK:-outputs/chk/bath48_after32_strict}"

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

echo "[bridge] recover=${START_RECOVER} -> ${BRIDGE_END}"
rm -rf "${BRIDGE_CHK}_cp"
"${MPIEXEC}" -np "${NP}" "${EXE}" \
  --recover="${START_RECOVER}" \
  -i "${BASE_INPUT}" "${SOFT_OVR}" "${QUIET_OVR}" \
  Executioner/nl_abs_tol=1e-3 \
  Executioner/end_time="${BRIDGE_END}" \
  Outputs/exodus=false \
  Outputs/solver_watch/file_base="${BRIDGE_BASE}_solver_watch" \
  Outputs/mesh_watch/file_base="${BRIDGE_BASE}_mesh_watch" \
  Outputs/chk/file_base="${BRIDGE_CHK}" \
  | tee "${BRIDGE_BASE}.log"

REC_BRIDGE="$(latest_recover_prefix "${BRIDGE_CHK}" || true)"
if [[ -z "${REC_BRIDGE}" ]]; then
  echo "[bridge] ERROR: no checkpoint under ${BRIDGE_CHK}_cp" >&2
  exit 2
fi
echo "[bridge] latest checkpoint: ${REC_BRIDGE}"

echo "[strict] recover=${REC_BRIDGE} -> ${FINAL_END}"
rm -rf "${STRICT_CHK}_cp"
"${MPIEXEC}" -np "${NP}" "${EXE}" \
  --recover="${REC_BRIDGE}" \
  -i "${BASE_INPUT}" "${STRICT_OVR}" "${QUIET_OVR}" \
  Executioner/end_time="${FINAL_END}" \
  Outputs/exodus=false \
  Outputs/solver_watch/file_base="${STRICT_BASE}_solver_watch" \
  Outputs/mesh_watch/file_base="${STRICT_BASE}_mesh_watch" \
  Outputs/chk/file_base="${STRICT_CHK}" \
  | tee "${STRICT_BASE}.log"

echo "[done] strict continuation finished"
