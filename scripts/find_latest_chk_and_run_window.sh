#!/usr/bin/env bash
set -euo pipefail

NP="${NP:-8}"
MPIEXEC="${MPIEXEC:-mpirun}"
EXE="${EXE:-./collie-opt}"
INPUT="${INPUT:-inputs/current/bath48_window_repro.i}"
SOLVER_OVR="${SOLVER_OVR:-inputs/solver_overrides/bath48_best.i}"
QUIET_OVR="${QUIET_OVR:-inputs/solver_overrides/quiet_console.i}"

# Default hint points to a known checkpoint family close to t~31.5.
CHK_HINT="${CHK_HINT:-outputs/chk/bath48_manual_bridge2}"

OUT_BASE="${OUT_BASE:-outputs/window_repro}"
LOG="${OUT_BASE}.log"

mkdir -p outputs outputs/chk

latest_recover_prefix() {
  local base="$1"
  local cp_dir="${base}"
  if [[ "${cp_dir}" != *_cp ]]; then
    cp_dir="${cp_dir}_cp"
  fi
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

pick_recover_prefix() {
  local rec=""

  if rec="$(latest_recover_prefix "${CHK_HINT}" 2>/dev/null)"; then
    printf '%s\n' "${rec}"
    return 0
  fi

  local cp
  for cp in $(ls -1dt outputs/chk/bath48*_cp 2>/dev/null); do
    if rec="$(latest_recover_prefix "${cp}" 2>/dev/null)"; then
      printf '%s\n' "${rec}"
      return 0
    fi
  done

  return 1
}

RECOVER_PREFIX="$(pick_recover_prefix || true)"
if [[ -z "${RECOVER_PREFIX}" ]]; then
  echo "[window] ERROR: no checkpoint prefix found (CHK_HINT=${CHK_HINT})." >&2
  exit 1
fi

echo "[window] recover prefix: ${RECOVER_PREFIX}"
echo "[window] run: ${INPUT} -> end_time=31.60"

"${MPIEXEC}" -np "${NP}" "${EXE}" --recover="${RECOVER_PREFIX}" -i "${INPUT}" "${SOLVER_OVR}" "${QUIET_OVR}" \
  Outputs/exodus=false \
  Outputs/solver_watch/file_base="${OUT_BASE}_solver_watch" \
  Outputs/mesh_watch/file_base="${OUT_BASE}_mesh_watch" \
  Outputs/chk/file_base="outputs/chk/window_repro" \
  | tee "${LOG}"

echo "[window] log: ${LOG}"
