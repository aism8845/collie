#!/usr/bin/env bash
set -euo pipefail

# Staged 48h run to bridge a late stiff window:
#   stage 1: baseline solver up to STAGE1_END
#   stage 2: rescue override through STAGE2_END
#   stage 3+: baseline to FINAL_END with automatic rescue re-bridging if needed
#
# Optional:
#   START_RECOVER=<checkpoint_prefix> to skip stage 1.

NP="${NP:-8}"
MPIEXEC="${MPIEXEC:-mpirun}"
EXE="${EXE:-./collie-opt}"
INPUT="${INPUT:-inputs/current/bath_calibration_48h.i}"
BASE_OVR="${BASE_OVR:-inputs/solver_overrides/bath48_best.i}"
RESCUE_OVR="${RESCUE_OVR:-inputs/solver_overrides/bath48_rescue_window.i}"
QUIET_OVR="${QUIET_OVR:-inputs/solver_overrides/quiet_console.i}"

STAGE1_END="${STAGE1_END:-31.25}"
STAGE2_END="${STAGE2_END:-31.35}"
FINAL_END="${FINAL_END:-48.0}"
START_RECOVER="${START_RECOVER:-}"
RESCUE_STEP="${RESCUE_STEP:-0.15}"
MAX_BRIDGES="${MAX_BRIDGES:-12}"
FINAL_EPS="${FINAL_EPS:-1e-6}"

OUT_PREFIX="${OUT_PREFIX:-outputs/bath48_staged}"

mkdir -p outputs outputs/chk

STAGE1_BASE="${OUT_PREFIX}_stage1"
STAGE2_BASE="${OUT_PREFIX}_rescue"
STAGE3_BASE="${OUT_PREFIX}_final"

CHK1_BASE="outputs/chk/bath48_staged_stage1"
CHK2_BASE="outputs/chk/bath48_staged_stage2"
CHK3_BASE="outputs/chk/bath48_staged_stage3"

float_gt() {
  local a="$1"
  local b="$2"
  python3 - "$a" "$b" <<'PY'
import sys
a = float(sys.argv[1])
b = float(sys.argv[2])
sys.exit(0 if a > b else 1)
PY
}

float_add() {
  local a="$1"
  local b="$2"
  python3 - "$a" "$b" <<'PY'
import sys
a = float(sys.argv[1])
b = float(sys.argv[2])
print(a + b)
PY
}

float_min() {
  local a="$1"
  local b="$2"
  python3 - "$a" "$b" <<'PY'
import sys
a = float(sys.argv[1])
b = float(sys.argv[2])
print(a if a < b else b)
PY
}

reached_target() {
  local t_last="$1"
  local t_target="$2"
  local eps="$3"
  python3 - "$t_last" "$t_target" "$eps" <<'PY'
import sys
t_last = float(sys.argv[1])
t_target = float(sys.argv[2])
eps = float(sys.argv[3])
sys.exit(0 if t_last >= (t_target - eps) else 1)
PY
}

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

check_recover_prefix() {
  local recover_prefix="$1"
  local mesh_entry="${recover_prefix}-mesh.cpa.gz"
  if [[ -f "${mesh_entry}" ]]; then
    return 0
  fi
  if [[ -f "${mesh_entry}/${NP}/header.gz" ]]; then
    return 0
  fi
  if [[ -d "${mesh_entry}" ]]; then
    local splits
    splits="$(find "${mesh_entry}" -maxdepth 1 -mindepth 1 -type d -printf '%f ' 2>/dev/null || true)"
    echo "[recover] ERROR: checkpoint split for NP=${NP} not found at ${mesh_entry}/${NP}/header.gz" >&2
    if [[ -n "${splits}" ]]; then
      echo "[recover] available splits: ${splits}" >&2
    fi
    return 1
  fi
  echo "[recover] ERROR: checkpoint mesh entry not found: ${mesh_entry}" >&2
  return 1
}

run_stage_no_recover() {
  local label="$1"
  local end_time="$2"
  local solver_ovr="$3"
  local out_base="$4"
  local chk_base="$5"

  rm -rf "${chk_base}_cp"
  echo "[${label}] run: 0 -> ${end_time}"
  set +e
  "${MPIEXEC}" -np "${NP}" "${EXE}" -i "${INPUT}" "${solver_ovr}" "${QUIET_OVR}" \
    Executioner/end_time="${end_time}" \
    Outputs/exodus=false \
    Outputs/solver_watch/file_base="${out_base}_solver_watch" \
    Outputs/mesh_watch/file_base="${out_base}_mesh_watch" \
    Outputs/chk/file_base="${chk_base}" \
    | tee "${out_base}.log"
  local rc=$?
  set -e
  return "${rc}"
}

run_stage_with_recover() {
  local label="$1"
  local recover_prefix="$2"
  local end_time="$3"
  local solver_ovr="$4"
  local out_base="$5"
  local chk_base="$6"

  if ! check_recover_prefix "${recover_prefix}"; then
    return 1
  fi

  rm -rf "${chk_base}_cp"
  echo "[${label}] run: recover ${recover_prefix} -> ${end_time}"
  set +e
  "${MPIEXEC}" -np "${NP}" "${EXE}" --recover="${recover_prefix}" -i "${INPUT}" "${solver_ovr}" "${QUIET_OVR}" \
    Executioner/end_time="${end_time}" \
    Outputs/exodus=false \
    Outputs/solver_watch/file_base="${out_base}_solver_watch" \
    Outputs/mesh_watch/file_base="${out_base}_mesh_watch" \
    Outputs/chk/file_base="${chk_base}" \
    | tee "${out_base}.log"
  local rc=$?
  set -e
  return "${rc}"
}

audit_stage() {
  local out_base="$1"
  if [[ ! -f "${out_base}_solver_watch.csv" ]]; then
    echo "[audit] skip: missing ${out_base}_solver_watch.csv"
    return 0
  fi
  python3 scripts/audit_one_solver_watch.py \
    "${out_base}_solver_watch.csv" \
    --mesh-watch "${out_base}_mesh_watch.csv" \
    --json-out "${out_base}_audit.json" \
    || true
}

last_time_from_solver_watch() {
  local csv_path="$1"
  python3 - "${csv_path}" <<'PY'
import csv
import math
import sys

path = sys.argv[1]
last_t = math.nan
try:
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            val = row.get("time", "")
            if val not in ("", None):
                last_t = float(val)
except Exception:
    sys.exit(1)

if math.isnan(last_t):
    sys.exit(1)
print(f"{last_t:.12g}")
PY
}

if [[ -n "${START_RECOVER}" ]]; then
  REC1="${START_RECOVER}"
  echo "[stage1] skipped; using START_RECOVER=${REC1}"
else
  run_stage_no_recover "stage1/base" "${STAGE1_END}" "${BASE_OVR}" "${STAGE1_BASE}" "${CHK1_BASE}"
  audit_stage "${STAGE1_BASE}"

  REC1="$(latest_recover_prefix "${CHK1_BASE}" || true)"
  if [[ -z "${REC1}" ]]; then
    echo "[stage1] ERROR: no checkpoint found under ${CHK1_BASE}_cp" >&2
    exit 1
  fi
fi

REC_CURRENT="${REC1}"
REC_FALLBACK="${REC_CURRENT}"

if float_gt "${STAGE2_END}" "${STAGE1_END}"; then
  run_stage_with_recover "stage2/rescue" "${REC_CURRENT}" "${STAGE2_END}" "${RESCUE_OVR}" "${STAGE2_BASE}" "${CHK2_BASE}" || true
  audit_stage "${STAGE2_BASE}"
  REC2="$(latest_recover_prefix "${CHK2_BASE}" || true)"
  if [[ -n "${REC2}" ]]; then
    REC_CURRENT="${REC2}"
    REC_FALLBACK="${REC_CURRENT}"
  fi
fi

echo "[stage3/base] attempt: recover ${REC_CURRENT} -> ${FINAL_END}"
if run_stage_with_recover "stage3/base" "${REC_CURRENT}" "${FINAL_END}" "${BASE_OVR}" "${STAGE3_BASE}" "${CHK3_BASE}"; then
  BASE_RC=0
else
  BASE_RC=$?
fi
audit_stage "${STAGE3_BASE}"

LAST_BASE_T="$(last_time_from_solver_watch "${STAGE3_BASE}_solver_watch.csv" || true)"
if [[ -z "${LAST_BASE_T}" ]]; then
  LAST_BASE_T="${STAGE2_END}"
fi

REC3="$(latest_recover_prefix "${CHK3_BASE}" || true)"
if [[ -n "${REC3}" ]]; then
  REC_FALLBACK="${REC3}"
fi

if [[ "${BASE_RC}" -eq 0 ]] && reached_target "${LAST_BASE_T}" "${FINAL_END}" "${FINAL_EPS}"; then
  echo "[done] staged run complete (single-bridge path)"
  exit 0
fi

echo "[loop] baseline did not reach FINAL_END (rc=${BASE_RC}, t_last=${LAST_BASE_T}); enabling auto re-bridge loop"

loop=1
while [[ "${loop}" -le "${MAX_BRIDGES}" ]]; do
  REC_BASE="$(latest_recover_prefix "${CHK3_BASE}" || true)"
  if [[ -z "${REC_BASE}" ]]; then
    REC_BASE="${REC_FALLBACK}"
    echo "[loop ${loop}] warning: baseline checkpoint missing, reusing ${REC_BASE}"
  fi

  T_START="$(last_time_from_solver_watch "${STAGE3_BASE}_solver_watch.csv" || true)"
  if [[ -z "${T_START}" ]]; then
    T_START="${LAST_BASE_T}"
  fi
  T_RESCUE_END="$(float_min "$(float_add "${T_START}" "${RESCUE_STEP}")" "${FINAL_END}")"

  LOOP_R_BASE="${OUT_PREFIX}_loop${loop}_rescue"
  LOOP_B_BASE="${OUT_PREFIX}_loop${loop}_base"
  LOOP_CHK_R="outputs/chk/bath48_staged_loop${loop}_rescue"
  LOOP_CHK_B="outputs/chk/bath48_staged_loop${loop}_base"

  echo "[loop ${loop}] rescue window ${T_START} -> ${T_RESCUE_END}"
  run_stage_with_recover "loop${loop}/rescue" "${REC_BASE}" "${T_RESCUE_END}" "${RESCUE_OVR}" "${LOOP_R_BASE}" "${LOOP_CHK_R}" || true
  audit_stage "${LOOP_R_BASE}"

  REC_R="$(latest_recover_prefix "${LOOP_CHK_R}" || true)"
  if [[ -z "${REC_R}" ]]; then
    echo "[loop ${loop}] ERROR: no rescue checkpoint found under ${LOOP_CHK_R}_cp" >&2
    exit 1
  fi
  REC_FALLBACK="${REC_R}"

  echo "[loop ${loop}] baseline retry ${T_RESCUE_END} -> ${FINAL_END}"
  if run_stage_with_recover "loop${loop}/base" "${REC_R}" "${FINAL_END}" "${BASE_OVR}" "${LOOP_B_BASE}" "${LOOP_CHK_B}"; then
    LOOP_BASE_RC=0
  else
    LOOP_BASE_RC=$?
  fi
  audit_stage "${LOOP_B_BASE}"

  LAST_BASE_T="$(last_time_from_solver_watch "${LOOP_B_BASE}_solver_watch.csv" || true)"
  if [[ -z "${LAST_BASE_T}" ]]; then
    LAST_BASE_T="${T_RESCUE_END}"
  fi

  STAGE3_BASE="${LOOP_B_BASE}"
  CHK3_BASE="${LOOP_CHK_B}"
  REC_B="$(latest_recover_prefix "${CHK3_BASE}" || true)"
  if [[ -n "${REC_B}" ]]; then
    REC_FALLBACK="${REC_B}"
  fi

  if [[ "${LOOP_BASE_RC}" -eq 0 ]] && reached_target "${LAST_BASE_T}" "${FINAL_END}" "${FINAL_EPS}"; then
    echo "[done] staged run complete after ${loop} additional bridge loop(s)"
    exit 0
  fi

  echo "[loop ${loop}] baseline retry not finished (rc=${LOOP_BASE_RC}, t_last=${LAST_BASE_T})"
  loop=$((loop + 1))
done

echo "[done] reached MAX_BRIDGES=${MAX_BRIDGES} without hitting FINAL_END=${FINAL_END}; last_t=${LAST_BASE_T}" >&2
exit 2
