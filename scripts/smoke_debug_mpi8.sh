#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${ROOT}"

# Canonical smoke + debug harness:
#   1) strict input check (MOOSE/front-end)
#   2) MPI solver/physics run with PETSc monitors
#   3) optional serial Jacobian spot-check
#
# Default MPI ranks remain 8. Override with SMOKE_NP.
SMOKE_NP="${SMOKE_NP:-8}"
SMOKE_INPUT="${SMOKE_INPUT:-inputs/current/RZ3_RD_AD_patch.i}"
SMOKE_END_TIME="${SMOKE_END_TIME:-0.03}"
SMOKE_REQUIRE_CURRENT_INPUT="${SMOKE_REQUIRE_CURRENT_INPUT:-1}"

SMOKE_RUN_JACOBIAN="${SMOKE_RUN_JACOBIAN:-1}"
SMOKE_JAC_INPUT="${SMOKE_JAC_INPUT:-${SMOKE_INPUT}}"
SMOKE_JAC_END_TIME="${SMOKE_JAC_END_TIME:-0.01}"
# Jacobian fail threshold defaults on for regression catching.
SMOKE_JAC_FAIL_THRESHOLD="${SMOKE_JAC_FAIL_THRESHOLD:-1e-3}"
SMOKE_JAC_AUTO_MESH_REDUCE="${SMOKE_JAC_AUTO_MESH_REDUCE:-1}"
SMOKE_JAC_NX="${SMOKE_JAC_NX:-40}"
SMOKE_JAC_NY="${SMOKE_JAC_NY:-12}"
SMOKE_JAC_OVERRIDES="${SMOKE_JAC_OVERRIDES:-}"
SMOKE_JAC_VIEW="${SMOKE_JAC_VIEW:-0}"

# Slowdown hotspot thresholds for report generation.
SMOKE_SLOW_NL_ITS="${SMOKE_SLOW_NL_ITS:-6}"
SMOKE_SLOW_LIN_ITS="${SMOKE_SLOW_LIN_ITS:-40}"
SMOKE_DT_DROP_FACTOR="${SMOKE_DT_DROP_FACTOR:-0.67}"

# Physics sanity thresholds (warning by default; fail when enabled).
SMOKE_FAIL_ON_PHYSICS="${SMOKE_FAIL_ON_PHYSICS:-0}"
SMOKE_FAIL_ON_NUMERICS="${SMOKE_FAIL_ON_NUMERICS:-0}"
SMOKE_N_MIN_FLOOR="${SMOKE_N_MIN_FLOOR:--1e-8}"
SMOKE_J_MIN_FLOOR="${SMOKE_J_MIN_FLOOR:-0.0}"
SMOKE_DPHYS_MIN_FLOOR="${SMOKE_DPHYS_MIN_FLOOR:-1e-12}"
SMOKE_GATE_MIN_FLOOR="${SMOKE_GATE_MIN_FLOOR:-0.0}"
SMOKE_GATE_MAX_CEIL="${SMOKE_GATE_MAX_CEIL:-1.000001}"

SMOKE_LOG_ROOT="${SMOKE_LOG_ROOT:-outputs/smoke_debug}"
STAMP="$(date +%Y%m%d_%H%M%S)"
RUN_DIR="${SMOKE_LOG_ROOT}/${STAMP}"
mkdir -p "${RUN_DIR}"

run_and_log() {
  local name="$1"
  shift
  local log="${RUN_DIR}/${name}.log"

  echo
  echo "==> ${name}"
  echo "+ $*"
  "$@" 2>&1 | tee "${log}"
}

run_and_log_allow_fail() {
  local name="$1"
  shift
  local log="${RUN_DIR}/${name}.log"
  local ec

  echo
  echo "==> ${name}"
  echo "+ $*"

  set +e
  "$@" 2>&1 | tee "${log}"
  ec=$?
  set -e

  return "${ec}"
}

is_mpi_sandbox_error() {
  local log="$1"
  grep -Eq \
    'xnet_pep_sock_create\(\).*Operation not permitted|shm_open error .*Permission denied|OFI endpoint open failed|HYDU_sock_listen .*Operation not permitted' \
    "${log}"
}

is_not_enough_slots_error() {
  local log="$1"
  grep -Eq 'not enough slots available|There are not enough slots available' "${log}"
}

looks_like_legacy_input() {
  local input_path="$1"
  [[ "${input_path}" == inputs/2DN/* ]] || \
    [[ "${input_path}" == inputs/legacy/* ]] || \
    [[ "${input_path}" == */inputs/2DN/* ]] || \
    [[ "${input_path}" == */inputs/legacy/* ]]
}

require_smoke_input() {
  local input_path="$1"
  local label="$2"

  if [[ ! -f "${input_path}" ]]; then
    echo "ERROR: ${label} file does not exist: ${input_path}" >&2
    exit 1
  fi

  if [[ "${SMOKE_REQUIRE_CURRENT_INPUT}" == "1" ]] && looks_like_legacy_input "${input_path}"; then
    echo "ERROR: ${label} is pointing to legacy inputs: ${input_path}" >&2
    echo "Set SMOKE_REQUIRE_CURRENT_INPUT=0 to bypass this guard intentionally." >&2
    exit 1
  fi
}

summarize_jacobian() {
  local log="$1"
  local ratio

  ratio="$(grep -Eo '\|\|J - Jfd\|\|_F/\|\|J\|\|_F = [0-9.eE+-]+' "${log}" | tail -n 1 | awk '{print $NF}' || true)"
  if [[ -z "${ratio}" ]]; then
    echo "Jacobian summary: no ratio line found in ${log}"
    return
  fi

  echo "Jacobian summary: ||J - Jfd||_F/||J||_F = ${ratio}"
  if [[ -n "${SMOKE_JAC_FAIL_THRESHOLD}" ]]; then
    if awk -v ratio="${ratio}" -v thresh="${SMOKE_JAC_FAIL_THRESHOLD}" 'BEGIN { exit !(ratio > thresh) }'; then
      echo "ERROR: Jacobian ratio ${ratio} exceeds threshold ${SMOKE_JAC_FAIL_THRESHOLD}" >&2
      exit 1
    fi
  fi
}

copy_watch_csvs() {
  local copied=0
  local src
  local dst

  for src in outputs/solver_watch.csv outputs/mesh_watch.csv; do
    if [[ -f "${src}" ]]; then
      dst="${RUN_DIR}/$(basename "${src}")"
      cp -f "${src}" "${dst}"
      copied=1
    fi
  done

  [[ "${copied}" -eq 1 ]]
}

build_solver_hotspots_csv() {
  local solver_csv="$1"
  local hotspots_csv="$2"

  awk -F, '
    function has(name) { return (name in idx) }
    function missing(raw, low) {
      low = tolower(raw)
      return (raw == "" || low == "na" || low == "nan" || low == "+inf" || low == "-inf" || low == "inf")
    }
    function sval(name, fallback, c) {
      c = idx[name]
      if (c)
        return $(c)
      return fallback
    }
    function nval(name, fallback, c, raw) {
      c = idx[name]
      if (!c)
        return fallback
      raw = $(c)
      if (missing(raw))
        return fallback
      return raw + 0
    }
    NR == 1 {
      for (i = 1; i <= NF; i++)
        idx[$i] = i

      if (!has("time") || !has("dt") || !has("nonlinear_its") || !has("linear_its")) {
        print "ERROR: solver_watch.csv missing one of required columns: time, dt, nonlinear_its, linear_its" > "/dev/stderr"
        exit 2
      }

      print "time,dt,dt_ratio,nonlinear_its,linear_its,mesh_distortion_warning,J_min_1,metric2_min_1,n_min,n_max,n_span,avg_Dphys,min_Dphys,avg_Dref,min_Dref,D_ref_phys_ratio,avg_gate_tot,min_gate_tot,max_gate_tot,avg_gp,avg_press,avg_fa,avg_ke,avg_kh,avg_eta,max_phi_cell"
      next
    }
    {
      time = nval("time", 0)
      dt = nval("dt", 0)
      nl = nval("nonlinear_its", 0)
      li = nval("linear_its", 0)
      dt_ratio = (prev_dt > 0) ? dt / prev_dt : 1.0

      printf "%.12g,%.12g,%.12g,%g,%g,%g,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n",
        time, dt, dt_ratio, nl, li,
        nval("mesh_distortion_warning", 0),
        sval("J_min_1", "NA"),
        sval("metric2_min_1", "NA"),
        sval("n_min", "NA"),
        sval("n_max", "NA"),
        sval("n_span", "NA"),
        sval("avg_Dphys", "NA"),
        sval("min_Dphys", "NA"),
        sval("avg_Dref", "NA"),
        sval("min_Dref", "NA"),
        sval("D_ref_phys_ratio", "NA"),
        sval("avg_gate_tot", "NA"),
        sval("min_gate_tot", "NA"),
        sval("max_gate_tot", "NA"),
        sval("avg_gp", "NA"),
        sval("avg_press", "NA"),
        sval("avg_fa", "NA"),
        sval("avg_ke", "NA"),
        sval("avg_kh", "NA"),
        sval("avg_eta", "NA"),
        sval("max_phi_cell", "NA")

      prev_dt = dt
    }
  ' "${solver_csv}" > "${hotspots_csv}"
}

write_slowdown_report() {
  local hotspots_csv="$1"
  local report_txt="$2"

  {
    echo "Slowdown hotspot criteria:"
    echo "  nonlinear_its >= ${SMOKE_SLOW_NL_ITS}"
    echo "  linear_its >= ${SMOKE_SLOW_LIN_ITS}"
    echo "  mesh_distortion_warning > 0"
    echo "  dt_ratio < ${SMOKE_DT_DROP_FACTOR}"
    echo
    echo "Hotspot rows (sorted by nonlinear_its then linear_its):"
    head -n 1 "${hotspots_csv}"
    awk -F, -v nl="${SMOKE_SLOW_NL_ITS}" -v lin="${SMOKE_SLOW_LIN_ITS}" -v dt_drop="${SMOKE_DT_DROP_FACTOR}" '
      NR == 1 { next }
      ($4 + 0 >= nl) || ($5 + 0 >= lin) || ($6 + 0 > 0) || ($3 + 0 < dt_drop) { print }
    ' "${hotspots_csv}" | sort -t, -k4,4nr -k5,5nr | head -n 16
    echo
    echo "Run extrema:"
    awk -F, '
      function missing(raw, low) {
        low = tolower(raw)
        return (raw == "" || low == "na" || low == "nan" || low == "+inf" || low == "-inf" || low == "inf")
      }
      NR == 1 {
        for (i = 1; i <= NF; i++)
          idx[$i] = i
        next
      }
      {
        t = $(idx["time"]) + 0
        dt = $(idx["dt"]) + 0
        nl = $(idx["nonlinear_its"]) + 0
        li = $(idx["linear_its"]) + 0
        dtr = $(idx["dt_ratio"]) + 0

        if (!seen || dt < min_dt) { min_dt = dt; min_dt_t = t }
        if (!seen || dtr < min_dtr) { min_dtr = dtr; min_dtr_t = t }
        if (!seen || nl > max_nl) { max_nl = nl; max_nl_t = t }
        if (!seen || li > max_li) { max_li = li; max_li_t = t }
        if ("J_min_1" in idx) {
          raw = $(idx["J_min_1"])
          if (!missing(raw)) {
            j = raw + 0
            if (!seen_j || j < min_j) { min_j = j; min_j_t = t }
            seen_j = 1
          }
        }
        if ("n_span" in idx) {
          raw = $(idx["n_span"])
          if (!missing(raw)) {
            ns = raw + 0
            if (!seen_n_span || ns > max_n_span) { max_n_span = ns; max_n_span_t = t }
            seen_n_span = 1
          }
        }
        if ("min_Dphys" in idx) {
          raw = $(idx["min_Dphys"])
          if (!missing(raw)) {
            dp = raw + 0
            if (!seen_dphys || dp < min_dphys) { min_dphys = dp; min_dphys_t = t }
            seen_dphys = 1
          }
        }
        if ("min_gate_tot" in idx) {
          raw = $(idx["min_gate_tot"])
          if (!missing(raw)) {
            gt = raw + 0
            if (!seen_gate || gt < min_gate) { min_gate = gt; min_gate_t = t }
            seen_gate = 1
          }
        }
        if ("max_phi_cell" in idx) {
          raw = $(idx["max_phi_cell"])
          if (!missing(raw)) {
            pc = raw + 0
            if (!seen_phi || pc > max_phi_cell) { max_phi_cell = pc; max_phi_cell_t = t }
            seen_phi = 1
          }
        }
        seen = 1
      }
      END {
        if (!seen) {
          print "  (no timestep rows found)"
          exit
        }
        printf "  max nonlinear_its: %g at time=%g\n", max_nl, max_nl_t
        printf "  max linear_its: %g at time=%g\n", max_li, max_li_t
        printf "  min dt: %g at time=%g\n", min_dt, min_dt_t
        printf "  min dt_ratio: %g at time=%g\n", min_dtr, min_dtr_t
        if (seen_j)
          printf "  min J_min_1: %g at time=%g\n", min_j, min_j_t
        if (seen_n_span)
          printf "  max n_span: %g at time=%g\n", max_n_span, max_n_span_t
        if (seen_dphys)
          printf "  min min_Dphys: %g at time=%g\n", min_dphys, min_dphys_t
        if (seen_gate)
          printf "  min min_gate_tot: %g at time=%g\n", min_gate, min_gate_t
        if (seen_phi)
          printf "  max max_phi_cell: %g at time=%g\n", max_phi_cell, max_phi_cell_t
      }
    ' "${hotspots_csv}"
  } > "${report_txt}"
}

write_physics_sanity_report() {
  local hotspots_csv="$1"
  local report_txt="$2"
  local rc

  set +e
  awk -F, \
    -v n_floor="${SMOKE_N_MIN_FLOOR}" \
    -v j_floor="${SMOKE_J_MIN_FLOOR}" \
    -v d_floor="${SMOKE_DPHYS_MIN_FLOOR}" \
    -v gate_floor="${SMOKE_GATE_MIN_FLOOR}" \
    -v gate_ceil="${SMOKE_GATE_MAX_CEIL}" \
    '
      function missing(raw, low) {
        low = tolower(raw)
        return (raw == "" || low == "na" || low == "nan" || low == "+inf" || low == "-inf" || low == "inf")
      }
      NR == 1 {
        for (i = 1; i <= NF; i++)
          idx[$i] = i
        next
      }
      NR > 1 {
        seen = 1
        if ("time" in idx && !missing($(idx["time"]))) {
          t = $(idx["time"]) + 0
          last_t = t
        } else
          t = 0

        if ("n_min" in idx) {
          raw = $(idx["n_min"])
          if (!missing(raw)) {
            n = raw + 0
            if (!seen_n || n < min_n) { min_n = n; min_n_t = t; seen_n = 1 }
          }
        }
        if ("J_min_1" in idx) {
          raw = $(idx["J_min_1"])
          if (!missing(raw)) {
            j = raw + 0
            if (!seen_j || j < min_j) { min_j = j; min_j_t = t; seen_j = 1 }
          }
        }
        if ("min_Dphys" in idx) {
          raw = $(idx["min_Dphys"])
          if (!missing(raw)) {
            d = raw + 0
            if (!seen_d || d < min_d) { min_d = d; min_d_t = t; seen_d = 1 }
          }
        }
        if ("min_gate_tot" in idx) {
          raw = $(idx["min_gate_tot"])
          if (!missing(raw)) {
            gmin = raw + 0
            if (!seen_gmin || gmin < min_g) { min_g = gmin; min_g_t = t; seen_gmin = 1 }
          }
        }
        if ("max_gate_tot" in idx) {
          raw = $(idx["max_gate_tot"])
          if (!missing(raw)) {
            gmax = raw + 0
            if (!seen_gmax || gmax > max_g) { max_g = gmax; max_g_t = t; seen_gmax = 1 }
          }
        }
      }
      END {
        if (!seen) {
          print "Physics sanity checks: no timestep rows found."
          exit 0
        }

        fail = 0
        print "Physics sanity checks:"
        printf "  final time: %g\n", last_t

        if (seen_n) {
          printf "  min n_min: %g at time=%g\n", min_n, min_n_t
          if (min_n < n_floor) {
            printf "  FAIL: n_min (%g) < floor (%g)\n", min_n, n_floor
            fail = 1
          }
        }
        if (seen_j) {
          printf "  min J_min_1: %g at time=%g\n", min_j, min_j_t
          if (min_j <= j_floor) {
            printf "  FAIL: J_min_1 (%g) <= floor (%g)\n", min_j, j_floor
            fail = 1
          }
        }
        if (seen_d) {
          printf "  min min_Dphys: %g at time=%g\n", min_d, min_d_t
          if (min_d < d_floor) {
            printf "  FAIL: min_Dphys (%g) < floor (%g)\n", min_d, d_floor
            fail = 1
          }
        }
        if (seen_gmin) {
          printf "  min min_gate_tot: %g at time=%g\n", min_g, min_g_t
          if (min_g < gate_floor) {
            printf "  FAIL: min_gate_tot (%g) < floor (%g)\n", min_g, gate_floor
            fail = 1
          }
        }
        if (seen_gmax) {
          printf "  max max_gate_tot: %g at time=%g\n", max_g, max_g_t
          if (max_g > gate_ceil) {
            printf "  FAIL: max_gate_tot (%g) > ceil (%g)\n", max_g, gate_ceil
            fail = 1
          }
        }

        if (!fail)
          print "  PASS: sanity thresholds satisfied."

        exit fail
      }
    ' "${hotspots_csv}" > "${report_txt}"
  rc=$?
  set -e
  return "${rc}"
}

write_failure_signature_report() {
  local report_txt="$1"
  shift
  local logs=("$@")
  local fail=0
  local patterns='SNES_DIVERGED_|KSP_DIVERGED_|DIVERGED_FNORM_NAN|DIVERGED_LINEAR_SOLVE|Linear solve did not converge|Nonlinear solve did not converge|floating point exception|segmentation fault|inverted element|negative jacobian|zero pivot|(^|[^A-Za-z])(NaN|nan|NAN|Inf|inf|INF)([^A-Za-z]|$)'
  local log
  local hits

  {
    echo "Numerical/physics failure signature scan:"
    echo "  patterns=${patterns}"
    for log in "${logs[@]}"; do
      [[ -f "${log}" ]] || continue
      echo
      echo "Log: ${log}"
      hits="$(grep -En "${patterns}" "${log}" | head -n 40 || true)"
      if [[ -n "${hits}" ]]; then
        fail=1
        echo "${hits}"
      else
        echo "  (no known signatures found)"
      fi
    done
  } > "${report_txt}"

  return "${fail}"
}

summarize_watch_diagnostics() {
  local solver_csv="${RUN_DIR}/solver_watch.csv"
  local hotspots_csv="${RUN_DIR}/solver_hotspots.csv"
  local report_txt="${RUN_DIR}/stiffness_summary.txt"
  local physics_report="${RUN_DIR}/physics_sanity.txt"

  if ! copy_watch_csvs; then
    echo "WARNING: no outputs/solver_watch.csv or outputs/mesh_watch.csv found to summarize"
    return
  fi

  if [[ ! -f "${solver_csv}" ]]; then
    echo "WARNING: solver_watch.csv missing; skipping slowdown summary"
    return
  fi

  if ! build_solver_hotspots_csv "${solver_csv}" "${hotspots_csv}"; then
    echo "WARNING: failed to build solver hotspot CSV from ${solver_csv}"
    return
  fi

  write_slowdown_report "${hotspots_csv}" "${report_txt}"
  echo "Stiffness summary: ${report_txt}"

  if write_physics_sanity_report "${hotspots_csv}" "${physics_report}"; then
    echo "Physics sanity: ${physics_report}"
  else
    echo "WARNING: physics sanity checks detected threshold violations. See ${physics_report}"
    if [[ "${SMOKE_FAIL_ON_PHYSICS}" == "1" ]]; then
      echo "ERROR: failing smoke due to physics sanity violations (SMOKE_FAIL_ON_PHYSICS=1)." >&2
      exit 1
    fi
  fi
}

echo "Smoke debug run directory: ${RUN_DIR}"
echo "SMOKE_NP=${SMOKE_NP}  SMOKE_INPUT=${SMOKE_INPUT}  SMOKE_END_TIME=${SMOKE_END_TIME}"
echo "SMOKE_SLOW_NL_ITS=${SMOKE_SLOW_NL_ITS}  SMOKE_SLOW_LIN_ITS=${SMOKE_SLOW_LIN_ITS}  SMOKE_DT_DROP_FACTOR=${SMOKE_DT_DROP_FACTOR}"
echo "SMOKE_REQUIRE_CURRENT_INPUT=${SMOKE_REQUIRE_CURRENT_INPUT}  SMOKE_FAIL_ON_PHYSICS=${SMOKE_FAIL_ON_PHYSICS}  SMOKE_FAIL_ON_NUMERICS=${SMOKE_FAIL_ON_NUMERICS}"

require_smoke_input "${SMOKE_INPUT}" "SMOKE_INPUT"
if [[ "${SMOKE_RUN_JACOBIAN}" != "0" ]]; then
  require_smoke_input "${SMOKE_JAC_INPUT}" "SMOKE_JAC_INPUT"
fi

if ! run_and_log_allow_fail input_check \
  ./collie-opt \
  --error \
  --error-unused \
  --check-input \
  -i "${SMOKE_INPUT}"; then
  if is_mpi_sandbox_error "${RUN_DIR}/input_check.log"; then
    echo "ERROR: MPI/libfabric initialization is blocked by current sandbox permissions." >&2
    echo "Hint: restart with permissions that allow local listener sockets and POSIX shared memory." >&2
    exit 2
  fi
  echo "ERROR: input_check failed; see ${RUN_DIR}/input_check.log" >&2
  exit 1
fi

if ! run_and_log_allow_fail mpi_solver_physics \
  mpiexec -n "${SMOKE_NP}" ./collie-opt \
  -i "${SMOKE_INPUT}" \
  --error \
  --error-unused \
  Executioner/end_time="${SMOKE_END_TIME}" \
  Outputs/exodus=false \
  Outputs/perf_graph=false \
  Debug/show_var_residual_norms=true \
  Executioner/petsc_options='-snes_monitor -snes_converged_reason -ksp_monitor_short -ksp_converged_reason'; then
  if is_not_enough_slots_error "${RUN_DIR}/mpi_solver_physics.log" && [[ "${SMOKE_NP}" != "1" ]]; then
    echo "WARNING: MPI slots were insufficient for SMOKE_NP=${SMOKE_NP}; retrying with SMOKE_NP=1"
    run_and_log mpi_solver_physics_retry_np1 \
      mpiexec -n 1 ./collie-opt \
      -i "${SMOKE_INPUT}" \
      --error \
      --error-unused \
      Executioner/end_time="${SMOKE_END_TIME}" \
      Outputs/exodus=false \
      Outputs/perf_graph=false \
      Debug/show_var_residual_norms=true \
      Executioner/petsc_options='-snes_monitor -snes_converged_reason -ksp_monitor_short -ksp_converged_reason'
  else
    if is_mpi_sandbox_error "${RUN_DIR}/mpi_solver_physics.log"; then
      echo "ERROR: MPI/libfabric initialization is blocked by current sandbox permissions." >&2
      echo "Hint: restart with permissions that allow local listener sockets and POSIX shared memory." >&2
      exit 2
    fi
    echo "ERROR: mpi_solver_physics failed; see ${RUN_DIR}/mpi_solver_physics.log" >&2
    exit 1
  fi
fi

summarize_watch_diagnostics

signature_logs=(
  "${RUN_DIR}/input_check.log"
  "${RUN_DIR}/mpi_solver_physics.log"
)
if [[ -f "${RUN_DIR}/mpi_solver_physics_retry_np1.log" ]]; then
  signature_logs+=("${RUN_DIR}/mpi_solver_physics_retry_np1.log")
fi

if [[ "${SMOKE_RUN_JACOBIAN}" != "0" ]]; then
  jac_opts='-snes_monitor -snes_converged_reason -ksp_converged_reason -snes_test_jacobian'
  jac_extra_args=()

  if [[ "${SMOKE_JAC_VIEW}" == "1" ]]; then
    jac_opts="${jac_opts} -snes_test_jacobian_view"
  fi

  if [[ -n "${SMOKE_JAC_OVERRIDES}" ]]; then
    read -r -a jac_extra_args <<< "${SMOKE_JAC_OVERRIDES}"
  elif [[ "${SMOKE_JAC_AUTO_MESH_REDUCE}" == "1" ]] && \
       grep -Eq 'type[[:space:]]*=[[:space:]]*GeneratedMeshGenerator' "${SMOKE_JAC_INPUT}"; then
    jac_extra_args=("Mesh/gm/nx=${SMOKE_JAC_NX}" "Mesh/gm/ny=${SMOKE_JAC_NY}")
  fi

  echo
  echo "SMOKE_RUN_JACOBIAN=${SMOKE_RUN_JACOBIAN}  SMOKE_JAC_INPUT=${SMOKE_JAC_INPUT}  SMOKE_JAC_END_TIME=${SMOKE_JAC_END_TIME}"
  if [[ "${#jac_extra_args[@]}" -gt 0 ]]; then
    echo "SMOKE_JAC_OVERRIDES=${jac_extra_args[*]}"
  fi
  run_and_log jacobian_spotcheck \
    ./collie-opt \
    -i "${SMOKE_JAC_INPUT}" \
    --error \
    --error-unused \
    Executioner/end_time="${SMOKE_JAC_END_TIME}" \
    Outputs/exodus=false \
    Outputs/csv=false \
    Outputs/perf_graph=false \
    "${jac_extra_args[@]}" \
    Executioner/petsc_options="${jac_opts}"
  summarize_jacobian "${RUN_DIR}/jacobian_spotcheck.log"
  signature_logs+=("${RUN_DIR}/jacobian_spotcheck.log")
fi

if write_failure_signature_report "${RUN_DIR}/numerics_signals.txt" "${signature_logs[@]}"; then
  echo "Numerical signature scan: ${RUN_DIR}/numerics_signals.txt (no known fail signatures)"
else
  echo "WARNING: numerical failure signatures found. See ${RUN_DIR}/numerics_signals.txt"
  if [[ "${SMOKE_FAIL_ON_NUMERICS}" == "1" ]]; then
    echo "ERROR: failing smoke due to numerical signature hits (SMOKE_FAIL_ON_NUMERICS=1)." >&2
    exit 1
  fi
fi

echo
echo "Smoke debug complete. Logs:"
echo "  ${RUN_DIR}/input_check.log"
echo "  ${RUN_DIR}/mpi_solver_physics.log"
if [[ -f "${RUN_DIR}/solver_watch.csv" ]]; then
  echo "  ${RUN_DIR}/solver_watch.csv"
fi
if [[ -f "${RUN_DIR}/mesh_watch.csv" ]]; then
  echo "  ${RUN_DIR}/mesh_watch.csv"
fi
if [[ -f "${RUN_DIR}/solver_hotspots.csv" ]]; then
  echo "  ${RUN_DIR}/solver_hotspots.csv"
fi
if [[ -f "${RUN_DIR}/stiffness_summary.txt" ]]; then
  echo "  ${RUN_DIR}/stiffness_summary.txt"
fi
if [[ -f "${RUN_DIR}/physics_sanity.txt" ]]; then
  echo "  ${RUN_DIR}/physics_sanity.txt"
fi
if [[ -f "${RUN_DIR}/numerics_signals.txt" ]]; then
  echo "  ${RUN_DIR}/numerics_signals.txt"
fi
if [[ "${SMOKE_RUN_JACOBIAN}" != "0" ]]; then
  echo "  ${RUN_DIR}/jacobian_spotcheck.log"
fi
