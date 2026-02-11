#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${ROOT}"

NP="${NP:-8}"
CASE_DIR="${CASE_DIR:-inputs/2DN/suites/jacobian_solver/permutations_rz3}"
OUT_ROOT="${OUT_ROOT:-outputs/suites/jacobian_solver/permutations_rz3}"
QUICK="${QUICK:-0}"
SUMMARY="${OUT_ROOT}/summary.txt"
RUN_OVERRIDES="${RUN_OVERRIDES:-}"
QUIET="${QUIET:-0}"

EXTRA_ARGS=()
if [[ -n "${RUN_OVERRIDES}" ]]; then
  # Space-delimited CLI overrides, e.g.
  # RUN_OVERRIDES="Executioner/end_time=30 Outputs/exodus=false"
  read -r -a EXTRA_ARGS <<< "${RUN_OVERRIDES}"
fi

mkdir -p "${OUT_ROOT}"

if [[ ! -d "${CASE_DIR}" ]]; then
  echo "ERROR: case directory not found: ${CASE_DIR}" >&2
  exit 2
fi

mapfile -t CASE_FILES < <(find "${CASE_DIR}" -maxdepth 1 -type f -name '*.i' | sort)

if [[ "${#CASE_FILES[@]}" -eq 0 ]]; then
  echo "ERROR: no .i files found in ${CASE_DIR}" >&2
  exit 2
fi

if [[ "${QUICK}" == "1" ]]; then
  mapfile -t CASE_FILES < <(
    for f in "${CASE_FILES[@]}"; do
      b="$(basename "${f}")"
      n="${b%%_*}"
      if [[ "${n}" =~ ^[0-9][0-9]$ ]] && ((10#${n} <= 7)); then
        echo "${f}"
      fi
    done
  )
fi

{
  echo "case_name | status | wall_time_s | last_time | min_elem_quality | min_volume_ratio | avg_volume_ratio | avg_phi_cell | avg_ke_total | corner_n | bulk_n | corner_ke | bulk_ke"
} > "${SUMMARY}"

echo "Running ${#CASE_FILES[@]} cases from ${CASE_DIR} (NP=${NP}, QUICK=${QUICK})"

for case_file in "${CASE_FILES[@]}"; do
  case_name="$(basename "${case_file}" .i)"
  case_out="${OUT_ROOT}/${case_name}"
  log_file="${case_out}/run.log"
  solver_csv="${case_out}/solver_watch.csv"

  mkdir -p "${case_out}"

  echo
  echo "=== [${case_name}] ==="
  start_ts="$(date +%s)"

  set +e
  if [[ "${QUIET}" == "1" ]]; then
    mpirun -np "${NP}" ./collie-opt -i "${case_file}" --output-formatter perf_graph "${EXTRA_ARGS[@]}" > "${log_file}" 2>&1
    rc=$?
  else
    mpirun -np "${NP}" ./collie-opt -i "${case_file}" --output-formatter perf_graph "${EXTRA_ARGS[@]}" 2>&1 | tee "${log_file}"
    rc=${PIPESTATUS[0]}
  fi
  set -e

  end_ts="$(date +%s)"
  wall_time=$((end_ts - start_ts))

  if [[ ${rc} -eq 0 ]]; then
    status="success"
  else
    status="fail"
  fi

  last_time="NA"
  min_elem_quality="NA"
  min_volume_ratio="NA"
  avg_volume_ratio="NA"
  avg_phi_cell="NA"
  avg_ke_total="NA"
  corner_n="NA"
  bulk_n="NA"
  corner_ke="NA"
  bulk_ke="NA"

  if [[ ! -f "${solver_csv}" ]]; then
    solver_base="$(
      awk '
        BEGIN { in_solver = 0 }
        /^\s*\[solver_watch\]\s*$/ { in_solver = 1; next }
        in_solver && /^\s*\[\]\s*$/ { in_solver = 0 }
        in_solver && /^\s*file_base\s*=/ {
          sub(/^[^=]*=/, "", $0)
          gsub(/^[[:space:]]+|[[:space:]]+$/, "", $0)
          gsub(/^'\''|'\''$/, "", $0)
          gsub(/^"|"$/, "", $0)
          print $0
          exit
        }
      ' "${case_file}"
    )"
    if [[ -n "${solver_base}" && -f "${solver_base}.csv" ]]; then
      solver_csv="${solver_base}.csv"
    fi
  fi

  if [[ -f "${solver_csv}" ]]; then
    read -r last_time min_elem_quality min_volume_ratio avg_volume_ratio avg_phi_cell avg_ke_total corner_n bulk_n corner_ke bulk_ke < <(
      awk -F, '
        function idx_of(name, i) {
          for (i = 1; i <= NF; ++i)
            if ($i == name)
              return i
          return 0
        }
        NR == 1 {
          i_time = idx_of("time")
          i_meq = idx_of("min_elem_quality")
          i_mvr = idx_of("min_volume_ratio")
          i_avgj = idx_of("avg_J")
          i_aphi = idx_of("avg_phi_cell")
          i_aket = idx_of("avg_ke_total")
          i_cn = idx_of("n_corner")
          i_bn = idx_of("n_bulk")
          i_ck = idx_of("ke_total_corner")
          i_bk = idx_of("ke_total_bulk")
          next
        }
        NR >= 2 {
          if (i_time) last_time = $i_time
          if (i_avgj) avgj = $i_avgj
          if (i_aphi) aphi = $i_aphi
          if (i_aket) aket = $i_aket
          if (i_cn) cn = $i_cn
          if (i_bn) bn = $i_bn
          if (i_ck) ck = $i_ck
          if (i_bk) bk = $i_bk

          if (i_meq && (min_meq == "" || $i_meq < min_meq)) min_meq = $i_meq
          if (i_mvr && (min_mvr == "" || $i_mvr < min_mvr)) min_mvr = $i_mvr
        }
        END {
          if (last_time == "") last_time = "NA"
          if (min_meq == "") min_meq = "NA"
          if (min_mvr == "") min_mvr = "NA"
          if (avgj == "") avgj = "NA"
          if (aphi == "") aphi = "NA"
          if (aket == "") aket = "NA"
          if (cn == "") cn = "NA"
          if (bn == "") bn = "NA"
          if (ck == "") ck = "NA"
          if (bk == "") bk = "NA"
          printf "%s %s %s %s %s %s %s %s %s %s\n", last_time, min_meq, min_mvr, avgj, aphi, aket, cn, bn, ck, bk
        }
      ' "${solver_csv}"
    )
  fi

  printf "%s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s\n" \
    "${case_name}" "${status}" "${wall_time}" "${last_time}" "${min_elem_quality}" "${min_volume_ratio}" \
    "${avg_volume_ratio}" "${avg_phi_cell}" "${avg_ke_total}" "${corner_n}" "${bulk_n}" "${corner_ke}" "${bulk_ke}" \
    >> "${SUMMARY}"

done

echo
echo "Summary written to ${SUMMARY}"
