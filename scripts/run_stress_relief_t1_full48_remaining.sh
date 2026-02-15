#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${ROOT}"

NP="${NP:-8}"
MPI_LAUNCHER="${MPI_LAUNCHER:-/home/amcgs/miniforge/envs/moose/bin/mpiexec}"
CASE_DIR="inputs/legacy/2DN/suites/stress_relief_t1/permutations_rz3"
OUT_ROOT="outputs/suites/stress_relief_t1"
QUICK_ANALYSIS="${OUT_ROOT}/quick30_mesh60/analysis_table.txt"
FULL_DIR="${OUT_ROOT}/full48"

STAMP="$(date +%Y%m%d_%H%M%S)"
BACKUP_DIR="${FULL_DIR}/backup_${STAMP}"
SUMMARY="${FULL_DIR}/summary_resume.txt"
ANALYSIS="${FULL_DIR}/analysis_table.txt"
CONSOLIDATED="${FULL_DIR}/consolidated_status_${STAMP}.txt"

mkdir -p "${FULL_DIR}" "${BACKUP_DIR}"

if [[ ! -f "${QUICK_ANALYSIS}" ]]; then
  echo "Missing quick analysis file: ${QUICK_ANALYSIS}" >&2
  exit 1
fi

if [[ ! -x "${MPI_LAUNCHER}" ]]; then
  MPI_LAUNCHER="$(command -v mpiexec || command -v mpirun || true)"
fi
if [[ -z "${MPI_LAUNCHER}" ]]; then
  echo "No MPI launcher found (checked MPI_LAUNCHER, mpiexec, mpirun)." >&2
  exit 1
fi

pick_top_two() {
  local analysis="$1"
  awk -F'\\|' '
    function trim(s) { gsub(/^[[:space:]]+|[[:space:]]+$/, "", s); return s }
    NR == 1 { next }
    {
      case_name = trim($1)
      status = trim($2)
      score = trim($21) + 0
      if (status == "success")
        printf "%s|%.6f\n", case_name, score
    }
  ' "${analysis}" | sort -t'|' -k2,2nr | head -n 2 | cut -d'|' -f1
}

parse_solver_csv() {
  local csv="$1"
  awk -F, '
    NR == 1 {
      for (i = 1; i <= NF; ++i)
        h[$i] = i
      next
    }
    {
      seen = 1

      t   = ("time" in h) ? $(h["time"]) + 0 : 0
      li  = ("linear_its" in h) ? $(h["linear_its"]) + 0 : 0
      ni  = ("nonlinear_its" in h) ? $(h["nonlinear_its"]) + 0 : 0
      meq = ("min_elem_quality" in h) ? $(h["min_elem_quality"]) + 0 : 0
      j1  = ("J_min_1" in h) ? $(h["J_min_1"]) + 0 : 0
      aj  = ("avg_J" in h) ? $(h["avg_J"]) + 0 : 0
      vc  = ("vol_change_pct" in h) ? $(h["vol_change_pct"]) + 0 : 0
      ap  = ("avg_pressure" in h) ? $(h["avg_pressure"]) + 0 : (("avg_press" in h) ? $(h["avg_press"]) + 0 : 0)
      agr = ("avg_gp_raw" in h) ? $(h["avg_gp_raw"]) + 0 : (("avg_gp" in h) ? $(h["avg_gp"]) + 0 : 0)
      agk = ("avg_gp_ke_used" in h) ? $(h["avg_gp_ke_used"]) + 0 : (("avg_gp" in h) ? $(h["avg_gp"]) + 0 : 0)
      ach = ("avg_chi" in h) ? $(h["avg_chi"]) + 0 : 0
      mchi = ("max_chi" in h) ? $(h["max_chi"]) + 0 : ach
      nmin = ("n_min" in h) ? $(h["n_min"]) + 0 : (("n_min_elem" in h) ? $(h["n_min_elem"]) + 0 : 0)
      nco  = ("n_corner" in h) ? $(h["n_corner"]) + 0 : 0
      nbu  = ("n_bulk" in h) ? $(h["n_bulk"]) + 0 : 0
      chco = ("chi_corner" in h) ? $(h["chi_corner"]) + 0 : 0
      chbu = ("chi_bulk" in h) ? $(h["chi_bulk"]) + 0 : 0

      if (!have_li || li > max_li) { max_li = li; have_li = 1 }
      if (li >= 30 && t_li30 == "") t_li30 = t
      if (!have_meq || meq < min_meq) { min_meq = meq; have_meq = 1 }
      if (!have_j1 || j1 < min_j1) { min_j1 = j1; have_j1 = 1 }
      if (!have_mchi || mchi > max_chi) { max_chi = mchi; t_max_chi = t; have_mchi = 1 }

      last_t = t
      last_li = li
      last_ni = ni
      last_avg_j = aj
      last_vc = vc
      last_ap = ap
      last_agr = agr
      last_agk = agk
      last_ach = ach
      last_nmin = nmin
      last_nco = nco
      last_nbu = nbu
      last_chco = chco
      last_chbu = chbu
    }
    END {
      if (!seen)
      {
        print "NA|NA|NA|NA|NA|NA|NA|NA|NA|NA|NA|NA|NA|NA|NA|NA|NA|NA|NA|NA"
        exit
      }

      if (t_li30 == "")
        t_li30 = "NA"

      printf "%.8g|%d|%s|%d|%d|%.12g|%.8g|%.8g|%.8g|%.8g|%.8g|%.8g|%.8g|%.8g|%.8g|%.8g|%.8g|%.8g|%.8g|%.8g\n", \
        last_t, max_li, t_li30, last_li, last_ni, min_meq, min_j1, last_avg_j, last_vc, last_ap, \
        last_agr, last_agk, max_chi, t_max_chi, last_ach, last_nmin, last_nco, last_nbu, last_chco, last_chbu
    }
  ' "${csv}"
}

init_summary() {
  local summary="$1"
  {
    echo "case_name | status | wall_s | final_time | max_li | t_first_li30 | final_li | final_ni | min_elem_quality | min_J | avg_J | vol_change_pct | avg_pressure | avg_gp_raw | avg_gp_ke_used | max_chi | t_max_chi | avg_chi | min_n | corner_n | bulk_n | chi_corner | chi_bulk"
  } > "${summary}"
}

build_analysis() {
  local summary="$1"
  local analysis="$2"

  awk -F'\\|' '
    function trim(s) { gsub(/^[[:space:]]+|[[:space:]]+$/, "", s); return s }
    NR == 1 {
      print "case_name | status | max_li | t_first_li30 | final_time | final_li | min_elem_quality | min_J | avg_J | vol_change_pct | avg_pressure | avg_gp_raw | avg_gp_ke_used | max_chi | avg_chi | min_n | corner_n | bulk_n | chi_corner | chi_bulk | score"
      next
    }
    {
      case_name = trim($1)
      status = trim($2)
      final_time = trim($4)
      max_li = trim($5)
      tli30 = trim($6)
      final_li = trim($7)
      min_eq = trim($9)
      min_j = trim($10)
      avg_j = trim($11)
      vc = trim($12)
      avg_p = trim($13)
      avg_gpr = trim($14)
      avg_gpk = trim($15)
      max_chi = trim($16)
      avg_chi = trim($18)
      min_n = trim($19)
      corner_n = trim($20)
      bulk_n = trim($21)
      chi_corner = trim($22)
      chi_bulk = trim($23)

      if (case_name == "00_baseline")
      {
        b_max_li = max_li + 0
        b_tli30 = (tli30 == "NA") ? 1e9 : (tli30 + 0)
      }

      score = -1e9
      if (status == "success" && final_time != "NA" && (final_time + 0) >= 29.9)
      {
        score = 0

        if (tli30 == "NA")
          score += 3
        else if ((tli30 + 0) > b_tli30)
          score += 1
        else
          score -= 2

        if ((max_li + 0) < b_max_li)
          score += 2
        if ((max_li + 0) <= 20)
          score += 1

        if ((min_n + 0) >= -1e-8 && (corner_n + 0) >= -1e-8)
          score += 2
        else
          score -= 4

        if ((vc + 0) >= 1200)
          score += 2
        else if ((vc + 0) >= 800)
          score += 1
        else
          score -= 2

        if ((min_eq + 0) >= 1e-4)
          score += 1
        else if ((min_eq + 0) < 1e-5)
          score -= 3
      }

      printf "%s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %.3f\n", \
             case_name, status, max_li, tli30, final_time, final_li, min_eq, min_j, avg_j, vc, avg_p, avg_gpr, avg_gpk, max_chi, avg_chi, min_n, corner_n, bulk_n, chi_corner, chi_bulk, score
    }
  ' "${summary}" > "${analysis}"
}

run_case() {
  local case_name="$1"
  local case_file="${CASE_DIR}/${case_name}.i"
  local case_out="${FULL_DIR}/${case_name}"
  local log_file="${case_out}/run.log"
  local solver_csv="${case_out}/solver_watch.csv"

  mkdir -p "${case_out}"

  local start_ts
  local end_ts
  local wall
  local rc
  local status
  local metrics

  local -a cmd=("${MPI_LAUNCHER}")
  cmd+=(
    -n "${NP}" ./collie-opt -i "${case_file}" --output-formatter perf_graph
    "Executioner/end_time=48"
    "Outputs/exodus=true"
    "Outputs/interval=20"
    "Outputs/file_base=${case_out}/${case_name}"
    "Outputs/mesh_watch/file_base=${case_out}/mesh_watch"
    "Outputs/solver_watch/file_base=${case_out}/solver_watch"
  )

  echo "=== [${case_name}] mode=full48 ==="
  start_ts="$(date +%s)"
  set +e
  "${cmd[@]}" 2>&1 | tee "${log_file}"
  rc="${PIPESTATUS[0]}"
  set -e
  end_ts="$(date +%s)"
  wall=$((end_ts - start_ts))

  if [[ "${rc}" -eq 0 ]]; then
    status="success"
  else
    status="fail"
  fi

  if [[ -f "${solver_csv}" ]]; then
    metrics="$(parse_solver_csv "${solver_csv}")"
  else
    metrics="NA|NA|NA|NA|NA|NA|NA|NA|NA|NA|NA|NA|NA|NA|NA|NA|NA|NA|NA|NA"
  fi

  printf "%s | %s | %s | %s\n" "${case_name}" "${status}" "${wall}" "${metrics}" >> "${SUMMARY}"
  return "${rc}"
}

write_consolidated_status() {
  local started_at="$1"
  local finished_at="$2"
  shift 2
  local cases=("$@")

  local ok_count
  local fail_count
  ok_count="$(awk -F'|' 'NR>1{gsub(/^[[:space:]]+|[[:space:]]+$/, "", $2); if ($2=="success") c++} END{print c+0}' "${SUMMARY}")"
  fail_count="$(awk -F'|' 'NR>1{gsub(/^[[:space:]]+|[[:space:]]+$/, "", $2); if ($2=="fail") c++} END{print c+0}' "${SUMMARY}")"

  {
    echo "stress_relief_t1 full48 completion status"
    echo "started_at: ${started_at}"
    echo "finished_at: ${finished_at}"
    echo "np: ${NP}"
    echo "mpi_launcher: ${MPI_LAUNCHER}"
    echo "root: ${ROOT}"
    echo "backup_dir: ${BACKUP_DIR}"
    echo
    echo "selected_cases:"
    for case_name in "${cases[@]}"; do
      echo "- ${case_name}"
    done
    echo
    echo "result_counts:"
    echo "- success: ${ok_count}"
    echo "- fail: ${fail_count}"
    echo
    echo "summary_file: ${SUMMARY}"
    echo "analysis_file: ${ANALYSIS}"
    echo
    echo "summary_rows:"
    tail -n +1 "${SUMMARY}"
  } > "${CONSOLIDATED}"

  cp -f "${SUMMARY}" "${FULL_DIR}/summary_resume_${STAMP}.txt"
  cp -f "${ANALYSIS}" "${FULL_DIR}/analysis_table_${STAMP}.txt"
}

mapfile -t TOP2 < <(pick_top_two "${QUICK_ANALYSIS}")
if [[ "${#TOP2[@]}" -eq 0 ]]; then
  echo "No successful quick cases available from ${QUICK_ANALYSIS}" >&2
  exit 1
fi

# Back up existing full48 bookkeeping and selected case folders.
for f in summary.txt summary_resume.txt analysis_table.txt resume_runner.log; do
  if [[ -f "${FULL_DIR}/${f}" ]]; then
    mv "${FULL_DIR}/${f}" "${BACKUP_DIR}/${f}"
  fi
done

for case_name in "${TOP2[@]}"; do
  if [[ -d "${FULL_DIR}/${case_name}" ]]; then
    mv "${FULL_DIR}/${case_name}" "${BACKUP_DIR}/${case_name}"
  fi
done

started_at="$(date -Iseconds)"
init_summary "${SUMMARY}"

overall_rc=0
for case_name in "${TOP2[@]}"; do
  set +e
  run_case "${case_name}"
  rc=$?
  set -e
  if [[ "${rc}" -ne 0 ]]; then
    overall_rc=1
  fi
done

build_analysis "${SUMMARY}" "${ANALYSIS}"
finished_at="$(date -Iseconds)"
write_consolidated_status "${started_at}" "${finished_at}" "${TOP2[@]}"

echo
echo "Summary:      ${SUMMARY}"
echo "Analysis:     ${ANALYSIS}"
echo "Consolidated: ${CONSOLIDATED}"
echo "Backup:       ${BACKUP_DIR}"

exit "${overall_rc}"
