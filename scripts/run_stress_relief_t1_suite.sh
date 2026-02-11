#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${ROOT}"

NP="${NP:-8}"
CASE_DIR="inputs/2DN/suites/stress_relief_t1/permutations_rz3"
OUT_ROOT="outputs/suites/stress_relief_t1"
QUICK_DIR="${OUT_ROOT}/quick30_mesh60"
FULL_DIR="${OUT_ROOT}/full48"

mkdir -p "${QUICK_DIR}" "${FULL_DIR}"

CASES=(
  00_baseline
  01_gp_to_ke_OFF
  02_gp_to_ke_LAGGED_filter_tau0p5h
  03_press_str_x0p5
  04_press_str_x2
  05_T1_on_kT1max5_chi0p22_beta120
  06_T1_on_kT1max5_PLUS_gpke_lag_filter_tau0p5h
)

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

run_case() {
  local mode_dir="$1"
  local case_name="$2"
  local end_time="$3"
  local exodus_flag="$4"
  local mesh_override="$5"
  local summary="$6"

  local case_file="${CASE_DIR}/${case_name}.i"
  local case_out="${mode_dir}/${case_name}"
  local log_file="${case_out}/run.log"
  local solver_csv="${case_out}/solver_watch.csv"

  mkdir -p "${case_out}"

  local start_ts
  local end_ts
  local wall
  local rc
  local status
  local metrics

  local -a cmd=(
    mpirun -np "${NP}" ./collie-opt -i "${case_file}" --output-formatter perf_graph
    "Executioner/end_time=${end_time}"
    "Outputs/exodus=${exodus_flag}"
    "Outputs/file_base=${case_out}/${case_name}"
    "Outputs/mesh_watch/file_base=${case_out}/mesh_watch"
    "Outputs/solver_watch/file_base=${case_out}/solver_watch"
  )

  if [[ "${mesh_override}" == "1" ]]; then
    cmd+=("Mesh/gm/nx=60" "Mesh/gm/ny=18")
  fi

  if [[ "${exodus_flag}" == "true" ]]; then
    cmd+=("Outputs/interval=20")
  fi

  echo "=== [${case_name}] mode=$(basename "${mode_dir}") ==="
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

  printf "%s | %s | %s | %s\n" "${case_name}" "${status}" "${wall}" "${metrics}" >> "${summary}"
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
      wall = trim($3)
      final_time = trim($4)
      max_li = trim($5)
      tli30 = trim($6)
      final_li = trim($7)
      final_ni = trim($8)
      min_eq = trim($9)
      min_j = trim($10)
      avg_j = trim($11)
      vc = trim($12)
      avg_p = trim($13)
      avg_gpr = trim($14)
      avg_gpk = trim($15)
      max_chi = trim($16)
      t_max_chi = trim($17)
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

        if ((t_max_chi + 0) >= 18.0)
          score += 1

        if ((chi_bulk + 0) > 1e-12)
        {
          ratio = (chi_corner + 0) / (chi_bulk + 0)
          if (ratio > 3.0)
            score -= 1
        }
      }

      printf "%s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %.3f\n", \
             case_name, status, max_li, tli30, final_time, final_li, min_eq, min_j, avg_j, vc, avg_p, avg_gpr, avg_gpk, max_chi, avg_chi, min_n, corner_n, bulk_n, chi_corner, chi_bulk, score
    }
  ' "${summary}" > "${analysis}"
}

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

QUICK_SUMMARY="${QUICK_DIR}/summary.txt"
QUICK_ANALYSIS="${QUICK_DIR}/analysis_table.txt"
FULL_SUMMARY="${FULL_DIR}/summary.txt"
FULL_ANALYSIS="${FULL_DIR}/analysis_table.txt"

init_summary "${QUICK_SUMMARY}"
for case_name in "${CASES[@]}"; do
  run_case "${QUICK_DIR}" "${case_name}" "30" "false" "1" "${QUICK_SUMMARY}"
done
build_analysis "${QUICK_SUMMARY}" "${QUICK_ANALYSIS}"

mapfile -t TOP2 < <(pick_top_two "${QUICK_ANALYSIS}")

init_summary "${FULL_SUMMARY}"
if [[ "${#TOP2[@]}" -eq 0 ]]; then
  echo "No successful quick cases; skipping full48."
else
  for case_name in "${TOP2[@]}"; do
    run_case "${FULL_DIR}" "${case_name}" "48" "true" "0" "${FULL_SUMMARY}"
  done
fi
build_analysis "${FULL_SUMMARY}" "${FULL_ANALYSIS}"

echo ""
echo "Quick summary:   ${QUICK_SUMMARY}"
echo "Quick analysis:  ${QUICK_ANALYSIS}"
echo "Full summary:    ${FULL_SUMMARY}"
echo "Full analysis:   ${FULL_ANALYSIS}"
