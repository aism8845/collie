#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
cd "${ROOT}"

STAMP="${DEBUG_STAMP:-$(date +%Y%m%d_%H%M%S)}"
OUT_ROOT="${DEBUG_OUT_ROOT:-outputs/debug_suite}"
mkdir -p "${OUT_ROOT}"
export DEBUG_STAMP="${STAMP}"

if [[ -n "${DEBUG_CASES:-}" ]]; then
  # Space-separated list of case scripts, e.g.:
  # DEBUG_CASES="scripts/suites/physics_material/cases/debug_case_01_nutrient_gating_off.sh scripts/suites/physics_material/cases/debug_case_05_newton_lu_solver.sh"
  IFS=' ' read -r -a CASE_SCRIPTS <<< "${DEBUG_CASES}"
else
  CASE_SCRIPTS=(
    scripts/suites/physics_material/cases/debug_case_01_nutrient_gating_off.sh
    scripts/suites/physics_material/cases/debug_case_02_pressure_gating_off.sh
    scripts/suites/physics_material/cases/debug_case_03_smooth_gate_transition.sh
    scripts/suites/physics_material/cases/debug_case_04_diffusivity_floor_up.sh
    scripts/suites/physics_material/cases/debug_case_05_newton_lu_solver.sh
    scripts/suites/physics_material/cases/debug_case_06_constant_small_dt.sh
  )
fi

SUMMARY_CSV="${OUT_ROOT}/${STAMP}_summary.csv"
SUMMARY_TXT="${OUT_ROOT}/${STAMP}_summary.txt"

echo "case,status,rc,max_nonlinear_its,max_linear_its,min_dt,min_J_min_1,min_metric2_min_1,min_gate_tot,min_Dphys,max_n_span,max_phi_cell" > "${SUMMARY_CSV}"

summarize_solver_watch() {
  local csv="$1"

  awk -F, '
    NR == 1 {
      for (i = 1; i <= NF; i++)
        idx[$i] = i
      next
    }
    {
      seen = 1
      nl = $(idx["nonlinear_its"]) + 0
      li = $(idx["linear_its"]) + 0
      dt = $(idx["dt"]) + 0

      if (!have_nl || nl > max_nl) { max_nl = nl; have_nl = 1 }
      if (!have_li || li > max_li) { max_li = li; have_li = 1 }
      if (!have_dt || dt < min_dt) { min_dt = dt; have_dt = 1 }

      if ("J_min_1" in idx) {
        j = $(idx["J_min_1"]) + 0
        if (!have_j || j < min_j) { min_j = j; have_j = 1 }
      }

      if ("metric2_min_1" in idx) {
        m2 = $(idx["metric2_min_1"]) + 0
        if (!have_m2 || m2 < min_m2) { min_m2 = m2; have_m2 = 1 }
      }

      if ("min_gate_tot" in idx) {
        gt = $(idx["min_gate_tot"]) + 0
        if (!have_gt || gt < min_gt) { min_gt = gt; have_gt = 1 }
      }

      if ("min_Dphys" in idx) {
        dp = $(idx["min_Dphys"]) + 0
        if (!have_dp || dp < min_dp) { min_dp = dp; have_dp = 1 }
      }

      if ("n_span" in idx) {
        ns = $(idx["n_span"]) + 0
        if (!have_ns || ns > max_ns) { max_ns = ns; have_ns = 1 }
      }

      if ("max_phi_cell" in idx) {
        pc = $(idx["max_phi_cell"]) + 0
        if (!have_pc || pc > max_pc) { max_pc = pc; have_pc = 1 }
      }
    }
    END {
      if (!seen) {
        print "NA,NA,NA,NA,NA,NA,NA,NA,NA"
        exit
      }

      printf "%g,%g,%g,", max_nl, max_li, min_dt
      if (have_j) printf "%g", min_j; else printf "NA"
      printf ","
      if (have_m2) printf "%g", min_m2; else printf "NA"
      printf ","
      if (have_gt) printf "%g", min_gt; else printf "NA"
      printf ","
      if (have_dp) printf "%g", min_dp; else printf "NA"
      printf ","
      if (have_ns) printf "%g", max_ns; else printf "NA"
      printf ","
      if (have_pc) printf "%g\n", max_pc; else printf "NA\n"
    }
  ' "${csv}"
}

OVERALL_RC=0

echo "Debug campaign stamp: ${STAMP}"
echo "Output root: ${OUT_ROOT}"
echo

for case_script in "${CASE_SCRIPTS[@]}"; do
  case_id="$(basename "${case_script}" .sh)"
  case_id="${case_id#debug_case_}"
  case_dir="${OUT_ROOT}/${STAMP}_${case_id}"
  solver_csv="${case_dir}/solver_watch.csv"

  echo "==> Running case ${case_id}"

  set +e
  "${case_script}"
  rc=$?
  set -e

  status="PASS"
  if [[ ${rc} -ne 0 ]]; then
    status="FAIL"
    OVERALL_RC=1
  fi

  if [[ -f "${solver_csv}" ]]; then
    metrics="$(summarize_solver_watch "${solver_csv}")"
  else
    metrics="NA,NA,NA,NA,NA,NA,NA,NA,NA"
  fi

  echo "${case_id},${status},${rc},${metrics}" >> "${SUMMARY_CSV}"
done

{
  echo "Debug campaign stamp: ${STAMP}"
  echo "Output root: ${OUT_ROOT}"
  echo
  echo "Case intent:"
  echo "  01_nutrient_gating_off: disable nutrient gate sensitivity (gamma_n0 -> 0)"
  echo "  02_pressure_gating_off: disable pressure gating (press_str -> 0)"
  echo "  03_smooth_gate_transition: widen smooth gate transition (smooth_eps_c -> 1e-3)"
  echo "  04_diffusivity_floor_up: raise diffusivity floor (D_floor -> 1e-6)"
  echo "  05_newton_lu_solver: switch to direct Newton+LU solve"
  echo "  06_constant_small_dt: reduce initial adaptive dt (TimeStepper/dt -> 0.005)"
  echo
  echo "Interpretation hints:"
  echo "  lower max_nonlinear_its and max_linear_its: easier nonlinear/linear solve"
  echo "  higher min_J_min_1 and min_metric2_min_1: less mesh distortion risk"
  echo "  higher min_Dphys: less near-degenerate diffusion tensor risk"
  echo "  lower max_n_span: weaker nutrient fronts / lower reaction-diffusion stiffness"
  echo
  if command -v column >/dev/null 2>&1; then
    column -s, -t "${SUMMARY_CSV}"
  else
    cat "${SUMMARY_CSV}"
  fi
} | tee "${SUMMARY_TXT}"

echo
echo "Summary CSV: ${SUMMARY_CSV}"
echo "Summary TXT: ${SUMMARY_TXT}"

exit "${OVERALL_RC}"
