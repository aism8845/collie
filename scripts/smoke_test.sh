#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
EXE="${ROOT}/collie-opt"
INP="${ROOT}/inputs/2DN/RZ_smoke.i"
SMOKE_DRY_RUN="${SMOKE_DRY_RUN:-0}"
SMOKE_DEBUG="${SMOKE_DEBUG:-0}"
SMOKE_NP="${SMOKE_NP:-8}"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --dry-run)
      SMOKE_DRY_RUN=1
      shift
      ;;
    --debug)
      SMOKE_DEBUG=1
      shift
      ;;
    --np)
      SMOKE_NP="${2:?--np requires a value}"
      shift 2
      ;;
    *)
      echo "[smoke] ERROR: unknown argument: $1"
      echo "[smoke] usage: $0 [--dry-run] [--debug] [--np N]"
      exit 2
      ;;
  esac
done

# Ensure conda works in non-interactive shells
if command -v conda >/dev/null 2>&1; then
  source "$(conda info --base)/etc/profile.d/conda.sh"
fi

set +u
conda activate moose
set -u

echo "[smoke] build"
cd "${ROOT}"
make -j8

if [[ ! -x "${EXE}" ]]; then
  echo "[smoke] ERROR: executable not found: ${EXE}"
  exit 1
fi

TMP="$(mktemp -d)"
trap 'echo "[smoke] (keeping tmp dir on failure): ${TMP}"; [[ -f "${TMP}/smoke.log" ]] && tail -n 80 "${TMP}/smoke.log";' ERR
echo "[smoke] run dir: ${TMP}"

cp "${INP}" "${TMP}/RZ_smoke.i"

cd "${TMP}"

launchers=()
if [[ -n "${SMOKE_MPIEXEC:-}" ]]; then
  launchers+=("${SMOKE_MPIEXEC}")
else
  command -v mpiexec >/dev/null 2>&1 && launchers+=("mpiexec")
  command -v mpirun  >/dev/null 2>&1 && launchers+=("mpirun")
fi

if [[ ${#launchers[@]} -eq 0 ]]; then
  echo "[smoke] ERROR: no MPI launcher found (set SMOKE_MPIEXEC to override)"
  exit 1
fi

echo "[smoke] MPI run: SMOKE_NP=${SMOKE_NP}"

if [[ "${SMOKE_DEBUG}" == "1" || "${SMOKE_DRY_RUN}" == "1" ]]; then
  echo "[smoke] debug: PATH=${PATH}"
  echo "[smoke] debug: launcher candidates: ${launchers[*]}"
  for l in "${launchers[@]}"; do
    echo "[smoke] debug: launcher '${l}' -> $(command -v "${l}" 2>/dev/null || echo 'not in PATH (custom override?)')"
  done
  echo "[smoke] debug: exe=${EXE}"
  echo "[smoke] debug: input=${TMP}/RZ_smoke.i"
fi

if [[ "${SMOKE_DRY_RUN}" == "1" ]]; then
  echo "[smoke] dry-run: would run one of:"
  for l in "${launchers[@]}"; do
    echo "[smoke] dry-run: ${l} -n ${SMOKE_NP} ${EXE} -i RZ_smoke.i"
  done
  exit 0
fi

rm -f smoke.log
run_ok=0
for l in "${launchers[@]}"; do
  echo "[smoke] trying launcher: ${l}" | tee -a smoke.log
  set +e
  "${l}" -n "${SMOKE_NP}" "${EXE}" -i RZ_smoke.i >> smoke.log 2>&1
  ec=$?
  set -e
  if [[ $ec -eq 0 ]]; then
    run_ok=1
    break
  fi
  echo "[smoke] launcher failed (${l}, exit=${ec}), trying next..." | tee -a smoke.log
done

if [[ "${run_ok}" != "1" ]]; then
  echo "[smoke] ERROR: all MPI launchers failed"
  exit 1
fi

echo "[smoke] OK"
tail -n 25 smoke.log
