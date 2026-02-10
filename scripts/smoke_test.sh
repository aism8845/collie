#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
EXE="${ROOT}/collie-opt"
INP="${ROOT}/inputs/2DN/RZ_smoke.i"

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
SMOKE_NP="${SMOKE_NP:-8}"
echo "[smoke] MPI run: SMOKE_NP=${SMOKE_NP}"
command -v mpiexec
mpiexec -n "${SMOKE_NP}" "${EXE}" -i RZ_smoke.i > smoke.log 2>&1

echo "[smoke] OK"
tail -n 25 smoke.log
