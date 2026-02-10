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
export MPICH_CH4_NETMOD="${MPICH_CH4_NETMOD:-shm}"
export MPIR_CVAR_CH4_NETMOD="${MPIR_CVAR_CH4_NETMOD:-shm}"
export FI_PROVIDER="${FI_PROVIDER:-shm}"
export FI_SHM_DISABLE_CMA="${FI_SHM_DISABLE_CMA:-1}"
export FI_SHM_ENABLE_XPMEM="${FI_SHM_ENABLE_XPMEM:-0}"
export MPICH_OFI_DISABLE_DATAGRAM="${MPICH_OFI_DISABLE_DATAGRAM:-1}"
if [[ -n "${SMOKE_NP:-}" ]]; then
  echo "[smoke] MPI enabled: SMOKE_NP=${SMOKE_NP}"
  command -v mpiexec
  mpiexec -version
  mpiexec -n "${SMOKE_NP}" "${EXE}" -i RZ_smoke.i > smoke.log 2>&1
else
  echo "[smoke] serial run (set SMOKE_NP for MPI)"
  "${EXE}" -i RZ_smoke.i > smoke.log 2>&1
fi

echo "[smoke] OK"
tail -n 25 smoke.log
