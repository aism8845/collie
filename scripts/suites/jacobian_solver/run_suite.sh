#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
cd "${ROOT}"

CASE_DIR="inputs/legacy/2DN/suites/jacobian_solver/permutations_rz3" \
OUT_ROOT="outputs/suites/jacobian_solver/permutations_rz3" \
./scripts/run_rz3_permutations.sh "$@"
