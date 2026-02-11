#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
cd "${ROOT}"

CASE_DIR="inputs/2DN/suites/physics_material/permutations_rz3" \
OUT_ROOT="outputs/suites/physics_material/permutations_rz3" \
./scripts/run_rz3_permutations.sh "$@"
