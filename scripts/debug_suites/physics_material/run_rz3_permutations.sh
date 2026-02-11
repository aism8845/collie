#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
cd "${ROOT}"
exec ./scripts/suites/physics_material/run_rz3_permutations.sh "$@"
