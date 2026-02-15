#!/usr/bin/env bash
set -euo pipefail

# Usage:
#   ./dump_and_map_empty_rows.sh /abs/or/rel/path/to/app-opt /abs/or/rel/path/to/case.i [run_dir]
#
APP_RAW="${1:?Need app binary path, e.g. ./collieTestApp-opt}"
INP_RAW="${2:?Need input file, e.g. ./inputs/current/RZ3_RD_AD_patch.i}"
RUNDIR="${3:-_zero_row_debug}"

# Resolve to absolute paths BEFORE we cd anywhere
APP="$(realpath "$APP_RAW")"
INP="$(realpath "$INP_RAW")"

if [[ ! -x "$APP" ]]; then
  echo "ERROR: App binary not found or not executable: $APP"
  echo "Tip: from repo root, run: ls -1 ./*-opt"
  exit 2
fi
if [[ ! -f "$INP" ]]; then
  echo "ERROR: Input file not found: $INP"
  exit 2
fi

mkdir -p "$RUNDIR"
cp -f "$INP" "$RUNDIR"/
INP_BASENAME="$(basename "$INP")"
cd "$RUNDIR"

# PETSc options file (passed through MOOSE/libMesh)
cat > petsc_opts.txt <<'EOF'
-snes_monitor
-snes_converged_reason
-ksp_monitor_short
-ksp_converged_reason

# Lightweight matrix summary in stdout
-mat_view ::ascii_info

# Dump matrix to PETSc binary for post-processing
-mat_view binary:J.petsc
EOF

echo "=== Running (serial recommended for easiest DOF mapping) ==="
echo "App:   $APP"
echo "Input: $INP_BASENAME"
echo "Dir:   $(pwd)"
echo

# Run (keep || true so we still try to map if the solve fails)
"$APP" -i "$INP_BASENAME" -options_file petsc_opts.txt || true

# Find a DOFMapOutput JSON:
# DOFMapOutput writes <file_base>.json (and maybe _1.json, _2.json...) :contentReference[oaicite:1]{index=1}
DOFMAP_JSON=""
# 1) Prefer anything that already has "dofmap" in the name (if you set a custom suffix)
DOFMAP_JSON="$(ls -1 *dofmap*.json 2>/dev/null | head -n 1 || true)"

# 2) Otherwise, scan *.json files and pick the first one that looks like DOFMapOutput schema
if [[ -z "$DOFMAP_JSON" ]]; then
  DOFMAP_JSON="$(python3 - <<'PY' || true
import glob, json, sys
for f in sorted(glob.glob("*.json")):
    try:
        d = json.load(open(f, "r"))
        if isinstance(d, dict) and "ndof" in d and "vars" in d:
            print(f)
            sys.exit(0)
    except Exception:
        pass
sys.exit(1)
PY
)"
fi

if [[ -z "$DOFMAP_JSON" ]]; then
  echo "ERROR: Could not find a DOFMapOutput JSON in $(pwd)."
  echo "JSON files present:"
  ls -1 *.json 2>/dev/null || true
  echo
  echo "Make sure your input includes a DOFMap output. Example:"
  echo "[Outputs]"
  echo "  [dofmap]"
  echo "    type = DOFMap"
  echo "    execute_on = 'INITIAL'"
  echo "  []"
  echo "[]"
  exit 2
fi

echo "DOF map JSON: $DOFMAP_JSON"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PY="${SCRIPT_DIR}/map_empty_rows_to_vars.py"

python3 "$PY" --mat J.petsc --dofmap "$DOFMAP_JSON" --out empty_rows_report.txt

echo
echo "Wrote: $(pwd)/empty_rows_report.txt"
