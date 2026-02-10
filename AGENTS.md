# Agent Rules for collie (MOOSE / SolidMechanics / AD)

## Non-negotiables
- SolidMechanics new system only (no TensorMechanics).
- AD everywhere: never mix AD and non-AD material properties with the same name.
- Do not change constitutive equations, sign conventions, or kinematics unless explicitly instructed.

## How to work
- Make minimal diffs.
- Prefer plumbing/refactors (params, member declarations, registration, includes) over physics edits.
- Explain what changed and why.

## Required checks before declaring success
- Run: make -j8
- Run: ./scripts/smoke_test.sh
- Smoke defaults to MPI via `mpiexec` with `SMOKE_NP=8` (override using `SMOKE_NP`).
- Canonical quick debug smoke command: `./scripts/smoke_debug_mpi8.sh`
- If a check fails: stop, report the error, propose the smallest fix.
