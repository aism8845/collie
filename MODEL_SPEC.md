# Collie MODEL_SPEC

This document is the engineering/physics contract for this repository.
Any code changes (human or AI) must satisfy the constraints below.

## 1) Non-negotiables

### 1.1 SolidMechanics workflow
- Use SolidMechanics *new system* only (no TensorMechanics objects).
- Finite strain. Use displaced mesh where appropriate.
- Displacement variable naming must match the input files used:
  - 3D: ux uy uz
  - RZ: u_r u_z
- Tractions/pressures must use SolidMechanics equivalents (avoid legacy objects that are not registered).

### 1.2 AD everywhere (no type-mixing)
- New kernels/materials must be AD (ADKernel/ADMaterial).
- A material property name must exist in **exactly one type** across the entire app:
  - If a property is ever needed in an AD context, it is AD everywhere.
  - Never retrieve/declare the same property name as both AD and non-AD.
  - Example failure mode we must avoid:
    "requested AD property 'volume_ratio' is already retrieved/declared as non-AD" (or vice versa).

### 1.3 Property naming conventions
- Use lower_snake_case for property names.
- Use suffixes to encode meaning (do not create duplicate names for the same quantity):
  - _old: old/state value (only if MOOSE provides it; don’t hand-roll)
  - _ref: reference/initial value
  - _trial: local variable only (do not store unless explicitly required)
- Keep units consistent and documented in comments where ambiguity exists.

### 1.4 Tensor types and includes
- Prefer ADRankTwoTensor (not ADSymmetricRankTwoTensor) unless symmetry is guaranteed by construction.
- Include the correct forward headers for AD tensor types when needed.

### 1.5 No “creative” physics edits
- Do not change constitutive equations, kinematics, or sign conventions unless explicitly requested.
- If a refactor changes algebraic structure, add a comment explaining equivalence and ensure tests pass.

## 2) Stateful variable update rules (MOOSE Materials)

MOOSE calls computeQpProperties() repeatedly during nonlinear iterations.
Stateful history must be treated carefully.

### 2.1 Declaring stateful properties
- Any history variable must be declared as stateful (e.g., declareADProperty<Real>("x") with stateful flag if applicable in your pattern).
- Old values must be read via the *_old material property handle that MOOSE provides for stateful props.

### 2.2 Writing rules (critical)
- Only write “new” values into the current qp slot (e.g., _x[_qp]).
- Never overwrite old values directly.
- Do not mutate stateful history in a way that depends on unconverged Newton iterations.

### 2.3 Trial/iteration-dependent quantities
- If you need iteration-local “trial” values, keep them as local variables in computeQpProperties()
  or store them in non-stateful properties explicitly labeled as iteration-local.

### 2.4 Consistency
- Residual/Jacobian consistency is mandatory:
  - If a quantity enters the residual, it must be AD-tracked (or have a provably consistent tangent).
  - Avoid hidden branching or discontinuities without smoothing/regularization.

## 3) Required tests before committing

Run these from the repo root.

### 3.1 Build
- make -j8

### 3.2 Smoke run (fastest representative input)
- ./scripts/smoke_test.sh

### 3.3 App test harness (when available)
- ./run_tests -j8

### 3.4 Jacobian check (when debugging Newton stalls)
- If supported by the executable/build, run the Jacobian checker on a minimal input and ensure no large discrepancies.
  (Use a small mesh / short run to keep it fast.)

## 4) AI-agent operating rules

If an AI agent is used:
- It must follow this MODEL_SPEC strictly.
- It must make minimal diffs, explain why each file changed, and never “invent” new physics.
- It must run the required tests and report pass/fail with the relevant log tail on failure.

## 3) Required tests before committing

Run these from the repo root before any substantial commit (or before an AI agent opens a PR):

```bash
make -j8
./scripts/smoke_test.sh
./run_tests -j8   # optional: only if you have harness tests configured

