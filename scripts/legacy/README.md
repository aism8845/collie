# Legacy Smoke Entrypoints

- `smoke_test_legacy.sh`: runs the canonical smoke harness against `inputs/legacy/2DN/RZ_smoke.i`.
- The wrapper disables `SMOKE_REQUIRE_CURRENT_INPUT` so legacy runs stay possible without changing current defaults.
- The active/default smoke entrypoint remains `scripts/smoke_test.sh`, which targets `inputs/current/RZ3_RD_AD_patch.i`.
