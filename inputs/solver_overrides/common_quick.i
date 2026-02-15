# Common short diagnostic window (no long runs).
# Bound work per case so benchmark wall time is predictable.
Executioner/end_time := 0.25
Executioner/num_steps := 5
Executioner/dtmin := 1e-4
Executioner/dtmax := 0.05

Executioner/TimeStepper/type := ConstantDT
Executioner/TimeStepper/dt := 0.05

# Keep quick diagnostics lightweight.
Outputs/exodus := false

# Coarse quick-benchmark mesh override.
Mesh/gm/nx := 60
Mesh/gm/ny := 20

# Avoid clobbering baseline calibration CSVs.
Outputs/solver_watch/file_base := outputs/bench_solver_watch
Outputs/mesh_watch/file_base := outputs/bench_mesh_watch
