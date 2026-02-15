# Quick solver diagnostics (short horizon, modest dt envelope).
Executioner/end_time := 0.25

Executioner/dtmax := 0.02
Executioner/dtmin := 1e-4
Executioner/TimeStepper/dt := 1e-2
Executioner/TimeStepper/growth_factor := 1.2
Executioner/TimeStepper/cutback_factor := 0.5
Executioner/TimeStepper/optimal_iterations := 18

# Avoid clobbering calibration outputs.
Outputs/solver_watch/file_base := outputs/bench_solver_watch
Outputs/mesh_watch/file_base := outputs/bench_mesh_watch
