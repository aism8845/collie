# -----------------------------------------------------------------------------
# bath48_window_repro.i
# -----------------------------------------------------------------------------
# Short restart-window reproducer for the late-time mechanics barrier.
# This file is intended to be used with --recover from a checkpoint near t~31.5.
#
# Notes:
# - Physics/materials/BCs are inherited unchanged from bath_calibration.i.
# - 48h wrapper behavior needed at this time window is re-applied here (V_bath and refill pulse).
# - This run keeps strict-ish stepping in a narrow window to reproduce failure fast.
!include bath_calibration.i

ScalarKernels/bath_mb/V_bath := 5000
ScalarKernels/bath_mb/n_feed := 1.0
ScalarKernels/bath_mb/k_refill := k_refill_pulse

Executioner/end_time := 31.60
Executioner/dtmax := 0.01
Executioner/dtmin := 1e-5
Executioner/TimeStepper/dt := 0.0025
Executioner/TimeStepper/growth_factor := 1.0
Executioner/TimeStepper/cutback_factor := 0.5
Executioner/TimeStepper/optimal_iterations := 20
Executioner/TimeStepper/iteration_window := 4

Outputs/solver_watch/file_base := outputs/window_repro_solver_watch
Outputs/mesh_watch/file_base := outputs/window_repro_mesh_watch
Outputs/chk/file_base := outputs/chk/window_repro

[Functions]
  [k_refill_pulse]
    type = ParsedFunction
    expression = '500.0*(0.5*(1+tanh((t-24.0)/0.02)) - 0.5*(1+tanh((t-24.1)/0.02)))'
  []
[]
