# Recommended late-window rescue override for bath48 (solver-only).
# Use for a restart window when the baseline stalls near t~31.

Preconditioning/active := 'pc_fsp3_schur_amg'

Executioner/scheme := implicit-euler
Executioner/nl_max_its := 160
Executioner/nl_abs_tol := 1e-4
Executioner/nl_rel_tol := 1e-5
Executioner/dtmax := 0.005
Executioner/dtmin := 1e-4

Executioner/petsc_options := '-snes_converged_reason -ksp_converged_reason'
Executioner/petsc_options_iname := '-snes_type -snes_linesearch_type -ksp_type -ksp_rtol -ksp_max_it -ksp_gmres_restart -pc_type -pc_fieldsplit_type -snes_lag_preconditioner -snes_lag_preconditioner_persists -snes_lag_jacobian -snes_lag_jacobian_persists'
Executioner/petsc_options_value := 'newtonls bt fgmres 1e-5 300 200 fieldsplit schur 2 true 2 true'

Executioner/TimeStepper/type := IterationAdaptiveDT
Executioner/TimeStepper/dt := 0.0025
Executioner/TimeStepper/growth_factor := 1.0
Executioner/TimeStepper/cutback_factor := 0.5
Executioner/TimeStepper/optimal_iterations := 20
Executioner/TimeStepper/iteration_window := 4
