# Best-known stable solver stack for bath_calibration:
# - 3-way Schur fieldsplit preconditioning (FSP AMG path)
# - robust globalization with newtonls + bt
# - Jacobian/preconditioner lagging for stability/throughput
#
# Use with:
#   inputs/solver_overrides/common_quick.i      (short diagnostics), or
#   no common_quick override                     (full runs)
# Optionally add:
#   inputs/solver_overrides/quiet_console.i

Preconditioning/active := 'pc_fsp3_schur_amg'

Executioner/petsc_options := '-snes_converged_reason -ksp_converged_reason'
Executioner/petsc_options_iname := '-snes_type -snes_linesearch_type -ksp_type -ksp_rtol -ksp_max_it -pc_type -pc_fieldsplit_type -snes_lag_preconditioner -snes_lag_preconditioner_persists -snes_lag_jacobian -snes_lag_jacobian_persists'
Executioner/petsc_options_value := 'newtonls bt fgmres 1e-6 200 fieldsplit schur 2 true 2 true'

