# Schur fieldsplit + robust globalization + Jacobian/PC lagging.
# NOTE: With PJFNK/mffd, SNES newtontr requires transpose ops that are unavailable.
# Use newtonls(bt) as the stable fallback globalization.
Preconditioning/active := 'pc_fsp3_schur_amg'

Executioner/petsc_options := '-snes_converged_reason -ksp_converged_reason'
Executioner/petsc_options_iname := '-snes_type -snes_linesearch_type -ksp_type -ksp_rtol -ksp_max_it -pc_type -pc_fieldsplit_type -snes_lag_preconditioner -snes_lag_preconditioner_persists -snes_lag_jacobian -snes_lag_jacobian_persists'
Executioner/petsc_options_value := 'newtonls bt fgmres 1e-6 200 fieldsplit schur 2 true 2 true'
