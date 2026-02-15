# Reproducible solver stack used for robust bath24_best-style runs.
# Solver-only knobs (no physics changes).

Preconditioning/active := 'pc_fsp3_schur_amg'

# Prevent over-solving in the stiff late window (keeps Newton off invalid mechanics states).
Executioner/nl_abs_tol := 1e-4
Executioner/nl_rel_tol := 1e-5
Executioner/nl_max_its := 120

Executioner/petsc_options := '-snes_converged_reason -ksp_converged_reason'
Executioner/petsc_options_iname := '-snes_type -snes_linesearch_type -ksp_type -ksp_rtol -ksp_max_it -pc_type -pc_fieldsplit_type -snes_lag_preconditioner -snes_lag_preconditioner_persists -snes_lag_jacobian -snes_lag_jacobian_persists'
Executioner/petsc_options_value := 'newtonls bt fgmres 1e-6 200 fieldsplit schur 2 true 2 true'
