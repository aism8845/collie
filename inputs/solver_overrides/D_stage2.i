# Strategy D fallback (stage 2): restore target consumption and continue
Executioner/scheme := implicit-euler
Executioner/end_time := 0.25

Materials/cell_gel/gamma_n0 := 500.0
Materials/nutrient_tl/gamma_n0 := 500.0

Executioner/petsc_options := '-snes_converged_reason -ksp_converged_reason -snes_monitor'
Executioner/petsc_options_iname := '-snes_type -snes_linesearch_type -ksp_type -pc_type -pc_factor_mat_solver_type -snes_rtol -snes_atol'
Executioner/petsc_options_value := 'newtonls bt preonly lu mumps 1e-4 1e-7'

