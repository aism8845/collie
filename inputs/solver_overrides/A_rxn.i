# Strategy A (IMEX split): reaction/mechanics full-step with subcycling
[Kernels]
  active := 'u_n_offdiag_x u_n_offdiag_y n_td n_rxn'
[]

Executioner/scheme := implicit-euler
Executioner/TimeStepper/type := ConstantDT
Executioner/TimeStepper/dt := 1e-3
# Executioner/end_time is provided by the driver at each split stage.

Executioner/petsc_options := '-snes_converged_reason -ksp_converged_reason'
Executioner/petsc_options_iname := '-snes_type -snes_linesearch_type -ksp_type -ksp_rtol -pc_type -pc_factor_mat_solver_type -snes_rtol -snes_atol'
Executioner/petsc_options_value := 'newtonls bt fgmres 1e-6 lu mumps 1e-4 1e-7'
