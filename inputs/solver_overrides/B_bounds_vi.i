# Strategy B: positivity bound using PETSc SNES-VI
[AuxVariables]
  [bounds_dummy]
    # Bounds aux must match the bounded variable FE type.
    family = LAGRANGE
    order = FIRST
  []
[]

[Bounds]
  [n_lower]
    type = ConstantBounds
    variable = bounds_dummy
    bounded_variable = n
    bound_type = lower
    bound_value = 0
    execute_on = 'initial timestep_begin'
  []
[]

Executioner/petsc_options := '-snes_converged_reason -ksp_converged_reason -snes_vi_monitor'
Executioner/petsc_options_iname := '-snes_type -ksp_type -pc_type -pc_factor_mat_solver_type -snes_rtol -snes_atol'
Executioner/petsc_options_value := 'vinewtonrsls preonly lu mumps 1e-4 1e-7'
