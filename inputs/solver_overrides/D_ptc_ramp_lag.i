# Strategy D: pseudo-transient damping + continuation + lagging
Executioner/scheme := implicit-euler

[Functions]
  [gamma_ramp]
    type = ParsedFunction
    expression = '500.0*min(1,t/0.05)'
  []
[]

[Controls]
  [gamma_ctrl]
    type = RealFunctionControl
    parameter = 'Materials/*/gamma_n0'
    function = gamma_ramp
    execute_on = 'initial timestep_begin'
  []
[]

Executioner/petsc_options := '-snes_converged_reason -ksp_converged_reason -snes_monitor'
Executioner/petsc_options_iname := '-snes_type -ksp_type -pc_type -pc_factor_mat_solver_type -snes_lag_preconditioner -snes_lag_preconditioner_persists -snes_lag_jacobian -snes_lag_jacobian_persists -snes_rtol -snes_atol'
Executioner/petsc_options_value := 'newtontr preonly lu mumps 2 true 2 true 1e-4 1e-7'

