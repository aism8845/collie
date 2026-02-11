# 2Dexp1_ideal_AB_growth.i
# Goal: nutrient-gated growth + uptake (Loops A/B), with:
#   - NO pressure inhibition (press_str huge)
#   - NO T1 rearrangements (k_T1_max = 0)
#   - supply only at TOP (right boundary no longer Dirichlet)

# Same as nutrient-only file above, except k_exp_max is restored.

[GlobalParams]
  displacements    = 'ux uy'
  large_kinematics = true
[]

[Mesh]
  parallel_type = distributed
  [rect]
    type = GeneratedMeshGenerator
    dim  = 2
    nx   = 200
    ny   = 100
    xmin = 0.0
    xmax = 2.0
    ymin = 0.0
    ymax = 0.5
  []
[]

[UserObjects]
  # Reads phi_ref_ic from the Exodus output of phi_ref_filter.i
  [phi_ref_solution]
    type = SolutionUserObject
    mesh = phi_ref_filter_out2.e
    system_variables = 'phi_ref_ic'
    timestep = LATEST
  []
[]

[Variables]
  [n]
    family = LAGRANGE
    order  = FIRST
  []
[]

[Physics]
  [SolidMechanics]
    [QuasiStatic]
      [all]
        add_variables = true
        displacements = 'ux uy'
        strain      = FINITE
        new_system  = true
        formulation = TOTAL
      []
    []
  []
[]

[BCs]
  [ux_left]
    type = DirichletBC
    variable = ux
    boundary = left
    value = 0.0
  []
  [uy_bottom]
    type = DirichletBC
    variable = uy
    boundary = bottom
    value = 0.0
  []
  [ux_right]
    type = DirichletBC
    variable = ux
    boundary = right
    value = 0.0
  []

  [n_top]
    type = DirichletBC
    variable = n
    boundary = top
    value = 1.0
  []
  [n_bottom]
    type = NeumannBC
    variable = n
    boundary = bottom
    value = 0.0
  []
  [n_left]
    type = NeumannBC
    variable = n
    boundary = left
    value = 0.0
  []
  [n_right]
    type = NeumannBC
    variable = n
    boundary = right
    value = 0.0
  []
[]

[AuxVariables]
  [phi_ref_ic]
    family = LAGRANGE
    order  = FIRST
  []
  [fa]
    family = MONOMIAL
    order  = CONSTANT
  []
  [phi_cell_aux]
    family = MONOMIAL
    order  = CONSTANT
  []
  [ke]
    family = MONOMIAL
    order  = CONSTANT
  []
[]

[Materials]
  [cell_gel]
    type = CellGelMixtureOpt

    G_cell = 1.0
    q_cell = -2.0
    G_gel  = 0.01

    # --- IDEALIZATION toggles ---
    k_exp_max = 0.7        # <-- GROWTH ON 
    k_T1_max  = 0.0        # <-- NO rearrangements
    press_str = 1e12       # <-- remove pressure inhibition 

    # growth envelope (keep your values)
    c1 = 10.0
    c2 = 2.0

    # other params kept for completeness
    chi_str  = 0.3
    beta_T1  = 10.0
    m_T1     = 2.0

    k_diss_0   = 0.0
    D_nutrient = 0.1

    phi_cell_0 = 0.1
    phi_ref_ic = phi_ref_ic
    epsilon    = 1e-8

    # step-like nutrient gate
    n_c1 = 0.05
    n_c2 = 12.0

    gamma_n0 = 20
    n = n
  []
[]

[Kernels]
  [n_td]
    type     = TimeDerivative
    variable = n
  []
  [n_diff]
    type        = MatAnisoDiffusion
    variable    = n
    diffusivity = D_eff
  []
  [n_rxn]
    type          = MatReaction
    variable      = n
    reaction_rate = n_source_ref
    args          = 'n'
  []
[]

[AuxKernels]
  [fa_aux]
    type     = MaterialRealAux
    variable = fa
    property = fa
    execute_on = 'initial timestep_end'
  []
  [phi_cell_aux_k]
    type     = MaterialRealAux
    variable = phi_cell_aux
    property = phi_cell
    execute_on = 'initial timestep_end'
  []
  [ke_aux]
    type     = MaterialRealAux
    variable = ke
    property = ke
    execute_on = 'initial timestep_end'
  []
[]

[ICs]
  [phi_ref_uniform]
    type     = ConstantIC
    variable = phi_ref_ic
    value    = 0.1
  []
  [n_ic]
    type     = ConstantIC
    variable = n
    value    = 0.0
  []
[]

[Executioner]
  type        = Transient
  scheme      = bdf2
  solve_type  = NEWTON
  line_search = bt

  dtmin    = 1e-4
  dtmax    = 0.1
  end_time = 34.0

  nl_rel_tol = 1e-5
  nl_abs_tol = 1e-8
  nl_max_its = 100

  petsc_options_iname = '-snes_type -snes_linesearch_type -ksp_type -ksp_rtol -pc_type -pc_hypre_type -snes_rtol -snes_atol'
  petsc_options_value =  'newtonls bt gmres 1e-5 hypre boomeramg 1e-4 1e-7'

  [TimeStepper]
    type               = IterationAdaptiveDT
    dt                 = 1e-2
    growth_factor      = 1.3
    cutback_factor     = 0.7
    optimal_iterations = 25
    iteration_window   = 2
  []
[]

[Outputs]
  exodus = true
  csv    = true
[]
