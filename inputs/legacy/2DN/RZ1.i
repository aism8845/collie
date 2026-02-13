# -----------------------------------------------------------------------------
# RZ_debug_coupled_min_FAILFAST_MEANINGFUL.i
# -----------------------------------------------------------------------------

[GlobalParams]
  displacements       = 'ux uy ur uz'
  use_displaced_mesh  = false
[]

[Mesh]
  coord_type    = RZ
  rz_coord_axis = Y
  parallel_type = replicated

  [gm]
    type = GeneratedMeshGenerator
    dim  = 2
    nx   = 40
    ny   = 12
    xmin = 0.0
    xmax = 2.0
    ymin = 0.0
    ymax = 0.5
  []
[]

[Variables]
  [n]
    family = LAGRANGE
    order  = FIRST
  []
[]

[AuxVariables]
  [phi_ref_ic]   family = LAGRANGE order = FIRST []

  [volume_ratio] family = MONOMIAL order = CONSTANT []
  [press]        family = MONOMIAL order = CONSTANT []
  [fa]           family = MONOMIAL order = CONSTANT []
  [phi_cell]     family = MONOMIAL order = CONSTANT []
  [D_phys_out]   family = MONOMIAL order = CONSTANT []
  [D_ref_out]    family = MONOMIAL order = CONSTANT []
[]

[ICs]
  active = 'phi_ref_rand n0'

  [phi_ref_rand]
    type     = RandomIC
    variable = phi_ref_ic
    min      = 0.01
    max      = 0.20
    seed     = 123
  []

  [n0]
    type     = ConstantIC
    variable = n
    value    = 1.0     # <-- DEBUG: start "on" to activate fa(n) immediately
  []
[]

[Physics]
  [SolidMechanics]
    [QuasiStatic]
      displacements = 'ux uy'
      add_variables = true
      strain        = FINITE
      new_system    = true
      formulation   = TOTAL
      [all] []
    []
  []
[]

[Kernels]
  [n_td]   type = TimeDerivative     variable = n []
  [n_diff] type = MatAnisoDiffusion  variable = n diffusivity = D_eff []
  [n_rxn]
    type          = MatReaction
    variable      = n
    reaction_rate = n_source_ref
    args          = n
  []
[]

[BCs]
  # Added ux_right to force nontrivial mechanics response early
  active = 'ux_axis uy_bottom ux_right n_top n_right'

  [ux_axis]    type = DirichletBC variable = ux boundary = left   value = 0.0 []
  [uy_bottom]  type = DirichletBC variable = uy boundary = bottom value = 0.0 []

  # DEBUG: radial confinement at outer boundary (forces stress/pressure when growth acts)
  [ux_right]   type = DirichletBC variable = ux boundary = right  value = 0.0 []

  [n_top]      type = DirichletBC variable = n  boundary = top    value = 1.0 []
  [n_right]    type = DirichletBC variable = n  boundary = right  value = 1.0 []
[]

[Materials]

    [nutrient_tl]
    type = ADNutrientTLTransport
    block = 0

    n      = n
    disp_r = ux
    disp_z = uy

    axisymmetric = true
    radial_coord = 0
    r_eps        = 1e-12

    phi_cell_ref = phi_cell_ref

    D0      = ${D0}
    D_floor = ${D_floor}

    gamma_n0 = ${gamma_n0}
    phi_max  = ${phi_max}

    n_c1 = ${n_c1}
    n_c2 = 12          # <-- NOTE: integer (not 12.0)
  []


  [cell_gel]
    type  = CellGelMixtureOpt
    block = 0

    # mechanics / growth
    G_cell    = 1.0
    q_cell    = -2.0
    k_exp_max = 0.2
    c1        = 10.0
    c2        = 2.0

    k_T1_max  = 0.0
    chi_str   = 1.0
    beta_T1   = 0.0
    m_T1      = 2.0

    press_str = 1.0
    G_gel     = 0.01
    k_diss_0  = 0.0

    # nutrient
    D_nutrient = 0.1
    phi_cell_0 = 0.10
    phi_ref_ic = phi_ref_ic
    epsilon    = 1e-6

    n_c1     = 0.05
    n_c2     = 12.0
    gamma_n0 = 20.0
    n        = n
  []
[]

[AuxKernels]
  active = 'vol_aux press_aux fa_aux phi_cell_aux D_phys_aux D_ref_aux'

  # DEBUG: execute all material-to-aux projections at INITIAL so t=0 is meaningful
  [vol_aux]
    type       = MaterialRealAux
    variable   = volume_ratio
    property   = volume_ratio
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [press_aux]
    type       = MaterialRealAux
    variable   = press
    property   = pressure
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [fa_aux]
    type       = MaterialRealAux
    variable   = fa
    property   = fa
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [phi_cell_aux]
    type       = MaterialRealAux
    variable   = phi_cell
    property   = phi_cell
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [D_phys_aux]
    type       = MaterialRealAux
    variable   = D_phys_out
    property   = D_phys
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [D_ref_aux]
    type       = MaterialRealAux
    variable   = D_ref_out
    property   = D_ref
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]

[Postprocessors]
  active = 'dt_pp avg_J J_min J_max press_min press_max n_min n_max'

  # DEBUG: execute on INITIAL so the "Time Step 0" table is not all zeros
  [dt_pp]
    type       = TimestepSize
    execute_on = 'INITIAL TIMESTEP_END'
  []

  [avg_J]
    type       = ElementAverageValue
    variable   = volume_ratio
    execute_on = 'INITIAL TIMESTEP_END'
  []

  [J_min]
    type       = ElementExtremeValue
    variable   = volume_ratio
    value_type = min
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [J_max]
    type       = ElementExtremeValue
    variable   = volume_ratio
    value_type = max
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [press_min]
    type       = ElementExtremeValue
    variable   = press
    value_type = min
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [press_max]
    type       = ElementExtremeValue
    variable   = press
    value_type = max
    execute_on = 'INITIAL TIMESTEP_END'
  []

  [n_min]
    type       = NodalExtremeValue
    variable   = n
    value_type = min
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [n_max]
    type       = NodalExtremeValue
    variable   = n
    value_type = max
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]

[Debug]
  show_var_residual_norms = true
[]

[Executioner]
  type        = Transient
  scheme      = bdf2
  solve_type  = NEWTON
  line_search = bt

  nl_rel_tol = 1e-6
  nl_abs_tol = 1e-9
  nl_max_its = 80

  automatic_scaling = false

  petsc_options       = '-snes_test_jacobian_view -snes_converged_reason'
  petsc_options_iname = '-snes_test_jacobian -ksp_type -pc_type'
  petsc_options_value = '1e-8 preonly lu'

  [TimeStepper]
    type = ConstantDT
    dt   = 1e-2
  []
  end_time = 0.2   # <-- DEBUG: keep short while iterating; raise later
[]

[Outputs]
  exodus     = true
  csv        = true
  perf_graph = true

  [dofmap]
    type            = DOFMap
    execute_on      = 'INITIAL'
    file_base       = 'RZ1'
    file_base_suffix = '_dofmap'
  []
[]