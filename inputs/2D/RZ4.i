# -----------------------------------------------------------------------------
# RZexp_ELM_nutrient_compare.i  (with YP3D-style debug aux fields + D_iso)
# -----------------------------------------------------------------------------
[GlobalParams]
  displacements = 'ux uy'
[]

[Mesh]
  coord_type    = RZ
  rz_coord_axis = Y
  parallel_type = distributed  # axial = y, radial = x  (so r = x, z = y)

  [gm]
    type = GeneratedMeshGenerator
    dim  = 2
    nx   = 100
    ny   = 60
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
  # Reference (IC) cell loading field used by CellGelMixtureOpt
  [phi_ref_ic]
    family = LAGRANGE
    order  = FIRST
  []

  # ---- YP3D-style debug fields (ParaView) ----
  [c11]
    family = MONOMIAL
    order  = CONSTANT
  []
  [volume_ratio]
    family = MONOMIAL
    order  = CONSTANT
  []
  [ke]
    family = MONOMIAL
    order  = CONSTANT
  []
  [kh]
    family = MONOMIAL
    order  = CONSTANT
  []
  [fa]
    family = MONOMIAL
    order  = CONSTANT
  []
  [press]
    family = MONOMIAL
    order  = CONSTANT
  []
  [gp]
    family = MONOMIAL
    order  = CONSTANT
  []
  [gate_tot]
    family = MONOMIAL
    order  = CONSTANT
  []
  [eta]
    family = MONOMIAL
    order  = CONSTANT
  []
  [phi_ref_aux]
    family = MONOMIAL
    order  = CONSTANT
  []
  [phi_cell]
    family = MONOMIAL
    order  = CONSTANT
  []

  # --- Diffusivity diagnostics (ParaView) ---
  # D_iso := (1/3) tr(D_eff)
  [Dxx]
    family = MONOMIAL
    order  = CONSTANT
  []
  [Dyy]
    family = MONOMIAL
    order  = CONSTANT
  []
  [Dzz]
    family = MONOMIAL
    order  = CONSTANT
  []
  [D_iso]
    family = MONOMIAL
    order  = CONSTANT
  []
[]

[UserObjects]
  [phi_ref_solution]
    type             = SolutionUserObject
    mesh             = phi_ref_filter_out2.e
    system_variables = 'phi_ref_ic'
    timestep         = LATEST
  []
[]

[ICs]
  active = 'phi_ref_from_file n0'

  [phi_ref_from_file]
    type          = SolutionIC
    variable      = phi_ref_ic
    solution_uo   = phi_ref_solution
    from_variable = phi_ref_ic
  []

  [phi_ref_ic_rand]
    type     = RandomIC
    variable = phi_ref_ic
    min      = 0.01
    max      = 0.20
    seed     = 123
  []

  [n0]
    type     = ConstantIC
    variable = n
    value    = 0.0
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

      [all]
      []
    []
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
    args          = n
  []
[]

[BCs]
  active = 'ux_axis uy_bottom n_top n_right'

  [ux_axis]
    type     = DirichletBC
    variable = ux
    boundary = left
    value    = 0.0
  []

  [uy_bottom]
    type     = DirichletBC
    variable = uy
    boundary = bottom
    value    = 0.0
  []

  [n_top]
    type     = DirichletBC
    variable = n
    boundary = top
    value    = 1.0
  []

  [n_right]
    type     = DirichletBC
    variable = n
    boundary = right
    value    = 1.0
  []
[]

[Materials]
  [cell_gel]
    type  = CellGelMixtureOpt
    block = 0

    # --- mechanics / growth ---
    G_cell    = 1.0
    q_cell    = -2.0

    k_exp_max = 0.1
    c1        = 10.0
    c2        = 2.0

    k_T1_max  = 0.0
    chi_str   = 1.0
    beta_T1   = 0.0
    m_T1      = 2.0

    # pressure gating ON
    press_str = 1.0

    G_gel     = 0.01
    k_diss_0  = 0.0

    # --- nutrient transport/consumption (coupled) ---
    D_nutrient = 0.1
    phi_cell_0 = 0.10

    phi_ref_ic = phi_ref_ic

    epsilon = 1e-6

    n_c1     = 0.05
    n_c2     = 12.0
    gamma_n0 = 20.0
    n        = n
  []
[]

[AuxKernels]
  # Stress component for ParaView
  [c11_aux]
    type            = RankTwoAux
    variable        = c11
    rank_two_tensor = cauchy_stress
    index_i         = 0
    index_j         = 0
  []

  # Material scalar diagnostics
  [vol_aux]
    type     = MaterialRealAux
    variable = volume_ratio
    property = volume_ratio
  []
  [ke_aux]
    type     = MaterialRealAux
    variable = ke
    property = ke
  []
  [kh_aux]
    type     = MaterialRealAux
    variable = kh
    property = kh
  []
  [fa_aux]
    type     = MaterialRealAux
    variable = fa
    property = fa
  []
  [press_aux]
    type     = MaterialRealAux
    variable = press
    property = pressure
  []
  [gp_aux]
    type     = MaterialRealAux
    variable = gp
    property = gp
  []
  [gate_tot_aux]
    type     = MaterialRealAux
    variable = gate_tot
    property = gate_tot
  []
  [eta_aux]
    type     = MaterialRealAux
    variable = eta
    property = eta
  []

  [phi_ref_aux_k]
    type       = MaterialRealAux
    variable   = phi_ref_aux
    property   = phi_cell_ref
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [phi_cell_aux]
    type       = MaterialRealAux
    variable   = phi_cell
    property   = phi_cell
    execute_on = 'INITIAL TIMESTEP_END'
  []

  # Diffusivity diagnostics: D_iso = tr(D_eff)/3
  [Dxx_aux]
    type      = MaterialRealTensorValueAux
    variable  = Dxx
    property  = D_eff
    component = 0
  []
  [Dyy_aux]
    type      = MaterialRealTensorValueAux
    variable  = Dyy
    property  = D_eff
    component = 4
  []
  [Dzz_aux]
    type      = MaterialRealTensorValueAux
    variable  = Dzz
    property  = D_eff
    component = 8
  []
  [D_iso_aux]
    type              = ParsedAux
    variable          = D_iso
    expression        = '(Dxx + Dyy + Dzz)/3'
    coupled_variables = 'Dxx Dyy Dzz'
  []
[]

[Postprocessors]
  [avg_J]
    type     = ElementAverageValue
    variable = volume_ratio
  []
  [vol_change_pct]
    type       = ParsedPostprocessor
    expression = '(avg_J - 1.0)*100.0'
    pp_names   = 'avg_J'
  []

  [avg_phi_cell]
    type     = ElementAverageValue
    variable = phi_cell
  []
  [avg_fa]
    type     = ElementAverageValue
    variable = fa
  []

  [avg_Diso]
    type     = ElementAverageValue
    variable = D_iso
  []
[]

[Preconditioning]
  [pc]
    type = SMP
    full = true
  []
[]

[Executioner]
  type        = Transient
  scheme      = bdf2
  solve_type  = PJFNK
  line_search = bt

  dtmin    = 1e-4
  dtmax    = 0.1
  end_time = 48.0

  nl_rel_tol = 1e-5
  nl_abs_tol = 1e-8
  nl_max_its = 100

  petsc_options_iname = '-snes_type -snes_linesearch_type -ksp_type -ksp_rtol -pc_type -pc_hypre_type -snes_rtol -snes_atol'
  petsc_options_value =  'newtonls bt gmres 1e-5 hypre boomeramg 1e-4 1e-7'

  automatic_scaling = true

  [TimeStepper]
    type               = IterationAdaptiveDT
    dt                 = 1e-2
    growth_factor      = 1.3
    cutback_factor     = 0.7
    optimal_iterations = 70
    iteration_window   = 2
  []
[]

[Outputs]
  exodus = true
  csv    = true
  perf_graph = true
[]
