# -----------------------------------------------------------------------------
# RZ2_no_nutrient_isolation.i
# Purpose: remove nutrients completely (no n variable; no diffusion/reaction),
#          keep heterogeneity so we test whether mixture/growth/localization
#          alone drives stalling.
# -----------------------------------------------------------------------------
[GlobalParams]
  displacements = 'ux uy'
[]

[Mesh]
  coord_type    = RZ
  rz_coord_axis = Y
  parallel_type = distributed

  [gm]
    type = GeneratedMeshGenerator
    dim  = 2
    nx   = 80
    ny   = 40
    xmin = 0.0
    xmax = 2.0
    ymin = 0.0
    ymax = 0.5
  []
[]

[AuxVariables]
  [phi_ref_ic]
    family = LAGRANGE
    order  = FIRST
  []

  [J]
    family = MONOMIAL
    order  = CONSTANT
  []
  [phi_cell_out]
    family = MONOMIAL
    order  = CONSTANT
  []
  [fa_out]
    family = MONOMIAL
    order  = CONSTANT
  []
  [press_out]
    family = MONOMIAL
    order  = CONSTANT
  []
[]

# --- OPTIONAL: file-seeded heterogeneity (only if you have an exodus to read) ---
#[UserObjects]
#  [phi_ref_solution]
#    type             = SolutionUserObject
#    mesh             = phi_ref_filter_out2.e   # <-- CHANGE THIS if using
#    system_variables = 'phi_ref_ic'
#    timestep         = LATEST
#  []
#[]

[ICs]
  # Pick one:
  active = 'phi_ref_ic_rand'

  #[phi_ref_from_file]
  #  type          = SolutionIC
  #  variable      = phi_ref_ic
  #  solution_uo   = phi_ref_solution
  #  from_variable = phi_ref_ic
  #[]

  [phi_ref_ic_rand]
    type     = RandomIC
    variable = phi_ref_ic
    min      = 0.01
    max      = 0.20
    seed     = 123
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

[BCs]
  active = 'ux_axis uy_bottom'

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
[]

[Materials]
  [cell_gel]
    type  = CellGelMixtureOpt
    block = 0

    # mechanics/growth
    G_cell    = 1.0
    q_cell    = -2.0
    k_exp_max = 0.9
    c1        = 10.0
    c2        = 2.0
    press_str = 1.0

    k_T1_max  = 0.0
    chi_str   = 1.0
    beta_T1   = 0.0
    m_T1      = 2.0

    G_gel     = 0.0001
    k_diss_0  = 0.0

    phi_cell_0 = 0.4
    phi_ref_ic = phi_ref_ic

    # FD/central-diff tangent step (you’re moving to central diff)
    epsilon = 1e-6

    # --- IMPORTANT: no nutrient coupling here: do NOT set "n = n"
    # gamma_n0 can be left 0
    gamma_n0 = 0.0
  []
[]

[AuxKernels]
  [J_aux]
    type     = MaterialRealAux
    variable = J
    property = volume_ratio
  []
  [phi_cell_aux]
    type     = MaterialRealAux
    variable = phi_cell_out
    property = phi_cell
  []
  [fa_aux]
    type     = MaterialRealAux
    variable = fa_out
    property = fa
  []
  [press_aux]
    type     = MaterialRealAux
    variable = press_out
    property = pressure
  []
[]

[Postprocessors]
  [avg_J]
    type     = ElementAverageValue
    variable = J
  []
  [vol_change_pct]
    type       = ParsedPostprocessor
    expression = '(avg_J - 1.0)*100.0'
    pp_names   = 'avg_J'
  []
  [avg_phi_cell]
    type     = ElementAverageValue
    variable = phi_cell_out
  []
  [avg_fa]
    type     = ElementAverageValue
    variable = fa_out
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
  end_time = 34.0

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
