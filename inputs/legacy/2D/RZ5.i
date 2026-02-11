# -----------------------------------------------------------------------------
# RZ_debug_coupled_min.i
# Minimal coupled mechanics + nutrient for Jacobian/tangent debugging.
#
# TOGGLES (edit these lines):
#   1) Solver:        set solve_type = NEWTON  OR  PJFNK
#   2) TimeStepper:   use ConstantDT (recommended for debug) or IterationAdaptiveDT
#   3) IC source:     [ICs] active = 'phi_ref_rand n0'  OR  'phi_ref_from_file n0'
#   4) Debug outputs: [AuxVariables]/[AuxKernels]/[Postprocessors] active = '...'
#   5) Jacobian test: add PETSc options -snes_test_jacobian ... (see bottom)
# -----------------------------------------------------------------------------

[GlobalParams]
  displacements = 'ux uy'
  use_displaced_mesh = false
[]

[Mesh]
  coord_type    = RZ
  rz_coord_axis = Y
  parallel_type = replicated

  [gm]
    type = GeneratedMeshGenerator
    dim  = 2

    # keep small for jacobian tests; scale up later
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
  # reference field used by your material (IC loaded or randomized)
  [phi_ref_ic]
    family = LAGRANGE
    order  = FIRST
  []

  # ----------------- DEBUG AUX (toggle via AuxKernels active list) -----------------
  [volume_ratio]  family = MONOMIAL  order = CONSTANT []
  [press]         family = MONOMIAL  order = CONSTANT []
  [fa]            family = MONOMIAL  order = CONSTANT []
  [gp]            family = MONOMIAL  order = CONSTANT []
  [gate_tot]      family = MONOMIAL  order = CONSTANT []
  [phi_cell]      family = MONOMIAL  order = CONSTANT []
  [D_phys_out]    family = MONOMIAL  order = CONSTANT []
  [D_ref_out]     family = MONOMIAL  order = CONSTANT []
[]

# OPTIONAL: only needed if you truly want phi_ref from an exodus solution
[UserObjects]
  [phi_ref_solution]
    type             = SolutionUserObject
    mesh             = phi_ref_filter_out2.e
    system_variables = 'phi_ref_ic'
    timestep         = LATEST
  []
[]

[ICs]
  # TOGGLE THIS:
  # active = 'phi_ref_from_file n0'
  active = 'phi_ref_rand n0'

  [phi_ref_from_file]
    type          = SolutionIC
    variable      = phi_ref_ic
    solution_uo   = phi_ref_solution
    from_variable = phi_ref_ic
  []

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
    value    = 0.5
  []
[]

[Physics]
  [SolidMechanics]
    [QuasiStatic]
      displacements      = 'ux uy'
      add_variables      = true
      strain             = FINITE
      new_system         = true
      formulation        = TOTAL
      [all] []
    []
  []
[]

[Kernels]
  # nutrient PDE
  [n_td]   type = TimeDerivative   variable = n []
  [n_diff] type = MatAnisoDiffusion variable = n diffusivity = D_eff []
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

    # --- nutrient transport/consumption ---
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
  # TOGGLE: comment/uncomment this active list for debug fields
  active = 'vol_aux press_aux fa_aux gp_aux gate_tot_aux phi_cell_aux D_phys_aux D_ref_aux'

  [vol_aux]      type = MaterialRealAux variable = volume_ratio property = volume_ratio []
  [press_aux]    type = MaterialRealAux variable = press        property = pressure     []
  [fa_aux]       type = MaterialRealAux variable = fa           property = fa           []
  [gp_aux]       type = MaterialRealAux variable = gp           property = gp           []
  [gate_tot_aux] type = MaterialRealAux variable = gate_tot     property = gate_tot     []
  [phi_cell_aux] type = MaterialRealAux variable = phi_cell     property = phi_cell     execute_on='INITIAL TIMESTEP_END' []

  [D_phys_aux]   type = MaterialRealAux variable = D_phys_out   property = D_phys []
  [D_ref_aux]    type = MaterialRealAux variable = D_ref_out    property = D_ref  []
[]

[Postprocessors]
  # TOGGLE: turn these on/off (they’re cheap, but keep minimal)
  active = 'avg_J vol_change_pct avg_phi_cell avg_fa'

  [avg_J]
    type     = ElementAverageValue
    variable = volume_ratio
  []
  [vol_change_pct]
    type       = ParsedPostprocessor
    expression = '(avg_J - 1.0)*100.0'
    pp_names   = 'avg_J'
  []
  [avg_phi_cell] type = ElementAverageValue variable = phi_cell []
  [avg_fa]       type = ElementAverageValue variable = fa       []
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

  # TOGGLE SOLVER:
  # solve_type  = PJFNK
  solve_type  = NEWTON

  line_search = bt

  dtmin    = 1e-4
  dtmax    = 1e-2
  end_time = 10.0   # keep short for debug; raise once stable

  nl_rel_tol = 1e-6
  nl_abs_tol = 1e-9
  nl_max_its = 50

  automatic_scaling = false   # recommended while doing jacobian checks

  # ----------------- PETSc options (choose one set) -----------------
  # Robust direct solve for small problems (best for Newton + jacobian debug)
  petsc_options_iname = '-snes_type -snes_linesearch_type -ksp_type -pc_type -snes_rtol -snes_atol'
  petsc_options_value =  'newtonls bt preonly lu 1e-6 1e-9'

  # If you want PJFNK:
  # petsc_options_iname = '-snes_type -snes_linesearch_type -ksp_type -ksp_rtol -pc_type -snes_rtol -snes_atol'
  # petsc_options_value =  'newtonls bt gmres 1e-6 hypre 1e-6 1e-9'

  # ----------------- TIME STEPPER (choose one) -----------------
  [TimeStepper]
    type = ConstantDT
    dt   = 1e-2
  []

  # [TimeStepper]
  #   type               = IterationAdaptiveDT
  #   dt                 = 1e-2
  #   growth_factor      = 1.3
  #   cutback_factor     = 0.7
  #   optimal_iterations = 12
  #   iteration_window   = 2
  # []
[]

[Outputs]
  exodus     = true
  csv        = true
  perf_graph = true
[]
