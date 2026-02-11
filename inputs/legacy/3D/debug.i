# YP3D3_cube_debug.i
# Minimal 3D debug case for CellGelMixtureOpt with nutrient + pressure gates.
# - Simple cube mesh (GeneratedMeshGenerator)
# - Minimal rigid-body pinning using 3 nodesets
# - Toggle gates/consumption by switching the active CellGelMixtureOpt block in [Materials]
#
# Time units: HOURS

[GlobalParams]
  displacements    = 'ux uy uz'
  large_kinematics = true
[]

[Mesh]
  # Replicated is simplest for single-core VS Code debugging.
  parallel_type = replicated

  [cube]
    type = GeneratedMeshGenerator
    dim  = 3
    nx   = 10
    ny   = 10
    nz   = 10
    xmin = 0.0
    xmax = 1.0
    ymin = 0.0
    ymax = 1.0
    zmin = 0.0
    zmax = 1.0
  []

  # --- Minimal rigid-body constraints via 3 pinned nodes (remove 6 rigid-body modes) ---
  # pin0: fixes ux,uy,uz at (0,0,0)
  # pin1: fixes uy,uz at (1,0,0)
  # pin2: fixes uz    at (0,1,0)
  [pin0]
    type             = ExtraNodesetGenerator
    input            = cube
    coord            = '0 0 0'
    use_closest_node = true
    new_boundary     = 'pin0'
  []
  [pin1]
    type             = ExtraNodesetGenerator
    input            = pin0
    coord            = '1 0 0'
    use_closest_node = true
    new_boundary     = 'pin1'
  []
  [pin2]
    type             = ExtraNodesetGenerator
    input            = pin1
    coord            = '0 1 0'
    use_closest_node = true
    new_boundary     = 'pin2'
  []

  final_generator = pin2
[]

[Variables]
  [n]
    family = LAGRANGE
    order  = FIRST
  []
[]

[AuxVariables]
  # Quenched heterogeneity field (optional)
  [phi_ref_ic]
    family = LAGRANGE
    order  = FIRST
  []

  # Diagnostics exported to Exodus
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
  [c11]
    family = MONOMIAL
    order  = CONSTANT
  []
[]

[Physics]
  [SolidMechanics]
    [QuasiStatic]
      [all]
        add_variables = true
        displacements = 'ux uy uz'
        strain        = FINITE
        new_system    = true
        formulation   = TOTAL
        use_displaced_mesh = true
      []
    []
  []
[]

[BCs]
  # --- Mechanical pins ---
  [pin0_ux]
    type     = DirichletBC
    preset   = true
    variable = ux
    boundary = pin0
    value    = 0.0
  []
  [pin0_uy]
    type     = DirichletBC
    preset   = true
    variable = uy
    boundary = pin0
    value    = 0.0
  []
  [pin0_uz]
    type     = DirichletBC
    preset   = true
    variable = uz
    boundary = pin0
    value    = 0.0
  []

  [pin1_uy]
    type     = DirichletBC
    preset   = true
    variable = uy
    boundary = pin1
    value    = 0.0
  []
  [pin1_uz]
    type     = DirichletBC
    preset   = true
    variable = uz
    boundary = pin1
    value    = 0.0
  []

  [pin2_uz]
    type     = DirichletBC
    preset   = true
    variable = uz
    boundary = pin2
    value    = 0.0
  []

  # --- Nutrient bath BCs on all outer faces ---
  [n_left]
    type     = DirichletBC
    boundary = left
    variable = n
    value    = 1.0
  []
  [n_right]
    type     = DirichletBC
    boundary = right
    variable = n
    value    = 1.0
  []
  [n_front]
    type     = DirichletBC
    boundary = front
    variable = n
    value    = 1.0
  []
  [n_back]
    type     = DirichletBC
    boundary = back
    variable = n
    value    = 1.0
  []
  [n_bottom]
    type     = DirichletBC
    boundary = bottom
    variable = n
    value    = 1.0
  []
  [n_top]
    type     = DirichletBC
    boundary = top
    variable = n
    value    = 1.0
  []
[]

# --- Optional: read a stored speckle field from an Exodus file ---
# If you want to use this, uncomment [UserObjects] + [AuxKernels][phi_ref_read],
# and remove/disable the RandomIC for phi_ref_ic in [ICs].
#
[UserObjects]
  [phi_ref_solution]
    type             = SolutionUserObject
    mesh             = phi_ref_filter3d_out.e
    system_variables = 'phi_ref_ic'
    timestep         = LATEST
    epsilon          = 1e-10
  []
[]

[Materials]
  # -----------------------------
  # TOGGLES (edit this one line):
  #   active = 'cell_gel_full'
  #   active = 'cell_gel_no_pressure_gate'
  #   active = 'cell_gel_no_nutrient_gate'
  #   active = 'cell_gel_no_consumption'
  #   active = 'cell_gel_mech_only'   # leaves n PDE active, but removes all gates + consumption
  # -----------------------------
  active = 'cell_gel_full'

  [cell_gel_full]
    type = CellGelMixtureOpt

    # cell phase
    G_cell    = 1.0
    q_cell    = -2.0

    # expansion (slow/stable for debugging)
    k_exp_max = .1
    c1        = 2.0
    c2        = 2.0

    # pressure gate (ON)
    press_str = 1.0

    # T1 (set k_T1_max=0 to disable)
    k_T1_max  = 0.5
    chi_str   = 0.3
    beta_T1   = 10.0
    m_T1      = 0.0

    # gel phase
    G_gel    = 0.01
    k_diss_0 = 0.05

    # nutrient coupling + activity gate (ON)
    n        = n
    D_nutrient = 0.1
    n_c1     = 0.05
    n_c2     = 2.0
    n_eps    = 1e-8

    # consumption (ON): sink is -J * gamma_n_local * n_+
    gamma_n0 = 20.0

    # initial reference fraction (overridden if phi_ref_ic is coupled)
    phi_cell_0 = 0.10
    phi_ref_ic = phi_ref_ic

    # FD tangent
    epsilon  = 1e-6
  []

  # Pressure gate OFF (gp≈1): make press_str huge
  [cell_gel_no_pressure_gate]
    type = CellGelMixtureOpt
    G_cell    = 1.0
    q_cell    = -2.0
    k_exp_max = .5
    c1        = 20.0
    c2        = 2.0
    press_str = 1e30
    k_T1_max  = 0.5
    chi_str   = 0.3
    beta_T1   = 10.0
    m_T1      = 0.0
    G_gel     = 0.01
    k_diss_0  = 0.05
    n         = n
    D_nutrient = 0.1
    n_c1      = 0.05
    n_c2      = 2.0
    n_eps     = 1e-8
    gamma_n0  = 20.0
    phi_cell_0 = 0.10
    phi_ref_ic = phi_ref_ic
    epsilon   = 1e-6
  []

  # Nutrient gate OFF (fa≈1): force n_c1→tiny AND increase n_eps so n_+ never reaches exact 0
  [cell_gel_no_nutrient_gate]
    type = CellGelMixtureOpt
    G_cell    = 1.0
    q_cell    = -2.0
    k_exp_max = .5
    c1        = 20.0
    c2        = 2.0
    press_str = 1.0
    k_T1_max  = 0.5
    chi_str   = 0.3
    beta_T1   = 10.0
    m_T1      = 0.0
    G_gel     = 0.01
    k_diss_0  = 0.05
    n         = n
    D_nutrient = 0.1
    n_c1      = 1e-30
    n_c2      = 2.0
    n_eps     = 1e-6
    gamma_n0  = 20.0
    phi_cell_0 = 0.10
    phi_ref_ic = phi_ref_ic
    epsilon   = 1e-6
  []

  # Consumption OFF (reaction kernel remains, but source becomes 0)
  [cell_gel_no_consumption]
    type = CellGelMixtureOpt
    G_cell    = 1.0
    q_cell    = -2.0
    k_exp_max = .5
    c1        = 20.0
    c2        = 2.0
    press_str = 1.0
    k_T1_max  = 0.5
    chi_str   = 0.3
    beta_T1   = 10.0
    m_T1      = 0.0
    G_gel     = 0.01
    k_diss_0  = 0.05
    n         = n
    D_nutrient = 0.1
    n_c1      = 0.05
    n_c2      = 2.0
    n_eps     = 1e-8
    gamma_n0  = 0.0
    phi_cell_0 = 0.10
    phi_ref_ic = phi_ref_ic
    epsilon   = 1e-6
  []

  # Mechanical-only growth (no nutrient gate, no pressure gate, no consumption, no T1)
  [cell_gel_mech_only]
    type = CellGelMixtureOpt
    G_cell    = 1.0
    q_cell    = -2.0
    k_exp_max = .5
    c1        = 20.0
    c2        = 2.0
    press_str = 1e30
    k_T1_max  = 0.0
    chi_str   = 0.3
    beta_T1   = 10.0
    m_T1      = 0.0
    G_gel     = 0.01
    k_diss_0  = 0.05
    n         = n
    D_nutrient = 0.1
    n_c1      = 1e-30
    n_c2      = 2.0
    n_eps     = 1e-6
    gamma_n0  = 0.0
    phi_cell_0 = 0.10
    phi_ref_ic = phi_ref_ic
    epsilon   = 1e-6
  []
[]

[ICs]
  active = 'n_ic'

  [n_ic]
    type     = ConstantIC
    variable = n
    value    = 0.0
  []

[]

[Kernels]
  # optional but recommended for stability (transient nutrient)
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
    args          = 'n'     # IMPORTANT: do NOT request ux uy uz unless you provide those derivatives
  []
[]

[AuxKernels]
  
  [phi_ref_read]
    type          = SolutionAux
    variable      = phi_ref_ic
    solution      = phi_ref_solution
    from_variable = phi_ref_ic
    direct        = true
    execute_on    = INITIAL
  []

  [c11_aux]
    type            = RankTwoAux
    variable        = c11
    rank_two_tensor = cauchy_stress
    index_i         = 0
    index_j         = 0
    execute_on      = 'INITIAL TIMESTEP_END'
  []

  [vol_aux]
    type      = MaterialRealAux
    variable  = volume_ratio
    property  = volume_ratio
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [ke_aux]
    type      = MaterialRealAux
    variable  = ke
    property  = ke
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [kh_aux]
    type      = MaterialRealAux
    variable  = kh
    property  = kh
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [fa_aux]
    type      = MaterialRealAux
    variable  = fa
    property  = fa
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [press_aux]
    type      = MaterialRealAux
    variable  = press
    property  = pressure
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [gp_aux]
    type      = MaterialRealAux
    variable  = gp
    property  = gp
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [gate_tot_aux]
    type      = MaterialRealAux
    variable  = gate_tot
    property  = gate_tot
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [eta_aux]
    type      = MaterialRealAux
    variable  = eta
    property  = eta
    execute_on = 'INITIAL TIMESTEP_END'
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
[]

[Executioner]
  type        = Transient
  scheme      = bdf2
  solve_type  = NEWTON
  line_search = bt

  # Debug default: run short; increase to 48.0 for full runs
  dtmin    = 1e-6
  dtmax    = 5e-2
  end_time = 24

  nl_max_its = 200
  l_max_its  = 200
  nl_rel_tol  = 1e-5
  nl_abs_tol  = 1e-8

  petsc_options = '-snes_monitor -snes_converged_reason -ksp_converged_reason'
  petsc_options_iname = '-snes_type -snes_linesearch_type -snes_linesearch_max_it -snes_max_it -snes_rtol -snes_atol -ksp_type -pc_type -pc_factor_mat_solver_type'
  petsc_options_value = 'newtonls bt 10 100 1e-4 2e-7 preonly lu mumps'


  automatic_scaling    = true
  compute_scaling_once = false

  [TimeStepper]
    type               = IterationAdaptiveDT
    dt                 = 1e-4
    optimal_iterations = 20
    growth_factor      = 1.15
    cutback_factor     = 0.7
  []
[]

[Preconditioning]
  [pc]
    type = SMP
    full = true
  []
[]

[Outputs]
  [exodus]
    type               = Exodus
    file_base          = YP3D3_cube_debug
    time_step_interval = 1
    execute_on         = 'INITIAL TIMESTEP_END'
    show               = 'ux uy uz n phi_ref_ic phi_ref_aux phi_cell press ke kh fa gp gate_tot eta volume_ratio c11'
  []

  [csv_out]
    type               = CSV
    file_base          = YP3D3_cube_debug_csv
    time_step_interval = 1
    execute_on         = 'TIMESTEP_END'
  []
[]
