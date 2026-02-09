# YP3D4_debugmesh2_pins.i
# - Uses an unstructured circular 2D mesh and extrudes to 3D
# - Removes rigid body modes with a minimal "gauge" pin set (6 scalar constraints)
# - Imports phi_ref_ic as an AuxVariable using SolutionUserObject + SolutionAux

[GlobalParams]
  displacements    = 'ux uy uz'
  large_kinematics = true
[]

[Mesh]
  parallel_type = distributed

  # Read the smooth debugmesh2 puck directly from the filter Exodus.
  # This ensures the mechanics mesh == filter mesh, so SolutionAux direct=true is valid.
  [base]
    type = FileMeshGenerator
    file = phi_ref_filter3d_simple_out.e
  []

  # Minimal rigid-body constraints (same as stable; just applied to the file mesh):
  [pin_center]
    type             = ExtraNodesetGenerator
    input            = base
    coord            = '0 0 0.5'
    use_closest_node = true
    new_boundary     = 'pin_center'
  []

  [pin_x]
    type             = ExtraNodesetGenerator
    input            = pin_center
    coord            = '1 0 0.5'
    use_closest_node = true
    new_boundary     = 'pin_x'
  []

  [pin_y]
    type             = ExtraNodesetGenerator
    input            = pin_x
    coord            = '0 1 0.5'
    use_closest_node = true
    new_boundary     = 'pin_y'
  []

  final_generator = pin_y
[]


[UserObjects]
  [phi_ref_solution]
    type                = SolutionUserObject
    mesh                = phi_ref_filter3d_simple_out.e
    system_variables    = 'phi_ref_ic'
    # Critical for your use case:
    timestep            = LATEST
    # Because your filter mesh is second-order (you converted before extrusion):
    nodal_variable_order = SECOND
    # Helps with fuzzy point location / near-boundary queries:
    # Start moderate; increase if you still see "Failed to access data at point ..."
    epsilon             = 1e-10
    execute_on          = 'INITIAL'
  []
[]


[Variables]
  [n]
    order = FIRST
    family = LAGRANGE
  []
[]

[Physics]
  [SolidMechanics]
    [QuasiStatic]
      [all]
        add_variables      = true
        displacements      = 'ux uy uz'
        strain             = FINITE
        new_system         = true
        formulation        = TOTAL
      []
    []
  []
[]

[AuxVariables]
  [phi_ref_ic]
    family = LAGRANGE
    order  = FIRST
  []

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
[]

[AuxKernels]

  [phi_ref_read]
    type          = SolutionAux
    variable      = phi_ref_ic
    solution      = phi_ref_solution
    from_variable = phi_ref_ic
    # Make this robust in distributed runs:
    direct        = true
    execute_on    = INITIAL
  []


  [c11_aux]
    type            = RankTwoAux
    variable        = c11
    rank_two_tensor = cauchy_stress
    index_i         = 0
    index_j         = 0
  []

  [vol_aux]
    type      = MaterialRealAux
    variable  = volume_ratio
    property  = volume_ratio
  []
  [ke_aux]
    type      = MaterialRealAux
    variable  = ke
    property  = ke
  []
  [kh_aux]
    type      = MaterialRealAux
    variable  = kh
    property  = kh
  []
  [fa_aux]
    type      = MaterialRealAux
    variable  = fa
    property  = fa
  []
  [press_aux]
    type      = MaterialRealAux
    variable  = press
    property  = pressure
  []
  [gp_aux]
    type      = MaterialRealAux
    variable  = gp
    property  = gp
  []
  [gate_tot_aux]
    type      = MaterialRealAux
    variable  = gate_tot
    property  = gate_tot
  []
  [eta_aux]
    type      = MaterialRealAux
    variable  = eta
    property  = eta
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

[ICs]
  [n_init]
    type     = ConstantIC
    variable = n
    value    = 0.02
  []
[]

[BCs]
  # Nutrient bath
  [n_dirichlet_all]
    type     = DirichletBC
    variable = n
    boundary = 'outer top bottom'
    value    = 1.0
  []

  # ---- Minimal rigid-body mode removal (6 scalar constraints; does not impose strain) ----
  # Node A (pin_center): remove translations
  [fix_A_x]
    type     = DirichletBC
    variable = ux
    boundary = pin_center
    value    = 0
  []
  [fix_A_y]
    type     = DirichletBC
    variable = uy
    boundary = pin_center
    value    = 0
  []
  [fix_A_z]
    type     = DirichletBC
    variable = uz
    boundary = pin_center
    value    = 0
  []

  # Node B (pin_x): remove two rotations
  [fix_B_y]
    type     = DirichletBC
    variable = uy
    boundary = pin_x
    value    = 0
  []
  [fix_B_z]
    type     = DirichletBC
    variable = uz
    boundary = pin_x
    value    = 0
  []

  # Node C (pin_y): remove final rotation
  [fix_C_z]
    type     = DirichletBC
    variable = uz
    boundary = pin_y
    value    = 0
  []
[]

[Materials]
  [cell_gel]
    type = CellGelMixtureOpt

    G_cell    = 1.0
    q_cell    = -2.0

    # ---- stability: slower/softer early growth ----
    k_exp_max = 0.1
    c1        = 20.0
    c2        = 2.0

    press_str = 1.0

    k_T1_max  = 0.5
    chi_str   = 0.3
    beta_T1   = 10.0
    m_T1      = 2.0

    G_gel      = 0.01

    # ---- stability: add a little damping in matrix network ----
    k_diss_0   = 0.000

    # ---- nutrient transport params ----
    D_nutrient = 0.1
    phi_cell_0 = 0.10

    # Read the speckle field into this AuxVariable at INITIAL
    phi_ref_ic = phi_ref_ic

    # ---- stability: FD tangent step (your material uses FD consistent tangent) ----
    epsilon    = 1e-8

    # ---- nutrient gate + consumption ----
    n_c1     = 0.5
    n_c2     = 10.0
    gamma_n0 = 20.0
    n        = n
  []
[]


[ICs]
  active = 'n_ic'

  [n_ic]
    type     = ConstantIC
    variable = n
    value    = 1e-5
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


[Executioner]
  type        = Transient
  scheme      = bdf2
  solve_type  = PJFNK
  line_search = bt

  # Time step limits
  dtmin    = 1e-6
  dtmax    = 1e-2
  end_time = 48.0
  nl_rel_tol  = 1e-5
  nl_abs_tol  = 1e-8

  # Nonlinear / linear iteration caps
  nl_max_its = 80
  l_max_its  = 200

  # Direct linear solves (robust, good for debugging)
  #petsc_options = '-snes_monitor -snes_converged_reason -ksp_converged_reason'
  #petsc_options_iname = '-snes_type -snes_linesearch_type -snes_linesearch_max_it -snes_max_it -snes_rtol -snes_atol -ksp_type -pc_type -pc_factor_mat_solver_type'
  #petsc_options_value = 'newtonls bt 10 80 1e-4 1e-7 preonly lu mumps'

  petsc_options_iname = '-snes_type -snes_linesearch_type -ksp_type -ksp_rtol -pc_type -pc_hypre_type -snes_rtol -snes_atol'
  petsc_options_value =  'newtonls bt gmres 1e-5 hypre boomeramg 1e-4 1e-7'

  # Automatic scaling is useful, but computing it every step is expensive
  automatic_scaling    = true
  compute_scaling_once = true

  [TimeStepper]
    type               = IterationAdaptiveDT
    dt                 = 1e-4
    optimal_iterations = 15
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
  exodus = true

  [out]
    type               = Nemesis
    file_base          = YP3D4_stable
    time_step_interval = 50
    execute_on         = 'TIMESTEP_END'
    show               = 'ux uy uz n phi_ref_ic phi_cell press ke kh fa gate_tot'
  []

  [csv_out]
    type               = CSV
    file_base          = YP3D4_stable_csv
    time_step_interval = 10
    execute_on         = 'TIMESTEP_END'
  []

  [chk]
    type               = Checkpoint
    time_step_interval = 50
    num_files          = 2
  []
[]
