# YP3D4_debugmesh2_nullspace.i
# - Uses the (unstructured) circular 2D mesh from debugmesh2.i and extrudes to 3D
# - Removes pin Dirichlet constraints
# - Removes rigid body modes via PETSc nullspace (3 translations + 3 rotations)
# - Imports phi_ref_ic as an AuxVariable using SolutionUserObject + SolutionAux (recommended pattern)

[GlobalParams]
  displacements    = 'ux uy uz'
  large_kinematics = true
[]

[Problem]
  # 3 translations + 3 rotations for a free 3D solid
  null_space_dimension      = 6
  near_null_space_dimension = 6
[]

[Mesh]
  # Most of the unstructured/circle generators are replicated-only.
  # If you want a distributed production run, generate this mesh once (mesh-only) and read it with FileMeshGenerator.
  parallel_type = replicated

  # ----- 2D circle curve -----
  R2  = 4.0
  tol = 0.08

  [circle_curve]
    type = ParsedCurveGenerator
    x_func = '2*cos(t)'
    y_func = '2*sin(t)'
    t_range = '0 6.283185307179586'
    num_points = 150
    closed = true
  []

  # ----- Unstructured (triangular) fill inside the curve -----
  [disk2d]
    type = XYDelaunayGenerator
    input = circle_curve
    boundary = 'circle'
    # Controls element count: ~ area / desired_area. For R=2, area ~ 12.57.
    desired_area = 0.002
    output_subdomain_name = 'puck'
    show_info = false
  []

  # Tag the external circular sides as a named sideset.

  [outer_sideset]
    type = ParsedGenerateSideset
    input = disk2d
    combinatorial_geometry = 'abs(x*x + y*y - R2) < tol'
    include_only_external_sides = true
    new_sideset_name = 'outer'
  []

  # Option B: upgrade element order so the boundary mapping is curved (less faceted) without increasing element count.
  [p2_2d]
    type = ElementOrderConversionGenerator
    input = outer_sideset
    conversion_type = SECOND_ORDER
  []

  # ----- Extrude to 3D puck -----
  [cyl]
    type = AdvancedExtruderGenerator
    input = p2_2d
    direction = '0 0 1'
    heights = '1.0'
    num_layers = '30'
    bottom_boundary = 'bottom'
    top_boundary = 'top'
  []

  final_generator = cyl
[]

# -------------------- Nullspace construction (rigid body modes) --------------------
# MOOSE allocates PETSc vectors named NullSpace_0.. and NearNullSpace_0.. when you set
# null_space_dimension / near_null_space_dimension. Then RigidBodyModes3D fills them.

[UserObjects]
  [rbm_null]
    type = RigidBodyModes3D
    block = 'puck'
    boundary = 'outer top bottom'

    disp_x = 'ux'
    disp_y = 'uy'
    disp_z = 'uz'

    subspace_name = 'NullSpace'
    subspace_indices = '0 1 2 3 4 5'
    modes = 'trans_x trans_y trans_z rot_x rot_y rot_z'

    execute_on = 'INITIAL'
  []

  [rbm_near]
    type = RigidBodyModes3D
    block = 'puck'
    boundary = 'outer top bottom'

    disp_x = 'ux'
    disp_y = 'uy'
    disp_z = 'uz'

    subspace_name = 'NearNullSpace'
    subspace_indices = '0 1 2 3 4 5'
    modes = 'trans_x trans_y trans_z rot_x rot_y rot_z'

    execute_on = 'INITIAL'
  []

  # Read phi_ref_ic from the filtering run output (must be on the same mesh, or use interpolation).
  [phi_ref_solution]
    type             = SolutionUserObject
    mesh             = phi_ref_filter3d_debugmesh2_out.e
    system_variables = 'phi_ref_ic'
    execute_on       = 'INITIAL'
  []
[]

[Functions]
  [gate_n_func]
    type = PiecewiseLinear
    x = '0 0.0032 0.0064 0.0096 0.0128 0.016'
    y = '0 0.2    0.4    0.6    0.8    1.0'
  []

  [pressure_gate]
    type = ParsedFunction
    expression = '0.5*(1.0 - tanh((p - 0.2)/0.05))'
    vars = 'p'
    vals = 'press'
  []

  [C1_gate]
    type = ParsedFunction
    expression = 'x*x*(3-2*x)'
    vars = 'x'
    vals = 'n'
  []

  [c1]
    type = ParsedFunction
    expression = 'max(0.0, min(1.0, (n / ncrit)))'
    vars = 'n ncrit'
    vals = 'n 0.01'
  []
[]

[Variables]
  [n]
    order = FIRST
    family = LAGRANGE
  []
[]

[AuxVariables]
  # Imported heterogeneity field (constant in time for this run)
  [phi_ref_ic]
    order  = FIRST
    family = LAGRANGE
  []

  [phi_ref_aux]
    order = FIRST
    family = LAGRANGE
  []

  [phi_cell]
    order = FIRST
    family = LAGRANGE
  []

  [press]
    order = FIRST
    family = LAGRANGE
  []

  [ke]
    order = FIRST
    family = LAGRANGE
  []

  [kh]
    order = FIRST
    family = LAGRANGE
  []

  [fa]
    order = FIRST
    family = LAGRANGE
  []

  [gate_tot]
    order = FIRST
    family = LAGRANGE
  []
[]

[AuxKernels]
  # Pull phi_ref_ic from the filter-output file into an AuxVariable (recommended pattern)
  [phi_ref_read]
    type          = SolutionAux
    variable      = phi_ref_ic
    solution      = phi_ref_solution
    from_variable = phi_ref_ic
    # If the mesh is identical, set direct = true for exact DOF mapping.
    direct        = true
    execute_on    = 'INITIAL'
    block         = 'puck'
  []

  [press_aux_k]
    type       = MaterialRealAux
    variable   = press
    property   = pbar
    execute_on = 'INITIAL TIMESTEP_END'
  []

  [ke_aux_k]
    type       = MaterialRealAux
    variable   = ke
    property   = kdiv
    execute_on = 'INITIAL TIMESTEP_END'
  []

  [kh_aux_k]
    type       = MaterialRealAux
    variable   = kh
    property   = kh
    execute_on = 'INITIAL TIMESTEP_END'
  []

  [fa_aux_k]
    type       = MaterialRealAux
    variable   = fa
    property   = fa
    execute_on = 'INITIAL TIMESTEP_END'
  []

  [gate_aux_k]
    type       = MaterialRealAux
    variable   = gate_tot
    property   = gate_tot
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

[ICs]
  [n_init]
    type     = ConstantIC
    variable = n
    value    = 0.02
  []
[]

[Kernels]
  [n_td]
    type     = TimeDerivative
    variable = n
  []

  [n_diff]
    type     = Diffusion
    variable = n
  []
[]

[BCs]
  # Nutrient bath
  [n_outer]
    type     = DirichletBC
    variable = n
    boundary = 'outer'
    value    = 0.02
  []
[]

[Materials]
  [cellgel]
    type = CellGelMixtureOpt
    block = puck

    # Coupled variables
    n = n
    phi_ref_ic = phi_ref_ic
    use_phi_ref = true

    # Switches
    use_nutrient         = true
    nutrient_gate_func   = gate_n_func
    use_pressure_gate    = true
    pressure_gate_func   = pressure_gate
    use_nutrient_gate2   = true
    nutrient_gate2_func  = c1

    # Base params
    phi_cell0      = 0.1
    Ec_ref         = 0.2
    Eg_ref         = 2.0
    kc0_ref        = 0.75
    kg0_ref        = 0.0
    mu_gel_ref     = 0.0005
    eta_gel_ref    = 0.0001
    eta_cell_ref   = 0.002

    # Growth / kinetics
    use_growth     = true
    k_exp_max      = 0.0433
    k_div_max      = 0.05
    kh_slope       = 2.0
    kh_offset      = 0.1

    # Numerical
    F_damp = 0.95

    # Outputs
    output_pbar            = true
    output_Jc              = true
    output_Jg              = true
    output_kdiv            = true
    output_kh              = true
    output_fa              = true
    output_phi_cell_ref    = true
    output_phi_cell        = true
    output_gate_tot        = true
    output_pk2_stress      = true
  []
[]

[Physics]
  [SolidMechanics]
    [QuasiStatic]
      new_system        = true
      strain            = FINITE
      use_displaced_mesh = true
      stress            = 'cellgel'
      add_variables     = true
      volumetric_locking_correction = true
    []
  []
[]

[Executioner]
  type        = Transient
  scheme      = bdf2
  solve_type  = NEWTON
  line_search = bt

  dtmin    = 1e-6
  dtmax    = 5e-3
  end_time = 48.0
  nl_rel_tol  = 1e-5
  nl_abs_tol  = 1e-8

  nl_max_its = 80
  l_max_its  = 200

  # Keep your LU+MUMPS defaults for now.
  # If you see singular-factorization complaints, switch to an iterative+AMG stack and keep NearNullSpace enabled.
  petsc_options = '-snes_monitor -snes_converged_reason -ksp_converged_reason'
  petsc_options_iname = '-snes_type -snes_linesearch_type -snes_linesearch_max_it -snes_max_it -snes_rtol -snes_atol -ksp_type -pc_type -pc_factor_mat_solver_type'
  petsc_options_value = 'newtonls bt 10 80 1e-4 1e-7 preonly lu mumps'

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
    type               = Exodus
    file_base          = YP3D4_debugmesh2_nullspace
    time_step_interval = 50
    execute_on         = 'TIMESTEP_END'
    show               = 'ux uy uz n phi_ref_ic phi_cell press ke kh fa gate_tot'
  []

  [csv_out]
    type               = CSV
    file_base          = YP3D4_debugmesh2_nullspace_csv
    time_step_interval = 10
    execute_on         = 'TIMESTEP_END'
  []

  [chk]
    type               = Checkpoint
    time_step_interval = 50
    num_files          = 2
  []
[]
