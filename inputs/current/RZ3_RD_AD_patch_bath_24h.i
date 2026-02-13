# -----------------------------------------------------------------------------
# RZexp_ELM_nutrient_compare.i  (with YP3D-style debug aux fields + D_iso)
# -----------------------------------------------------------------------------
# Bath calibration mapping for this input: 1 sim second = 1 hour.
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
    nx   = 140
    ny   = 40
    xmin = 0.0
    xmax = 2.0
    ymin = 0.0
    ymax = 0.5
    bias_y = .9
    bias_x = .99
  []
[]

[Variables]
  [n]
    family = LAGRANGE
    order  = FIRST
  []
  [bath]
    family = SCALAR
    order  = FIRST
    initial_condition = 1.0
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

  [elem_quality]
    family = MONOMIAL
    order  = CONSTANT
  []

  # --- Diffusivity diagnostics (ParaView) ---
  [D_iso]
    family = MONOMIAL
    order  = CONSTANT
  []

  [D_phys_out]
    family = MONOMIAL
    order  = CONSTANT
  []
  [D_ref_out]
    family = MONOMIAL
    order  = CONSTANT
  []
[]

[UserObjects]
  [phi_ref_solution]
    type             = SolutionUserObject
    mesh             = phi_ref_filter2d_out.e
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
    value    = 0.5
  []
[]

[Functions]
  [k_refill_zero]
    type = ParsedFunction
    expression = '0'
  []
[]

[Physics]
  [SolidMechanics]
    [QuasiStatic]
      new_system = true
      strain = FINITE
      displacements = 'ux uy'
      add_variables = true
      formulation = TOTAL
      use_automatic_differentiation = false
      [all] []
    []
  []
[]

[Kernels]
  [u_n_offdiag_x]
    type = TLStressDivergenceNutrientOffDiag
    variable = ux
    component = 0
    displacements = 'ux uy'
    large_kinematics = true
    n = n
  []
  [u_n_offdiag_y]
    type = TLStressDivergenceNutrientOffDiag
    variable = uy
    component = 1
    displacements = 'ux uy'
    large_kinematics = true
    n = n
  []
  [n_td]
    type = ADJnTimeDerivative
    variable = n
    coef = J_nutr
    coef_dot = Jdot_nutr
    include_jdot = true
  []
  [n_diff]
    type = ADMatTensorDiffusion
    variable = n
    diffusivity = D_eff_nutr
  []
  [n_rxn]
    type = ADMatReactionSigned
    variable = n
    reaction_rate = n_source_ref_nutr
  []
[]

[ScalarKernels]
  [bath_td]
    type = ODETimeDerivative
    variable = bath
  []
  [bath_mb]
    type = BathMassBalance
    variable = bath
    # BathMassBalance contributes R = -(flux_pp / V_bath + ...),
    # so using outward flux makes bath decrease when nutrient enters tissue.
    flux_pp = flux_out
    # 10 mL bath = 10,000 mm^3 (consistent with mm mesh-length convention).
    V_bath = 10000
    k_refill = k_refill_zero
    n_feed = 1.0
  []
[]


[BCs]
  active = 'ux_axis uy_bottom n_bath'

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

  [n_bath]
    type     = ScalarCoupledDirichletBC
    variable = n
    boundary = 'top right'
    bath = bath
  []
[]

[Materials]
  [cell_gel]
    type  = CellGelMixtureOpt
    block = 0

    # Nutrient normalization: n = C / Cb with Cb = 111 mM (20 mg/mL glucose).
    # So n=1 at the Dirichlet boundary corresponds to 111 mM glucose.

    # --- mechanics / growth ---
    G_cell    = 1.0
    q_cell    = -2.0

    k_exp_max = 0.175
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
    D_nutrient = 1.7
    phi_cell_0 = 0.60

    phi_ref_ic = phi_ref_ic

    epsilon = 1e-6

    n_c1     = 0.005
    n_c2     = 3.0
    gamma_n0 = 10.0
    smooth_eps_c = 1e-6
    n        = n
  []
  [nutrient_tl]
    type  = ADNutrientTLTransport
    block = 0

    # TEMP: keep these only to satisfy current required params
    disp_r = ux
    disp_z = uy

    # keep the rest
    n            = n
    phi_cell_ref = phi_cell_ref
    D0      = 0.5
    D_floor = 1e-12
    gamma_n0 = 1.0
    phi_max  = 0.65
    n_c1     = 0.005
    n_c2     = 3
    smooth_eps_c = 1e-6
  []
  
  [nutr_kin]
    type = ADComputeAxisymmetricRZFiniteStrain
    block = 0
    displacements = 'ux uy'
    base_name = nutr_
  []
[]

[AuxKernels]
  [elem_quality_aux]
    type = ElementQualityAux
    variable = elem_quality
    metric = JACOBIAN
    execute_on = 'initial nonlinear timestep_end'
  []

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
  [D_phys_aux]
    type     = MaterialRealAux
    variable = D_phys_out
    property = D_phys
  []
  [D_ref_aux]
    type     = MaterialRealAux
    variable = D_ref_out
    property = D_ref
  []

[]

[VectorPostprocessors]
  [mesh_distortion_watch]
    type = MeshDistortionWatch
    deformation_gradient_property = deformation_gradient
    exclude_radius_factor = 2.0
    execute_on = 'timestep_end'
    outputs = 'none'
  []
[]

[Postprocessors]
  [flux_out]
    type = SideTensorDiffusiveFluxIntegral
    variable = n
    # This postprocessor is non-AD; use the non-AD transport tensor property.
    diffusivity = D_eff
    boundary = 'top right'
    execute_on = 'timestep_end'
    outputs = 'solver_watch'
  []
  [flux_in]
    type = ParsedPostprocessor
    expression = '-flux_out'
    pp_names = 'flux_out'
    execute_on = 'timestep_end'
    outputs = 'solver_watch'
  []
  [bath_value]
    type = ScalarVariable
    variable = bath
    execute_on = 'initial timestep_end'
    outputs = 'solver_watch'
  []
  [bath_mgml]
    type = ParsedPostprocessor
    expression = '20.0 * bath_value'
    pp_names = 'bath_value'
    execute_on = 'timestep_end'
    outputs = 'solver_watch'
  []
  [bath_mM]
    type = ParsedPostprocessor
    expression = '111.0 * bath_value'
    pp_names = 'bath_value'
    execute_on = 'timestep_end'
    outputs = 'solver_watch'
  []

  [dt]
    type = TimestepSize
    execute_on = 'timestep_end'
    outputs = 'mesh_watch solver_watch'
  []

  [nonlinear_its]
    type = NumNonlinearIterations
    execute_on = 'timestep_end'
    outputs = 'solver_watch'
  []
  [linear_its]
    type = NumLinearIterations
    execute_on = 'timestep_end'
    outputs = 'solver_watch'
  []

  [min_elem_quality]
    type = ElementExtremeValue
    variable = elem_quality
    value_type = min
    execute_on = 'initial timestep_end'
    outputs = 'none'
  []

  [J_min_1]
    type = VectorPostprocessorComponent
    vectorpostprocessor = mesh_distortion_watch
    vector_name = J_min_1
    index = 0
    execute_on = 'timestep_end'
    outputs = 'mesh_watch solver_watch'
  []
  [J_min_1_x]
    type = VectorPostprocessorComponent
    vectorpostprocessor = mesh_distortion_watch
    vector_name = J_min_1_x
    index = 0
    execute_on = 'timestep_end'
    outputs = 'mesh_watch'
  []
  [J_min_1_y]
    type = VectorPostprocessorComponent
    vectorpostprocessor = mesh_distortion_watch
    vector_name = J_min_1_y
    index = 0
    execute_on = 'timestep_end'
    outputs = 'mesh_watch'
  []
  [J_min_2]
    type = VectorPostprocessorComponent
    vectorpostprocessor = mesh_distortion_watch
    vector_name = J_min_2
    index = 0
    execute_on = 'timestep_end'
    outputs = 'mesh_watch'
  []
  [J_min_2_x]
    type = VectorPostprocessorComponent
    vectorpostprocessor = mesh_distortion_watch
    vector_name = J_min_2_x
    index = 0
    execute_on = 'timestep_end'
    outputs = 'mesh_watch'
  []
  [J_min_2_y]
    type = VectorPostprocessorComponent
    vectorpostprocessor = mesh_distortion_watch
    vector_name = J_min_2_y
    index = 0
    execute_on = 'timestep_end'
    outputs = 'mesh_watch'
  []

  [metric2_min_1]
    type = VectorPostprocessorComponent
    vectorpostprocessor = mesh_distortion_watch
    vector_name = metric2_min_1
    index = 0
    execute_on = 'timestep_end'
    outputs = 'mesh_watch solver_watch'
  []
  [metric2_min_1_x]
    type = VectorPostprocessorComponent
    vectorpostprocessor = mesh_distortion_watch
    vector_name = metric2_min_1_x
    index = 0
    execute_on = 'timestep_end'
    outputs = 'mesh_watch'
  []
  [metric2_min_1_y]
    type = VectorPostprocessorComponent
    vectorpostprocessor = mesh_distortion_watch
    vector_name = metric2_min_1_y
    index = 0
    execute_on = 'timestep_end'
    outputs = 'mesh_watch'
  []
  [metric2_min_2]
    type = VectorPostprocessorComponent
    vectorpostprocessor = mesh_distortion_watch
    vector_name = metric2_min_2
    index = 0
    execute_on = 'timestep_end'
    outputs = 'mesh_watch'
  []
  [metric2_min_2_x]
    type = VectorPostprocessorComponent
    vectorpostprocessor = mesh_distortion_watch
    vector_name = metric2_min_2_x
    index = 0
    execute_on = 'timestep_end'
    outputs = 'mesh_watch'
  []
  [metric2_min_2_y]
    type = VectorPostprocessorComponent
    vectorpostprocessor = mesh_distortion_watch
    vector_name = metric2_min_2_y
    index = 0
    execute_on = 'timestep_end'
    outputs = 'mesh_watch'
  []

  [dt_limit_from_J]
    type = ParsedPostprocessor
    pp_names = 'J_min_1'
    expression = 'if(J_min_1 < 0.20, 1e-4, 1e99)'
    execute_on = 'timestep_end'
    outputs = 'none'
  []
  [mesh_distortion_warning]
    type = ParsedPostprocessor
    pp_names = 'J_min_1'
    expression = 'if(J_min_1 < 0.20, 1, 0)'
    execute_on = 'timestep_end'
    outputs = 'console solver_watch'
  []

  [avg_J]
    type     = ElementAverageValue
    variable = volume_ratio
    outputs = 'solver_watch'
  []
  [vol_change_pct]
    type       = ParsedPostprocessor
    expression = '(avg_J - 1.0)*100.0'
    pp_names   = 'avg_J'
    outputs = 'solver_watch'
  []

  [avg_phi_cell]
    type     = ElementAverageValue
    variable = phi_cell
    outputs = 'solver_watch'
  []
  [avg_fa]
    type     = ElementAverageValue
    variable = fa
    outputs = 'solver_watch'
  []
  [avg_ke]
    type     = ElementAverageValue
    variable = ke
    outputs = 'solver_watch'
  []
  [avg_kh]
    type     = ElementAverageValue
    variable = kh
    outputs = 'solver_watch'
  []
  [avg_gp]
    type     = ElementAverageValue
    variable = gp
    outputs = 'solver_watch'
  []
  [avg_press]
    type     = ElementAverageValue
    variable = press
    outputs = 'solver_watch'
  []
  [avg_eta]
    type     = ElementAverageValue
    variable = eta
    outputs = 'solver_watch'
  []
  [avg_gate_tot]
    type     = ElementAverageValue
    variable = gate_tot
    outputs = 'solver_watch'
  []
  [min_gate_tot]
    type       = ElementExtremeValue
    variable   = gate_tot
    value_type = min
    outputs = 'solver_watch'
  []
  [max_gate_tot]
    type       = ElementExtremeValue
    variable   = gate_tot
    value_type = max
    outputs = 'solver_watch'
  []

  [n_min]
    type       = NodalExtremeValue
    variable   = n
    value_type = min
    outputs = 'solver_watch'
  []
  [n_max]
    type       = NodalExtremeValue
    variable   = n
    value_type = max
    outputs = 'solver_watch'
  []
  [n_span]
    type       = ParsedPostprocessor
    expression = 'n_max - n_min'
    pp_names   = 'n_max n_min'
    outputs = 'solver_watch'
  []

  [max_phi_cell]
    type       = ElementExtremeValue
    variable   = phi_cell
    value_type = max
    outputs = 'solver_watch'
  []

  [avg_Diso]
    type     = ElementAverageValue
    variable = D_iso
    outputs = 'solver_watch'
  []

  [avg_Dphys]
    type     = ElementAverageValue
    variable = D_phys_out
    outputs = 'solver_watch'
  []
  [min_Dphys]
    type       = ElementExtremeValue
    variable   = D_phys_out
    value_type = min
    outputs = 'solver_watch'
  []
  [avg_Dref]
    type     = ElementAverageValue
    variable = D_ref_out
    outputs = 'solver_watch'
  []
  [min_Dref]
    type       = ElementExtremeValue
    variable   = D_ref_out
    value_type = min
    outputs = 'solver_watch'
  []
  [D_ref_phys_ratio]
    type       = ParsedPostprocessor
    expression = 'avg_Dref / avg_Dphys'
    pp_names   = 'avg_Dref avg_Dphys'
    outputs = 'solver_watch'
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
  dtmax    = 0.05
  end_time = 24.0

  nl_rel_tol = 1e-5
  nl_abs_tol = 1e-8
  nl_max_its = 100

  petsc_options_iname = '-snes_type -snes_linesearch_type -ksp_type -ksp_rtol -pc_type -sub_pc_type -snes_rtol -snes_atol'
  petsc_options_value =  'newtonls bt fgmres 2e-5 bjacobi ilu 1e-4 1e-7'

  automatic_scaling = true

  [TimeStepper]
    type               = IterationAdaptiveDT
    dt                 = 1e-2
    growth_factor      = 1.3
    cutback_factor     = 0.5
    optimal_iterations = 25
    iteration_window   = 2
  []
[]

[Outputs]
  exodus = false
  perf_graph = false

  [mesh_watch]
    type = CSV
    file_base = outputs/bath_24h_mesh_watch
    execute_on = 'TIMESTEP_END'
    show = 'dt J_min_1 J_min_1_x J_min_1_y metric2_min_1 metric2_min_1_x metric2_min_1_y J_min_2 J_min_2_x J_min_2_y metric2_min_2 metric2_min_2_x metric2_min_2_y'
  []

  [solver_watch]
    type = CSV
    file_base = outputs/bath_24h_solver_watch
    execute_on = 'TIMESTEP_END'
    show = 'dt nonlinear_its linear_its J_min_1 metric2_min_1 mesh_distortion_warning bath_value bath_mgml bath_mM flux_out flux_in avg_J vol_change_pct avg_phi_cell max_phi_cell n_min n_max n_span avg_Dphys min_Dphys avg_Dref min_Dref D_ref_phys_ratio avg_gate_tot min_gate_tot max_gate_tot avg_gp avg_press avg_fa avg_ke avg_kh avg_eta'
  []
[]
