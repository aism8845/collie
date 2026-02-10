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
    D_nutrient = 0.1
    phi_cell_0 = 0.60

    phi_ref_ic = phi_ref_ic

    epsilon = 1e-6

    n_c1     = 0.05
    n_c2     = 12.0
    gamma_n0 = 20.0
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
    D0      = 0.1
    D_floor = 1e-12
    gamma_n0 = 20.0
    phi_max  = 0.65
    n_c1     = 0.05
    n_c2     = 12.0
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
  []
[]

[Postprocessors]
  [dt]
    type = TimestepSize
    execute_on = 'timestep_end'
    outputs = 'mesh_watch'
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
    outputs = 'mesh_watch'
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
    outputs = 'mesh_watch'
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
    outputs = 'console'
  []

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

  [avg_Dphys]
    type     = ElementAverageValue
    variable = D_phys_out
  []
  [avg_Dref]
    type     = ElementAverageValue
    variable = D_ref_out
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
  petsc_options_value =  'newtonls bt gmres 2e-5 hypre boomeramg 1e-4 1e-7'

  automatic_scaling = true

  [TimeStepper]
    type               = IterationAdaptiveDT
    dt                 = 1e-2
    growth_factor      = 1.3
    cutback_factor     = 0.7
    optimal_iterations = 80
    iteration_window   = 2
  []
[]

[Outputs]
  exodus = true
  perf_graph = true

  [mesh_watch]
    type = CSV
    file_base = outputs/mesh_watch
    execute_on = 'TIMESTEP_END'
    show = 'dt J_min_1 J_min_1_x J_min_1_y metric2_min_1 metric2_min_1_x metric2_min_1_y J_min_2 J_min_2_x J_min_2_y metric2_min_2 metric2_min_2_x metric2_min_2_y'
  []
[]
