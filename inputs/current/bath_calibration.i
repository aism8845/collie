# -----------------------------------------------------------------------------
# bath_calibration.i  (finite-bath depletion calibration, 24 h run)
# -----------------------------------------------------------------------------
# Scaling / normalizers (experiment-normalized nondimensionalization):
#   t0 = 1 hour        => 1 sim second = 1 hour
#   L0 = 1 mm          => mesh coordinates are in mm
#   C0 = 20 mg/mL      => ~111 mM glucose, n = C/C0, bath = C_bath/C0
#   G0 = 0.7 MPa       => p* = p_phys / G0
# Geometry mapping:
#   R_phys = 2 mm  -> xmax = 2.0
#   H_phys = 1 mm  -> this input uses a half-height model with ymax = 0.5
#                    (symmetry about the puck midplane).
# Bath volume mapping:
#   V_bath_phys = 10 mL = 10,000 mm^3 -> V_bath* = 10000 (with L0 = 1 mm)
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

  # ---- Mechanics admissibility diagnostics ----
  [J_mech_aux]
    family = MONOMIAL
    order  = CONSTANT
  []
  [det_bE_cell_aux]
    family = MONOMIAL
    order  = CONSTANT
  []
  [min_eig_bE_cell_aux]
    family = MONOMIAL
    order  = CONSTANT
  []
  [det_bE_pmat_aux]
    family = MONOMIAL
    order  = CONSTANT
  []
  [min_eig_bE_pmat_aux]
    family = MONOMIAL
    order  = CONSTANT
  []
  [be_cell_nonspd_aux]
    family = MONOMIAL
    order  = CONSTANT
  []
  [be_pmat_nonspd_aux]
    family = MONOMIAL
    order  = CONSTANT
  []
  [be_cell_projected_aux]
    family = MONOMIAL
    order  = CONSTANT
  []
  [be_pmat_projected_aux]
    family = MONOMIAL
    order  = CONSTANT
  []
  [be_cell_projection_amount_aux]
    family = MONOMIAL
    order  = CONSTANT
  []
  [be_pmat_projection_amount_aux]
    family = MONOMIAL
    order  = CONSTANT
  []
  [F_inv_guard_nutr_aux]
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
    # Property provided by ADNutrientTLTransport (audited in src/materials/ADNutrientTLTransport.C).
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
    V_bath = 5000
    k_refill = k_refill_zero
    n_feed = 1.0
  []
[]


[BCs]
  active = 'ux_axis uy_bottom n_sym_bottom n_bath'

  # Axisymmetry for mechanics
  [ux_axis]
    type     = DirichletBC
    variable = ux
    boundary = left
    value    = 0.0
  []

  # Midplane symmetry for mechanics: uy=0, ux free
  [uy_bottom]
    type     = DirichletBC
    variable = uy
    boundary = bottom
    value    = 0.0
  []

  # Symmetry (no diffusive flux) for nutrient on the midplane
  [n_sym_bottom]
    type     = NeumannBC
    variable = n
    boundary = bottom
    value    = 0.0
  []

  # Bath Dirichlet only on exposed surfaces (top + right)
  [n_bath]
    type     = ScalarCoupledDirichletBC
    variable = n
    boundary = 'top right'
    bath     = bath
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
    # press_str = P_str/G0; with G0=0.7 MPa, press_str=0.64 -> P_str ~= 0.45 MPa.
    press_str = 0.64
    # Tier-1 regularization: C2-smooth pressure gate blend around p_cell = 0.
    press_gate_smooth = 0.005

    G_gel     = 0.01 
    k_diss_0  = 0.0

    # --- nutrient transport/consumption (coupled) ---
    D_nutrient = 1.0
    phi_cell_0 = 0.60

    phi_ref_ic = phi_ref_ic

    epsilon = 1e-6

    n_c1     = 0.005
    # Nutrient-gate regularization (fix #4): soften Hill steepness at low n.
    n_c2     = 2.0
    # Keep synced with nutrient_tl/gamma_n0 for consistent reporting, although n_rxn uses n_source_ref_nutr.
    gamma_n0 = 500.0
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
    # D0 carries the nutrient diffusivity scale used by D_eff_nutr in the nutrient PDE.
    D0      = 2.7
    D_floor = 1e-12
    gamma_n0 = 500.0
    phi_max  = 0.65
    n_c1     = 0.005
    n_c2     = 2
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

  [J_mech_aux_k]
    type     = MaterialRealAux
    variable = J_mech_aux
    property = J_mech
  []
  [det_bE_cell_aux_k]
    type     = MaterialRealAux
    variable = det_bE_cell_aux
    property = det_bE_cell
  []
  [min_eig_bE_cell_aux_k]
    type     = MaterialRealAux
    variable = min_eig_bE_cell_aux
    property = min_eig_bE_cell
  []
  [det_bE_pmat_aux_k]
    type     = MaterialRealAux
    variable = det_bE_pmat_aux
    property = det_bE_pmat
  []
  [min_eig_bE_pmat_aux_k]
    type     = MaterialRealAux
    variable = min_eig_bE_pmat_aux
    property = min_eig_bE_pmat
  []
  [be_cell_nonspd_aux_k]
    type     = MaterialRealAux
    variable = be_cell_nonspd_aux
    property = be_cell_nonspd_flag
  []
  [be_pmat_nonspd_aux_k]
    type     = MaterialRealAux
    variable = be_pmat_nonspd_aux
    property = be_pmat_nonspd_flag
  []
  [be_cell_projected_aux_k]
    type     = MaterialRealAux
    variable = be_cell_projected_aux
    property = be_cell_projected_flag
  []
  [be_pmat_projected_aux_k]
    type     = MaterialRealAux
    variable = be_pmat_projected_aux
    property = be_pmat_projected_flag
  []
  [be_cell_projection_amount_aux_k]
    type     = MaterialRealAux
    variable = be_cell_projection_amount_aux
    property = be_cell_projection_amount
  []
  [be_pmat_projection_amount_aux_k]
    type     = MaterialRealAux
    variable = be_pmat_projection_amount_aux
    property = be_pmat_projection_amount
  []
  [F_inv_guard_nutr_aux_k]
    type     = MaterialRealAux
    variable = F_inv_guard_nutr_aux
    property = F_inv_guard_flag_nutr
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
    # SideTensorDiffusiveFluxIntegral retrieves non-AD tensor properties.
    # D_eff_nutr is AD-only here, so this must use the non-AD tensor D_eff.
    diffusivity = D_eff
    boundary = 'top right'
    execute_on = 'timestep_end'
    outputs = 'solver_watch'
  []
  [flux_top]
    type = SideTensorDiffusiveFluxIntegral
    variable = n
    diffusivity = D_eff
    boundary = top
    execute_on = 'timestep_end'
    outputs = 'solver_watch'
  []
  [flux_right]
    type = SideTensorDiffusiveFluxIntegral
    variable = n
    diffusivity = D_eff
    boundary = right
    execute_on = 'timestep_end'
    outputs = 'solver_watch'
  []
  [flux_full]
    type = ParsedPostprocessor
    expression = 'flux_top + 2.0*flux_right'
    pp_names = 'flux_top flux_right'
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
  [NumFixedPointIterations]
    type = NumFixedPointIterations
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

  [min_detF]
    type       = ElementExtremeValue
    variable   = J_mech_aux
    value_type = min
    outputs = 'solver_watch'
  []
  [min_det_bE_cell]
    type       = ElementExtremeValue
    variable   = det_bE_cell_aux
    value_type = min
    outputs = 'solver_watch'
  []
  [min_min_eig_bE_cell]
    type       = ElementExtremeValue
    variable   = min_eig_bE_cell_aux
    value_type = min
    outputs = 'solver_watch'
  []
  [min_det_bE_pmat]
    type       = ElementExtremeValue
    variable   = det_bE_pmat_aux
    value_type = min
    outputs = 'solver_watch'
  []
  [min_min_eig_bE_pmat]
    type       = ElementExtremeValue
    variable   = min_eig_bE_pmat_aux
    value_type = min
    outputs = 'solver_watch'
  []
  [max_be_cell_nonspd]
    type       = ElementExtremeValue
    variable   = be_cell_nonspd_aux
    value_type = max
    outputs = 'solver_watch'
  []
  [max_be_pmat_nonspd]
    type       = ElementExtremeValue
    variable   = be_pmat_nonspd_aux
    value_type = max
    outputs = 'solver_watch'
  []
  [max_be_cell_projected]
    type       = ElementExtremeValue
    variable   = be_cell_projected_aux
    value_type = max
    outputs = 'solver_watch'
  []
  [max_be_pmat_projected]
    type       = ElementExtremeValue
    variable   = be_pmat_projected_aux
    value_type = max
    outputs = 'solver_watch'
  []
  [max_be_cell_projection_amount]
    type       = ElementExtremeValue
    variable   = be_cell_projection_amount_aux
    value_type = max
    outputs = 'solver_watch'
  []
  [max_be_pmat_projection_amount]
    type       = ElementExtremeValue
    variable   = be_pmat_projection_amount_aux
    value_type = max
    outputs = 'solver_watch'
  []
  [max_F_inv_guard_nutr]
    type       = ElementExtremeValue
    variable   = F_inv_guard_nutr_aux
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
  # Default baseline remains SMP/LU.
  # Use override files to switch Preconditioning/active between:
  #   pc_smp_lu, pc_fsp3_schur_amg, pc_fsp3_schur_gamg
  active = 'pc_smp_lu'

  [pc_smp_lu]
    type = SMP
    full = true
  []

  [pc_fsp3_schur_amg]
    type = FSP
    topsplit = 'mn_bath'

    # PETSc Schur requires two fields at each split level.
    # Use nested Schur splits to keep 3 physics blocks:
    #   top: (mech+nutr) vs bath
    #   nested inside (mech+nutr): mech vs nutr
    [mn_bath]
      splitting = 'mn bath'
      splitting_type = schur
      # Recommended Schur controls for robust factorization/preconditioning.
      petsc_options_iname = '-pc_fieldsplit_schur_fact_type -pc_fieldsplit_schur_precondition'
      petsc_options_value = 'full selfp'
    []

    [mn]
      vars = 'ux uy n'
      splitting = 'mech nutr'
      splitting_type = schur
      petsc_options_iname = '-pc_fieldsplit_schur_fact_type -pc_fieldsplit_schur_precondition'
      petsc_options_value = 'full selfp'
    []

    [mech]
      vars = 'ux uy'
      petsc_options_iname = '-ksp_type -pc_type -pc_factor_mat_solver_type'
      petsc_options_value = 'preonly lu mumps'
    []

    [nutr]
      vars = 'n'
      petsc_options_iname = '-ksp_type -pc_type -pc_hypre_type'
      petsc_options_value = 'preonly hypre boomeramg'
    []

    [bath]
      vars = 'bath'
      petsc_options_iname = '-ksp_type -pc_type'
      petsc_options_value = 'preonly lu'
    []
  []

  [pc_fsp3_schur_gamg]
    type = FSP
    topsplit = 'mn_bath'

    [mn_bath]
      splitting = 'mn bath'
      splitting_type = schur
      petsc_options_iname = '-pc_fieldsplit_schur_fact_type -pc_fieldsplit_schur_precondition'
      petsc_options_value = 'full selfp'
    []

    [mn]
      vars = 'ux uy n'
      splitting = 'mech nutr'
      splitting_type = schur
      petsc_options_iname = '-pc_fieldsplit_schur_fact_type -pc_fieldsplit_schur_precondition'
      petsc_options_value = 'full selfp'
    []

    [mech]
      vars = 'ux uy'
      petsc_options_iname = '-ksp_type -pc_type -pc_factor_mat_solver_type'
      petsc_options_value = 'preonly lu mumps'
    []

    [nutr]
      vars = 'n'
      petsc_options_iname = '-ksp_type -pc_type'
      petsc_options_value = 'preonly gamg'
    []

    [bath]
      vars = 'bath'
      petsc_options_iname = '-ksp_type -pc_type'
      petsc_options_value = 'preonly lu'
    []
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
  # Allow the nonlinear solve to work through the late-time stiff window
  # before forcing aggressive global cutbacks.
  nl_max_its = 160

  petsc_options = '-snes_converged_reason -ksp_converged_reason'
  petsc_options_iname = '-snes_type -snes_linesearch_type -ksp_type -pc_type -pc_factor_mat_solver_type -snes_rtol -snes_atol'
  petsc_options_value =  'newtonls bt preonly lu mumps 1e-4 1e-7'

  automatic_scaling = true

  [TimeStepper]
    type               = IterationAdaptiveDT
    # Start slightly smaller and ramp more gently to reduce early overshoot.
    dt                 = 5e-3
    growth_factor      = 1.2
    cutback_factor     = 0.5
    # Robustness tuning: keep dt at 0.05 longer unless iterations are very high.
    optimal_iterations = 40
    iteration_window   = 6
  []
[]

[Outputs]
  exodus = true
  perf_graph = false

  [mesh_watch]
    type = CSV
    file_base = outputs/bath_calibration_mesh_watch
    execute_on = 'TIMESTEP_END'
    show = 'dt J_min_1 J_min_1_x J_min_1_y metric2_min_1 metric2_min_1_x metric2_min_1_y J_min_2 J_min_2_x J_min_2_y metric2_min_2 metric2_min_2_x metric2_min_2_y'
  []

  [solver_watch]
    type = CSV
    file_base = outputs/bath_calibration_solver_watch
    execute_on = 'TIMESTEP_END'
    show = 'dt nonlinear_its linear_its J_min_1 bath_value n_min n_max flux_top flux_right flux_full NumFixedPointIterations min_detF min_det_bE_cell min_min_eig_bE_cell min_det_bE_pmat min_min_eig_bE_pmat max_be_cell_nonspd max_be_pmat_nonspd max_be_cell_projected max_be_pmat_projected max_be_cell_projection_amount max_be_pmat_projection_amount max_F_inv_guard_nutr'
  []

  [chk]
    type = Checkpoint
    file_base = outputs/bath_calibration_chk
    num_files = 3
    interval = 1
  []
[]
