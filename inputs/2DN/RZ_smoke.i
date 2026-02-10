# -----------------------------------------------------------------------------
# RZ_smoke.i  (fast integration test for AD + SolidMechanics new system)
# Exercises:
#   - SolidMechanics QuasiStatic new_system + finite strain + AD
#   - CellGelMixtureOpt + ADNutrientTLTransport
#   - ADJnTimeDerivative, ADMatTensorDiffusion, ADMatReactionSigned
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

    # Small mesh for speed (smoke test)
    nx   = 30
    ny   = 10
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
  # reference loading field (coupled into CellGelMixtureOpt)
  [phi_ref_ic]
    family = LAGRANGE
    order  = FIRST
  []

  # Minimal diagnostics
  [volume_ratio]
    family = MONOMIAL
    order  = CONSTANT
  []
  [phi_cell]
    family = MONOMIAL
    order  = CONSTANT
  []
  [D_phys_out]
    family = MONOMIAL
    order  = CONSTANT
  []
[]

[ICs]
  active = 'phi_ref_ic_rand n0'

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

    # nutrient coupling
    D_nutrient = 0.1
    phi_cell_0 = 0.10
    phi_ref_ic = phi_ref_ic
    epsilon    = 1e-6

    n_c1     = 0.05
    n_c2     = 12.0
    gamma_n0 = 20.0
    n        = n
  []

  [nutrient_tl]
    type  = ADNutrientTLTransport
    block = 0

    n = n

    # TEMP: keep these only to satisfy current required params
    disp_r = ux
    disp_z = uy

    axisymmetric = true
    radial_coord = 0

    # expects the property produced by the mechanics material
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
  [vol_aux]
    type     = MaterialRealAux
    variable = volume_ratio
    property = volume_ratio
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [phi_cell_aux]
    type     = MaterialRealAux
    variable = phi_cell
    property = phi_cell
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [D_phys_aux]
    type     = MaterialRealAux
    variable = D_phys_out
    property = D_phys
    execute_on = 'INITIAL TIMESTEP_END'
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
  scheme      = implicit-euler
  solve_type  = NEWTON
  line_search = bt

  # Very short run
  dt       = 1e-2
  end_time = 5e-2

  nl_rel_tol = 1e-6
  nl_abs_tol = 1e-9
  nl_max_its = 25

  # Robust small-problem solve
  petsc_options_iname = '-snes_type -snes_linesearch_type -ksp_type -pc_type -snes_rtol -snes_atol'
  petsc_options_value =  'newtonls bt preonly lu 1e-8 1e-10'
[]

[Outputs]
  exodus = true
  csv    = true
  execute_on = 'FINAL'
[]
