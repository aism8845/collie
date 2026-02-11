# -----------------------------------------------------------------------------
# RZexp_ELM_nutrient_compare.i
#
# Axisymmetric (RZ) 2D yeast-puck growth model using your CellGelMixtureOpt.
# Goal: reproduce Figure_ELM Nutrient-style volume change trends:
#   V/V0 vs (i) yeast loading (phi_cell_0), (ii) boundary glucose (n_bc),
#   (iii) puck thickness (ymax), for a fixed puck radius (xmax).
#
# Coordinate convention (MOOSE axisymmetric): displacements ordered (u_r, u_z)
# so here: ux ≡ u_r, uy ≡ u_z.
#
# Geometry scaling (matches your 2Dexp files):
#   Interpret length unit L0 ≈ 2 mm.
#   Then: radius xmax=1.0 -> 2 mm (diameter 4 mm); thickness ymax=0.5 -> 1 mm.
#
# Sweep knobs (recommended overrides):
#   - yeast loading:     Materials/cell_gel/phi_cell_0 = 0.02,0.05,0.10,0.20,...
#   - bath glucose:      BCs/n_top/value = 0.01,0.1,1,10  (if n normalized by 20 mg/mL)
#   - thickness:         Mesh/gm/ymax = 0.25,0.5,1.0      (-> 0.5,1,2 mm if L0=2mm)
#   - gel stiffness:     Materials/cell_gel/G_gel = ...
# -----------------------------------------------------------------------------
[GlobalParams]
  displacements = 'ux uy'
[]


[Mesh]
  # Prefer Mesh/coord_type (Problem/coord_type is deprecated in newer MOOSE)
  coord_type    = RZ
  rz_coord_axis = Y
  parallel_type = distributed  # matches 2Dexp1; OK for SolutionUserObject   # axial = y, radial = x  (so r = x, z = y)

  [gm]
    type = GeneratedMeshGenerator
    dim  = 2
    nx   = 80
    ny   = 40
    xmin = 0.0
    xmax = 2.0     # radius (matches 2Dexp1 exodus)
    ymin = 0.0
    ymax = 0.5     # thickness (matches 2Dexp1 exodus)
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

  # Diagnostics pulled from material properties
  [J]            # "volume_ratio" in the material (det(F))
    family = MONOMIAL
    order  = CONSTANT
  []
  [phi_cell_out] # current cell volume fraction from material
    family = MONOMIAL
    order  = CONSTANT
  []
  [fa_out]       # nutrient gate f_a(n) from material
    family = MONOMIAL
    order  = CONSTANT
  []
  [press_out]    # pressure-like measure from material (for debugging)
    family = MONOMIAL
    order  = CONSTANT
  []
[]


[UserObjects]
  # Reads phi_ref_ic from the Exodus output used in 2Dexp1 (phi_ref_filter_out2.e)
  # Put phi_ref_filter_out2.e in your run directory (or use an absolute path).
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
      # defaults (propagated to subblocks)
      displacements = 'ux uy'
      add_variables = true
      strain        = FINITE

      # If you are using the *new Lagrangian kernel system*:
      new_system   = true
      formulation  = TOTAL   # IMPORTANT for axisym (see below)

      [all]
        # empty is OK: presence of a subblock triggers kernel/material construction
      []
    []
  []
[]

[Kernels]
  # Nutrient diffusion + consumption (transient)
  [n_td]
    type     = TimeDerivative
    variable = n
  []
  [n_diff]
    type        = MatAnisoDiffusion
    variable    = n
    diffusivity = D_eff          # provided by CellGelMixtureOpt
  []
  [n_rxn]
    type          = MatReaction
    variable      = n
    reaction_rate = n_source_ref  # provided by CellGelMixtureOpt
    args          = n
  []
[]

[BCs]
  # Activate 'n_right' too if you want side-feeding in addition to top-feeding.
  active = 'ux_axis uy_bottom n_top n_right'

  # Axisymmetry: u_r(r=0) = 0
  [ux_axis]
    type     = DirichletBC
    variable = ux
    boundary = left
    value    = 0.0
  []

  # Anchor the puck axially at the bottom (prevents rigid translation in z)
  [uy_bottom]
    type     = DirichletBC
    variable = uy
    boundary = bottom
    value    = 0.0
  []

  # Nutrient bath at the top face: n = n_bc
  [n_top]
    type     = DirichletBC
    variable = n
    boundary = top
    value    = 1.0    # <-- default normalized bath concentration; override as needed
  []

  # Optional: side-feeding (outer radius boundary)
  [n_right]
    type     = DirichletBC
    variable = n
    boundary = right
    value    = 1.0
  []
[]

[Materials]
  [cell_gel]
    type = CellGelMixtureOpt
    block = 0

    # --- mechanics / growth ---
    G_cell    = 1.0
    q_cell    = -2.0

    k_exp_max = 0.2
    c1        = 10.0
    c2        = 2.0

    press_str = 1.0

    k_T1_max  = 0.0
    chi_str   = 1.0
    beta_T1   = 0.0
    m_T1      = 2.0

    # Gel stiffness (relative); tune to match compressive modulus groups
    G_gel     = 0.01

    # Dissolution/cleavage (off here)
    k_diss_0  = 0.0

    # --- nutrient transport/consumption (coupled) ---
    D_nutrient = 0.1
    phi_cell_0 = 0.10

    # Initial cell fraction field (speckle, etc) – here uniform constant IC
    phi_ref_ic = phi_ref_ic

    # Numerical regularization
    epsilon = 1e-6

    # Nutrient gate: fa(n)= n^m/(n^m + n_c1^m), with m = n_c2
    n_c1     = 0.05
    n_c2     = 12.0
    gamma_n0 = 10.0
    n        = n
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
  # In RZ, PP integrals account for the axisymmetric coordinate system.
  # avg_J = (1/V0)∫ J dV0 = V/V0 exactly.
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
  end_time = 48.0

  nl_rel_tol = 1e-5
  nl_abs_tol = 1e-8
  nl_max_its = 100

  petsc_options_iname = '-snes_type -snes_linesearch_type -ksp_type -ksp_rtol -pc_type -pc_hypre_type -snes_rtol -snes_atol'
  petsc_options_value =  'newtonls bt gmres 1e-5 hypre boomeramg 1e-4 1e-7'

   automatic_scaling    = true

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
