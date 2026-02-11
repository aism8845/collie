# 2DExp1.i
# 2D plane-strain strip: CellGelMixtureOpt + SolidMechanics (new system) + quasi-steady nutrient gate
# Time units: HOURS
# Goal: fast 2D experimentation with the same material/solver structure as the 3D debug case.

[GlobalParams]
  displacements    = 'ux uy'
  large_kinematics = true
[]

[Mesh]
  parallel_type = replicated  # safest for SolutionUserObject / SolutionIC
  [rect]
    type = GeneratedMeshGenerator
    dim  = 2
    nx   = 200
    ny   = 100
    xmin = 0.0
    xmax = 2.0
    ymin = 0.0
    ymax = 0.5
  []
[]

[UserObjects]
  # Reads phi_ref_ic from the Exodus output of phi_ref_filter.i
  [phi_ref_solution]
    type = SolutionUserObject
    mesh = phi_ref_filter_out2.e
    system_variables = 'phi_ref_ic'
    timestep = LATEST
  []
[]


[Variables]
  # SolidMechanics will create ux, uy
  [n]
    family = LAGRANGE
    order  = FIRST
  []
[]

[Physics]
  [SolidMechanics]
    [QuasiStatic]
      [all]
        add_variables = true
        displacements = 'ux uy'

        strain      = FINITE
        new_system  = true
        formulation = TOTAL
        use_displaced_mesh = true
      []
    []
  []
[]


[BCs]
  # Infinite sheet surrogate in x
  [ux_left]
    type = DirichletBC
    variable = ux
    boundary = left
    value = 0.0
  []

  # Axial clamp (this is the “Euler analogue” ingredient)
  [uy_bottom]
    type = DirichletBC
    variable = uy
    boundary = bottom
    value = 0.0
  []

  [n_top]
    type = DirichletBC
    variable = n
    boundary = top
    value = 1.0
  []

  # no-flux bottom + sides
  [n_bottom]
    type = NeumannBC
    variable = n
    boundary = bottom
    value = 0.0
  []
  [n_left]
    type = NeumannBC
    variable = n
    boundary = left
    value = 0.0
  []
  [n_right]
    type = DirichletBC
    variable = n
    boundary = right
    value = 1.0
  []



[]


[Materials]

  # Nutrient diffusion and (average) consumption rate.
  #
  # Time unit = 1 h, length unit = 1 → 2 mm.
  # D_phys(glucose in 10% PAAm) ~ (3–5)×10^-10 m^2/s ≈ 1–2 mm^2/h
  #   → D* = D_phys t0 / L^2 ≈ 0.25–0.5; choose 0.5 for slightly faster
  #     mixing and numerical robustness.
  #
  # For ~10^9 cells/mL with q_s ~ 5–15 mmol/(gDW·h) and ~3×10^-11 gDW/cell,
  #   volumetric glucose uptake → |gamma_n_phys| ~ 0.02–0.1 h^-1 depending
  #   on exact cell density; we pick an intermediate value and scale it to

  [cell_gel]
    type = CellGelMixtureOpt

    # ----------------- Cell phase elasticity -----------------
    # G_cell sets the stress scale (we normalize by this).
    G_cell    = 1.0
    q_cell    = -2.0

    # ----------------- Growth / expansion drive ---------------
    # Max volumetric expansion rate from unconstrained doubling time Td = 1 h:
    #   k_exp_max = ln(2)/Td ≈ 0.693 h^-1.
    # The C++ uses an envelope k_exp0(t; c1, c2) so that growth ramps up
    # over the first ~10 h and then saturates.
    k_exp_max = 0.7
    c1        = 10.0
    c2        = 2.0

    # ----------------- Pressure gate g_p(p) -------------------
    # press0 is the characteristic hydrostatic cell stress where growth
    # is half-suppressed. With G_cell=1, press0=1 sets this at O(G_cell).
    press_str = 1.0

    # ----------------- T1 neighbor-swapping gate --------------
    # k_T1_max comparable to division rate → rearrangements on O(2 h).
    k_T1_max  = .5
    chi_str   = 0.3
    beta_T1   = 10.0
    m_T1      = 2.0

    # ----------------- Gel network ----------------------------
    # G_gel << G_cell; with G_gel=0.02 the gel is ~50x softer than cells,
    # consistent with ~10–50 kPa gel vs 0.5–2 MPa effective cell stiffness.
    G_gel     = 0.01
    k_diss_0  = 0.0     # turn >0 to add viscoelastic stress relaxation
    D_nutrient = 0.1  # nutrient diffusivity in gel (length^2/time)

    # ----------------- Initial cell fraction ------------------
    # phi_cell_0 ~ 0.1 ↔ O(10^9 cells/mL) given typical yeast volumes,
    # consistent with ~6 wt% dried yeast in the pre-gel solution.
    phi_cell_0 = 0.1

    # Couple a pre-seeded speckle field into the material (used ONLY at t=0 to set phi_cell_ref)
    phi_ref_ic = phi_ref_ic

    # FD tangent eps (used internally by your material)
    epsilon    = 1e-8

    # ----------------- Nutrient gate f_a(n) -------------------
    # f_a(n) = n^n_c2 / (n^n_c2 + n_c1^n_c2)
    # Monod Ks for glucose in S. cerevisiae is typically Ks <= 0.15 g/L;
    # with a 0.1D bath c0 ~ 2 g/L, Ks/c0 ~ 0.05–0.1. We choose n_c1=0.05
    # so growth drops when n < ~5% of bath concentration.
    n_c1 = 0.5
    n_c2 = 2.0
    
    # Baseline nutrient consumption (1/h in your nondimensional time units)
    gamma_n0 = 20

    # Couple nutrient field into CellGel
    n = n
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

# [Dampers]
#   [n_nonneg]
#     type        = BoundingValueNodalDamper
#     variable    = n
#     min_value   = 0.0
#     min_damping = 1e-10
#   []
# []


[AuxVariables]
  [c11]
    family = MONOMIAL
    order  = CONSTANT
  []

  [volume_ratio]   # J = det(F)
    family = MONOMIAL
    order  = CONSTANT
  []

  [ke]             # local volumetric growth rate
    family = MONOMIAL
    order  = CONSTANT
  []

  [kh]             # volumetric strain rate / divergence-like term
    family = MONOMIAL
    order  = CONSTANT
  []

  [fa]             # nutrient gate 0–1
    family = MONOMIAL
    order  = CONSTANT
  []

  [press]          # hydrostatic mixture pressure
    family = MONOMIAL
    order  = CONSTANT
  []

  [gp]             # pressure gate 0–1
    family = MONOMIAL
    order  = CONSTANT
  []

  [gate_tot]       # combined gate f_a * g_p
    family = MONOMIAL
    order  = CONSTANT
  []

  [eta]            # cell number / growth phase (from CellGel)
    family = MONOMIAL
    order  = CONSTANT
  []

  [phi_ref_ic]        # nodal field used only to seed phi_cell_ref (speckle)
    family = LAGRANGE
    order  = FIRST
  []

  [phi_ref_aux]
    family = MONOMIAL
    order  = CONSTANT
  []

  [phi_ref_ic_cell]
    family = MONOMIAL
    order  = CONSTANT
  []

  [phi_cell_aux]
    family = MONOMIAL
    order  = CONSTANT
  []

[]

[AuxKernels]
  # Cauchy σ_11 for quick stress checks
  [c11_aux]
    type            = RankTwoAux
    variable        = c11
    rank_two_tensor = cauchy_stress
    index_i         = 0
    index_j         = 0
  []

  # Direct scalar outputs from CellGel
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
    type     = MaterialRealAux
    variable = phi_ref_aux
    property = phi_cell_ref
    execute_on = 'INITIAL TIMESTEP_END'
  []


  [phi_ref_ic_cell_k]
    type              = ParsedAux
    variable          = phi_ref_ic_cell
    coupled_variables = 'phi_ref_ic'
    expression        = 'phi_ref_ic'
    execute_on        = 'INITIAL TIMESTEP_END'
  []

  [phi_cell_aux_k]
    type       = MaterialRealAux
    variable   = phi_cell_aux
    property   = phi_cell
    execute_on = 'INITIAL TIMESTEP_END'
  []


[]

[ICs]
  active = 'phi_ref_from_file n_ic'

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
    max      = 0.2
    seed     = 123
  []

  [n_ic]
    type     = ConstantIC
    variable = n
    value    = 0.0
  []
[]


[Executioner]
  type        = Transient
  scheme      = bdf2
  solve_type  = NEWTON
  line_search = bt

  # Time is in hours (consistent with your growth/consumption parameters)
  dt       = 5e-3
  dtmin    = 1e-6
  dtmax    = 5e-2
  end_time = 24.0

  # Important for mixed-mechanics + transport scaling
  automatic_scaling = true

  
  nl_rel_tol  = 1e-5
  nl_abs_tol  = 2e-7
  l_max_its   = 50

  [TimeStepper]
    type               = IterationAdaptiveDT
    optimal_iterations = 20
    dt = .001
    growth_factor      = 1.15
    cutback_factor     = 0.7
  []

  # Robust direct solve (matches your 3D debug settings)
  petsc_options = '-snes_monitor -snes_converged_reason -ksp_converged_reason'
  petsc_options_iname  = '-snes_type -snes_linesearch_type -snes_linesearch_max_it -snes_max_it -snes_rtol -snes_atol -ksp_type -pc_type -pc_factor_mat_solver_type'
  petsc_options_value  = 'newtonls bt 10 80 1e-4 2e-7 preonly lu mumps'
[]

[Preconditioning]
  [pc]
    type = SMP
    full = true
  []
[]

[Postprocessors]
  # --- Gating and growth averages ---
  [avg_fa]
    type       = ElementAverageValue
    variable   = fa
    execute_on = 'initial timestep_end'
  []

  [avg_gp]
    type       = ElementAverageValue
    variable   = gp
    execute_on = 'initial timestep_end'
  []

  [avg_gate_tot]
    type       = ElementAverageValue
    variable   = gate_tot
    execute_on = 'initial timestep_end'
  []

  [avg_ke]
    type       = ElementAverageValue
    variable   = ke
    execute_on = 'initial timestep_end'
  []

  [avg_eta]
    type       = ElementAverageValue
    variable   = eta
    execute_on = 'initial timestep_end'
  []

  # Product of averages <f_a>*<g_p> (for comparison with <f_a g_p>)
  [avg_fa_times_avg_gp]
    type       = ParsedPostprocessor
    pp_names   = 'avg_fa avg_gp'
    expression = 'avg_fa * avg_gp'
    execute_on = 'initial timestep_end'
  []

  # --- Volume & mechanics ---
  # Average J = det(F) over the reference domain; equals V(t)/V0 exactly.
  [avg_J]
    type       = ElementAverageValue
    variable   = volume_ratio
    execute_on = 'initial timestep_end'
  []

  [avg_press]
    type       = ElementAverageValue
    variable   = press
    execute_on = 'initial timestep_end'
  []

  [avg_n]
    type       = ElementAverageValue
    variable   = n
    execute_on = 'initial timestep_end'
  []

  [V]
    type              = VolumePostprocessor
    use_displaced_mesh = true    # <-- fix this name
    execute_on        = 'initial timestep_end'
  []

  [gate_min]
    type = ElementExtremeValue
    variable = gate_tot
    value_type = min
    execute_on = 'initial timestep_end'
  []

  [gate_max]
    type = ElementExtremeValue
    variable = gate_tot
    value_type = max
    execute_on = 'initial timestep_end'
  []

  [press_min]
    type = ElementExtremeValue
    variable = press
    value_type = min
    execute_on = 'initial timestep_end'
  []

  [press_max]
    type = ElementExtremeValue
    variable = press
    value_type = max
    execute_on = 'initial timestep_end'
  []


[]

[Outputs]
  exodus    = true
  csv       = true
  file_base = 2DExp2
  interval  = 10
[]
