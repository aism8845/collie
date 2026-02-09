# phi_ref_filter3d_debugmesh2.i
# Correlated random field filter on the SAME mesh topology as YP3D4_debugmesh2_nullspace.i


phi0 = 0.10
Aphi = 0.03
lc   = 0.001

phi_min = ${fparse phi0 - Aphi}
phi_max = ${fparse phi0 + Aphi}

Dflt  = 1
t_end = ${fparse lc*lc / Dflt}

[Mesh]
  parallel_type = replicated

  # 1) Circular boundary as a curve mesh (closed loop)
  [outer_curve]
    type = ParsedCurveGenerator
    x_formula = 'R*cos(t)'
    y_formula = 'R*sin(t)'
    section_bounding_t_values = '0 ${fparse 2*pi}'
    constant_names = 'pi R'
    constant_expressions = '${fparse pi} 2.0'
    nums_segments = '150'
    is_closed_loop = true
  []

  # 2) Triangulate inside the curve
  [tri2d]
    type = XYDelaunayGenerator
    boundary = outer_curve
    desired_area = 0.01
    refine_boundary = false
    output_subdomain_name = 'puck'
  []

  # 3) Tag the outer boundary as "outer" (robust ring test; avoids abs() parser issues)
  [outer_sideset]
    type = ParsedGenerateSideset
    input = tri2d

    combinatorial_geometry = 'abs(x*x + y*y - R2) < tol'
    constant_names = 'R2 tol'
    constant_expressions = '4.0 0.08'

    new_sideset_name = 'outer'
    include_only_external_sides = true
  []

  # 4) Optional: boundary correction for circular fidelity/area preservation
  [circ_fix]
    type = CircularBoundaryCorrectionGenerator
    input = outer_sideset
    input_mesh_circular_boundaries = 'outer'
    custom_circular_tolerance = 1e-8
  []

  # 5) Option B: convert 2D elements to second order before extrusion
  [p2_2d]
    type = ElementOrderConversionGenerator
    input = circ_fix
    conversion_type = SECOND_ORDER
  []

  # 6) Extrude to 3D puck (triangles -> prisms)
  [cyl]
    type = AdvancedExtruderGenerator
    input = p2_2d
    direction = '0 0 1'
    heights = '1.0'
    num_layers = '40'
    bottom_boundary = 'bottom'
    top_boundary = 'top'
  []

  final_generator = cyl
[]


[Variables]
  [phi_ref_ic]
    family = LAGRANGE
    order  = FIRST
  []
[]

[ICs]
  [phi_ref_rand]
    type     = RandomIC
    variable = phi_ref_ic
    min      = ${phi_min}
    max      = ${phi_max}
    seed     = 12345
  []
[]

[Kernels]
  [td]
    type     = TimeDerivative
    variable = phi_ref_ic
  []

  [diff]
    type     = Diffusion
    variable = phi_ref_ic
  []
[]

[Dampers]
  [bounds]
    type     = BoundingValueNodalDamper
    variable = phi_ref_ic
    min_value = 1e-8
    max_value = 0.9999
  []
[]

[Executioner]
  type        = Transient
  scheme      = bdf2
  solve_type  = NEWTON

  dt       = ${fparse t_end/20}
  dtmin    = ${fparse t_end/200}
  dtmax    = ${fparse t_end/10}
  end_time = ${t_end}

  nl_rel_tol  = 1e-10
  nl_abs_tol  = 1e-12

  petsc_options_iname = '-ksp_type -pc_type -ksp_rtol'
  petsc_options_value = 'cg hypre 1e-12'
[]

[Outputs]
  [out]
    type      = Exodus
    file_base = phi_ref_filter3d_debugmesh4_out
    time_step_interval  = 1
  []
[]
