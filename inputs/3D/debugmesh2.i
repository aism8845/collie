[Mesh]
  parallel_type = replicated

  # ---- Parameters (keep your R=2 and height=1) ----
  # Outer radius R=2, so R^2=4.

  # 1) Circular boundary as a curve mesh (closed loop)
  [outer_curve]
    type = ParsedCurveGenerator
    x_formula = 'R*cos(t)'
    y_formula = 'R*sin(t)'
    section_bounding_t_values = '0 ${fparse 2*pi}'
    constant_names = 'pi R'
    constant_expressions = '${fparse pi} 2.0'
    nums_segments = '80'          # increase for smoother rim
    is_closed_loop = true
  []

  # 2) Triangulate inside the curve
  [tri2d]
    type = XYDelaunayGenerator
    boundary = outer_curve
    desired_area = 0.005           # tune ~ element size; smaller = more elements
    refine_boundary = false
    output_subdomain_name = 'puck' # gives you a named block for materials
  []

  # 3) Explicitly tag the outer boundary as "outer"
  #    (helps downstream generators expecting an external boundary name)
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

  # 5) Option B: convert 2D elements to second order before extrusion (recommended)
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
    num_layers = '30'
    bottom_boundary = 'bottom'
    top_boundary = 'top'
  []

  # Pins (same idea as your file)
  [pin_center]
    type = ExtraNodesetGenerator
    input = cyl
    coord = '0 0 0.5'
    use_closest_node = true
    new_boundary = 'pin_center'
  []
  [pin_x]
    type = ExtraNodesetGenerator
    input = pin_center
    coord = '1 0 0.5'
    use_closest_node = true
    new_boundary = 'pin_x'
  []
  [pin_y]
    type = ExtraNodesetGenerator
    input = pin_x
    coord = '0 1 0.5'
    use_closest_node = true
    new_boundary = 'pin_y'
  []

  final_generator = pin_y
[]


[Executioner]
  type = Steady
[]

[Outputs]
  exodus = true
  file_base = simplemesh_
[]
