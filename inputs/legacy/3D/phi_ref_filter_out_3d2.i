# phi_ref_filter_3d_cartcut.i
# Correlated random field for phi_ref_ic on a 3D puck mesh built from:
# Cartesian square -> mark disk -> delete outside -> extrude

phi0 = 0.10
Aphi = 0.08

lc   = 0.001

phi_min = ${fparse phi0 - Aphi}
phi_max = ${fparse phi0 + Aphi}

Dflt  = 1
t_end = ${fparse lc*lc / Dflt}

[Mesh]
  parallel_type = replicated

  # --- 2D base mesh: uniform Cartesian square covering the disk radius R=2 ---
  # Choose nx=ny=80 so dx=dy=4/80=0.05, giving ~O(5k) quads in the disk (similar scale to your annular case).
  [square]
    type      = GeneratedMeshGenerator
    dim       = 2
    nx        = 100
    ny        = 100
    xmin      = -2.0
    xmax      =  2.0
    ymin      = -2.0
    ymax      =  2.0
    show_info = false
  []

  # Mark elements whose centroids lie inside the disk as block 1 ("puck")
  [mark_disk]
    type                   = ParsedSubdomainMeshGenerator
    input                  = square
    combinatorial_geometry = 'x*x + y*y <= 4.0'   # R^2 = 2^2 = 4
    block_id               = 1
    block_name             = 'puck'
  []

  # Delete the outside (original block 0), and name the cut boundary "outer"
  [delete_outside]
    type         = BlockDeletionGenerator
    input        = mark_disk
    block        = '0'
    new_boundary = 'outer'
  []

  # --- Extrude to 3D puck ---
  [cyl]
    type            = AdvancedExtruderGenerator
    input           = delete_outside
    direction       = '0 0 1'
    heights         = '1.0'
    num_layers      = '30'
    bottom_boundary = 'bottom'
    top_boundary    = 'top'
  []

  # Pins (same as your file)
  [pin_center]
    type             = ExtraNodesetGenerator
    input            = cyl
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
    type      = BoundingValueNodalDamper
    variable  = phi_ref_ic
    min_value = 1e-8
    max_value = 0.9999
  []
[]

[Executioner]
  type       = Transient
  scheme     = bdf2
  solve_type = NEWTON

  dt       = ${fparse t_end/20}
  dtmin    = ${fparse t_end/200}
  dtmax    = ${fparse t_end/10}
  end_time = ${t_end}

  nl_rel_tol = 1e-10
  nl_abs_tol = 1e-12

  petsc_options_iname = '-ksp_type -pc_type -ksp_rtol'
  petsc_options_value = 'cg       hypre    1e-12'
[]

[Outputs]
  [out]
    type      = Exodus
    file_base = /home/amcgs/projects/ralphie/inputs/3D/phi_ref_filter3d_out
    interval  = 1
  []
[]
