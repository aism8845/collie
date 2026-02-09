# phi_ref_filter_3d.i
# Correlated random field for phi_ref_ic on the SAME 3D cylinder mesh as YP3D.i
# Correlation length: lc ~ sqrt(Dflt * t_end)

phi0 = 0.10
Aphi = 0.1

# IMPORTANT:
# Your current YP3D mesh is fairly coarse (nr=10, nt=10, layers=8), so lc=0.01 (from 2D)
# is far below element size and will not visibly change anything.
# Start with lc ~ O(element size). For your mesh, a reasonable first pass is 0.1–0.3.
lc   = 0.001

phi_min = ${fparse phi0 - Aphi}
phi_max = ${fparse phi0 + Aphi}

Dflt  = 1
t_end = ${fparse lc*lc / Dflt}

[Mesh]
  parallel_type = replicated

  # Same as YP3D.i
  [disk]
    type                 = AnnularMeshGenerator
    rmin                 = 0.0
    rmax                 = 2.0
    nr                   = 60
    nt                   = 80
    boundary_name_prefix = 'disk_'
    show_info            = false
  []

  [cyl]
    type            = AdvancedExtruderGenerator
    input           = disk
    direction       = '0 0 1'
    heights         = '1.0'
    num_layers      = '60'
    bottom_boundary = 'bottom'
    top_boundary    = 'top'
  []

  # Optional, but keeping these makes the mesh match YP3D more closely
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
