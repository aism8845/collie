# phi_ref_filter2d_RZ_generatedmesh.i
# Correlated random field (RandomIC + diffusion filter) on the SAME RZ GeneratedMesh you provided.

phi0 = 0.085
Aphi = 0.03
lc   = 0.001      # NOTE: choose lc >= a few * local element size to actually see "splotches"

phi_min = ${fparse phi0 - Aphi}
phi_max = ${fparse phi0 + Aphi}

Dflt  = 1
t_end = ${fparse lc*lc / Dflt}

[Mesh]
  coord_type    = RZ
  rz_coord_axis = Y
  parallel_type = distributed   # axial = y, radial = x (r=x, z=y)

  [gm]
    type = GeneratedMeshGenerator
    dim  = 2
    nx   = 140
    ny   = 40
    xmin = 0.0
    xmax = 2.0
    ymin = 0.0
    ymax = 0.5
    bias_y = 1.0
    bias_x = 1.0
  []
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
  solve_type = LINEAR    # diffusion filter is linear

  dt       = ${fparse t_end/20}
  dtmin    = ${fparse t_end/200}
  dtmax    = ${fparse t_end/10}
  end_time = ${t_end}

  petsc_options_iname  = '-ksp_type -pc_type -ksp_rtol'
  petsc_options_value  = 'cg hypre 1e-12'
[]

[Outputs]
  [out]
    type      = Exodus
    file_base = inputs/current/phi_ref_filter2d_out
    time_step_interval = 1
  []
[]
