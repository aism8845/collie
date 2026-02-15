# Strategy B fallback: small isotropic artificial nutrient diffusion
[Materials]
  [art]
    type = ADGenericConstantMaterial
    prop_names = 'D_art'
    prop_values = '1e-4'
  []
[]

[Kernels]
  [n_artdiff]
    type = ADMatDiffusion
    variable = n
    diffusivity = D_art
  []
[]

