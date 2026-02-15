# 3-way Schur fieldsplit with GAMG fallback for nutrient block.
Preconditioning/active := 'pc_fsp3_schur_gamg'

Executioner/petsc_options := '-snes_converged_reason -ksp_converged_reason'
Executioner/petsc_options_iname := '-snes_type -snes_linesearch_type -ksp_type -ksp_rtol -ksp_max_it -pc_type -pc_fieldsplit_type'
Executioner/petsc_options_value := 'newtonls bt fgmres 1e-6 200 fieldsplit schur'
