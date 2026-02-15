# Enable Tier-1 SPD projection for bE tensors in CellGelMixtureOpt.
# Defaults in C++ are off; this override is for barrier-window diagnostics only.
Materials/cell_gel/be_spd_project := true
Materials/cell_gel/be_eig_floor := 1e-12
Materials/cell_gel/be_project_mode := both
