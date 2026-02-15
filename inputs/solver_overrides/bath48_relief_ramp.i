# Optional late-time stress-relief tweak for debugging the detF barrier.
# Keeps baseline unchanged unless this override is included.
# Note: RealFunctionControl on this material parameter is not enabled in this app,
# so this override uses a small constant k_T1_max during the restart window.

Materials/cell_gel/enable_T1 := true
Materials/cell_gel/chi_T1_smooth := 0.02
Materials/cell_gel/k_T1_max := 0.03
