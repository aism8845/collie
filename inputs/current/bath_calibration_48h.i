# -----------------------------------------------------------------------------
# bath_calibration_48h.i
# -----------------------------------------------------------------------------
# 48-hour wrapper around bath_calibration.i with a smooth 24-hour media refresh.
#
# Scaling / normalizers (experiment-normalized nondimensionalization):
#   t0 = 1 hour        => 1 sim second = 1 hour
#   L0 = 1 mm          => mesh coordinates are in mm
#   C0 = 20 mg/mL      => ~111 mM glucose, n = C/C0, bath = C_bath/C0
#   G0 = 0.7 MPa       => p* = p_phys / G0
#
# Geometry mapping:
#   R_phys = 2 mm  -> xmax = 2.0
#   H_phys = 1 mm  -> half-height model uses ymax = 0.5 (midplane symmetry)
#
# Bath accounting (half-height symmetry intent):
#   Use V_bath = 5000 (effective half-bath moles) with bath BC on 'top right'.
#
# Refresh model:
#   k_refill(t) = k_max * [s_up(t) - s_down(t)]
#   s_up   = 0.5*(1+tanh((t-24)/w))
#   s_down = 0.5*(1+tanh((t-(24+dur))/w))
#   defaults: k_max=500 1/h, w=0.02 h, dur=0.10 h.
!include bath_calibration.i

Executioner/end_time := 48.0

ScalarKernels/bath_mb/V_bath := 5000
ScalarKernels/bath_mb/n_feed := 1.0
ScalarKernels/bath_mb/k_refill := k_refill_pulse

Outputs/solver_watch/file_base := outputs/bath48_solver_watch
Outputs/mesh_watch/file_base := outputs/bath48_mesh_watch
Outputs/chk/file_base := outputs/chk/bath48
# Checkpoint output schedule is inherited from bath_calibration.i [Outputs/chk].

[Functions]
  [k_refill_pulse]
    type = ParsedFunction
    expression = '500.0*(0.5*(1+tanh((t-24.0)/0.02)) - 0.5*(1+tanh((t-24.1)/0.02)))'
  []
[]
