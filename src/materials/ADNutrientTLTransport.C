#include "ADNutrientTLTransport.h"

#include "MathUtils.h"
#include "MooseEnum.h"

#include <cmath>

registerMooseObject("collieApp", ADNutrientTLTransport);

namespace
{
template <typename T>
T
smooth_max(const T & x, const T & a, const Real eps)
{
  return 0.5 * (x + a + std::sqrt((x - a) * (x - a) + eps * eps));
}

template <typename T>
T
smooth_min(const T & x, const T & b, const Real eps)
{
  return 0.5 * (x + b - std::sqrt((x - b) * (x - b) + eps * eps));
}

template <typename T>
T
smooth_clamp(const T & x, const T & a, const T & b, const Real eps)
{
  return smooth_min(smooth_max(x, a, eps), b, eps);
}
} // namespace

InputParameters
ADNutrientTLTransport::validParams()
{
  InputParameters params = ADMaterial::validParams();

  params.addRequiredCoupledVar("n", "Nutrient variable.");

  // In RZ with SolidMechanics new system: ux=u_r, uy=u_z typically
  params.addRequiredCoupledVar("disp_r", "Radial displacement variable (u_r).");
  params.addRequiredCoupledVar("disp_z", "Axial displacement variable (u_z).");

  params.addParam<bool>("axisymmetric", true,
                        "If true, include hoop stretch F_theta_theta = 1 + u_r/r.");
  params.addParam<unsigned int>("radial_coord", 0,
                                "Which coordinate index is the radial coordinate r (0->x, 1->y, 2->z).");
  params.addParam<Real>("r_eps", 1e-12, "Small radius cutoff to avoid u_r/r singularity.");

  params.addRequiredParam<MaterialPropertyName>(
      "phi_cell_ref",
      "Name of the (non-AD) stateful material property for referential cell volume fraction.");

  params.addRequiredParam<Real>("D0", "Reference diffusivity.");
  params.addParam<Real>("D_floor", 1e-12, "Lower bound for diffusivity.");

  params.addRequiredParam<Real>("gamma_n0", "Nutrient consumption prefactor.");
  params.addRequiredParam<Real>("phi_max", "Maximum cell volume fraction.");
  MooseEnum crowding_model("maxwell bruggeman percolation", "maxwell");
  params.addParam<MooseEnum>("crowding_model", crowding_model,
                             "Crowding model for physical diffusivity.");
  params.addParam<Real>("crowd_exp", 2.0, "Exponent for bruggeman/percolation crowding models.");
  params.addRequiredParam<Real>("n_c1", "Half-activation concentration for f_a(n).");
  params.addRequiredParam<unsigned int>("n_c2", "Hill exponent for f_a(n) (integer).");
  params.addParam<Real>("smooth_eps_c", 1e-12, "Smoothing epsilon for concentration/phi clamps.");
  params.addParam<Real>("smooth_eps_D", 1e-12, "Smoothing epsilon for diffusivity floor.");

  params.addClassDescription(
      "AD TL nutrient transport coefficients with full coupling to displacements (builds F from AD displacement gradients).");

  return params;
}

ADNutrientTLTransport::ADNutrientTLTransport(const InputParameters & parameters)
  : ADMaterial(parameters),
    _n(adCoupledValue("n")),

    _ur(adCoupledValue("disp_r")),
    _uz(adCoupledValue("disp_z")),
    _grad_ur(adCoupledGradient("disp_r")),
    _grad_uz(adCoupledGradient("disp_z")),

    _ur_old(coupledValueOld("disp_r")),
    _uz_old(coupledValueOld("disp_z")),
    _grad_ur_old(coupledGradientOld("disp_r")),
    _grad_uz_old(coupledGradientOld("disp_z")),

    _ur_older(coupledValueOlder("disp_r")),
    _uz_older(coupledValueOlder("disp_z")),
    _grad_ur_older(coupledGradientOlder("disp_r")),
    _grad_uz_older(coupledGradientOlder("disp_z")),

    _axisymmetric(getParam<bool>("axisymmetric")),
    _radial_coord(getParam<unsigned int>("radial_coord")),
    _r_eps(getParam<Real>("r_eps")),
    _axial_coord((_radial_coord == 0) ? 1 : 0),

    _phi_cell_ref(getMaterialProperty<Real>(getParam<MaterialPropertyName>("phi_cell_ref"))),

    _D0(getParam<Real>("D0")),
    _D_floor(getParam<Real>("D_floor")),
    _gamma_n0(getParam<Real>("gamma_n0")),
    _phi_max(getParam<Real>("phi_max")),
    _crowding_model(getParam<MooseEnum>("crowding_model")),
    _crowd_exp(getParam<Real>("crowd_exp")),
    _n_c1(getParam<Real>("n_c1")),
    _n_c2(getParam<unsigned int>("n_c2")),
    _n_c1_pow(std::pow(_n_c1, static_cast<Real>(_n_c2))),
    _smooth_eps_c(getParam<Real>("smooth_eps_c")),
    _smooth_eps_D(getParam<Real>("smooth_eps_D")),

    _J_nutr(declareADProperty<Real>("J_nutr")),
    _Jdot_nutr(declareADProperty<Real>("Jdot_nutr")),
    _D_eff_nutr(declareADProperty<RealTensorValue>("D_eff_nutr")),
    _n_source_ref_nutr(declareADProperty<Real>("n_source_ref_nutr")),
    _jdot_rate_nutr(declareADProperty<Real>("jdot_rate_nutr")),

    _D_phys_nutr(declareADProperty<Real>("D_phys_nutr")),
    _D_ref_nutr(declareADProperty<Real>("D_ref_nutr")),
    _D_iso_nutr(declareADProperty<Real>("D_iso_nutr"))
{
  if (_radial_coord > 2)
    mooseError("ADNutrientTLTransport: 'radial_coord' must be 0, 1, or 2.");
}

void
ADNutrientTLTransport::computeQpProperties()
{
  const Real r = _q_point[_qp](_radial_coord);

  // --- Current gradients (AD) ---
  const auto & gur = _grad_ur[_qp];
  const auto & guz = _grad_uz[_qp];

  const ADReal dur_dr = gur(_radial_coord);
  const ADReal dur_dz = gur(_axial_coord);

  const ADReal duz_dr = guz(_radial_coord);
  const ADReal duz_dz = guz(_axial_coord);

  // --- Build F in cylindrical RZ basis: [r,z,theta] ---
  ADRankTwoTensor F;
  F.zero();

  // r-z block
  F(0,0) = 1.0 + dur_dr;  // F_rr
  F(0,1) =        dur_dz; // F_rz
  F(1,0) =        duz_dr; // F_zr
  F(1,1) = 1.0 + duz_dz;  // F_zz

  // theta-theta
  F(2,2) = 1.0;
  if (_axisymmetric)
  {
    const ADReal ur_here = _ur[_qp];
    if (r > _r_eps)
      F(2,2) = 1.0 + ur_here / r;
    else
      // on axis with u_r=0 BC: limit u_r/r -> du_r/dr
      F(2,2) = 1.0 + dur_dr;
  }

  const ADReal J = F.det();
  _J_nutr[_qp] = J;

  // --- Helper: J from OLD states (Real) ---
  auto J_from_old = [&](const VariableValue & urv,
                        const VariableValue & uzv,
                        const VariableGradient & gurv,
                        const VariableGradient & guzv) -> Real
  {
    const auto & gur_old = gurv[_qp];
    const auto & guz_old = guzv[_qp];

    const Real dur_dr_o = gur_old(_radial_coord);
    const Real dur_dz_o = gur_old(_axial_coord);
    const Real duz_dr_o = guz_old(_radial_coord);
    const Real duz_dz_o = guz_old(_axial_coord);

    RankTwoTensor Fold;
    Fold.zero();

    Fold(0,0) = 1.0 + dur_dr_o;
    Fold(0,1) =        dur_dz_o;
    Fold(1,0) =        duz_dr_o;
    Fold(1,1) = 1.0 + duz_dz_o;

    Fold(2,2) = 1.0;
    if (_axisymmetric)
    {
      const Real ur_o = urv[_qp];
      if (r > _r_eps)
        Fold(2,2) = 1.0 + ur_o / r;
      else
        Fold(2,2) = 1.0 + dur_dr_o;
    }

    return Fold.det();
  };

  const Real J_old   = (_t_step >= 1) ? J_from_old(_ur_old,   _uz_old,   _grad_ur_old,   _grad_uz_old)   : 1.0;
  const Real J_older = (_t_step >= 2) ? J_from_old(_ur_older, _uz_older, _grad_ur_older, _grad_uz_older) : 1.0;

  // --- Jdot: BDF1 at step 1; variable-step BDF2 afterward ---
  ADReal Jdot = 0.0;
  if (_t_step >= 1)
  {
    if (_t_step == 1 || _dt_old <= 0.0)
      Jdot = (J - J_old) / _dt;
    else
    {
      const Real rdt = _dt / _dt_old; // dt_n / dt_{n-1}
      const Real a0  = (1.0 + 2.0 * rdt) / (1.0 + rdt);
      const Real a1  = -(1.0 + rdt);
      const Real a2  = (rdt * rdt) / (1.0 + rdt);

      Jdot = (a0 * J + a1 * J_old + a2 * J_older) / _dt;
    }
  }

  _Jdot_nutr[_qp] = Jdot;

  // ADMatReaction adds (-rate*n). To add +Jdot*n, set rate = -Jdot.
  _jdot_rate_nutr[_qp] = -Jdot;

  // --- phi(J,phi_ref) mapping ---
  const Real phi_ref = _phi_cell_ref[_qp];
  const ADReal denom = (J - 1.0) * phi_ref + 1.0;
  ADReal phi = (J * phi_ref) / denom;

  phi = smooth_clamp(phi, ADReal(0.0), ADReal(_phi_max), _smooth_eps_c);

  // --- D_phys(phi) ---
  ADReal D_phys = _D0;
  if (_crowding_model == "maxwell")
    D_phys = _D0 * (1.0 - phi) / (1.0 + 0.5 * phi);
  else if (_crowding_model == "bruggeman")
  {
    const ADReal one_minus_phi = smooth_clamp(ADReal(1.0) - phi, ADReal(0.0), ADReal(1.0), _smooth_eps_c);
    D_phys = _D0 * std::pow(one_minus_phi, _crowd_exp);
  }
  else if (_crowding_model == "percolation")
  {
    const ADReal phi_max = std::max(_phi_max, 1e-16);
    const ADReal s = smooth_clamp(ADReal(1.0) - phi / phi_max, ADReal(0.0), ADReal(1.0), _smooth_eps_c);
    D_phys = _D0 * std::pow(s, _crowd_exp);
  }
  else
    mooseError("Unsupported crowding_model: ", _crowding_model);

  D_phys = smooth_max(D_phys, ADReal(_D_floor), _smooth_eps_D);

  // diagnostics (unique names)
  _D_phys_nutr[_qp] = D_phys;
  _D_ref_nutr[_qp]  = J * D_phys;

  // --- K = J * D_phys * F^{-1} * F^{-T} ---
  const ADRankTwoTensor Finv = F.inverse();
  const ADRankTwoTensor A    = Finv * Finv.transpose();

  ADRealTensorValue K;
  for (unsigned int i = 0; i < 3; ++i)
    for (unsigned int j = 0; j < 3; ++j)
      K(i,j) = J * D_phys * A(i,j);

  _D_eff_nutr[_qp] = K;

  _D_iso_nutr[_qp] = (K(0,0) + K(1,1) + K(2,2)) / 3.0;

  // --- f_a(n): Hill gate ---
  ADReal nloc = _n[_qp];
  nloc = smooth_max(nloc, ADReal(0.0), _smooth_eps_c);

  const ADReal num  = MathUtils::pow(nloc, static_cast<int>(_n_c2));
  const ADReal den2 = num + ADReal(_n_c1_pow) + 1e-16;
  const ADReal fa   = num / den2;

  const ADReal gamma_local = _gamma_n0 * phi * fa;

  // ADMatReaction adds (-rate*n); to add +J*gamma*n, set rate = -J*gamma
  _n_source_ref_nutr[_qp] = -J * gamma_local;
}
