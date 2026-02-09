#include "CellGelMixtureOpt.h"
#include "RankTwoTensor.h"
#include "RankFourTensor.h"
#include "MooseMesh.h"
#include <cmath>

registerMooseObject("collieApp", CellGelMixtureOpt);

InputParameters
CellGelMixtureOpt::validParams()
{
  InputParameters params = ComputeLagrangianStressPK2::validParams();
  params.addClassDescription("CellGelMixtureOpt: mixture cell+gel material with nutrient and stress gates.");

  params.addParam<Real>("phi_cell_0", 0.1, "Initial/reference cell volume fraction");
  params.addParam<Real>("G_cell", 1.0, "Cell shear modulus");
  params.addParam<Real>("G_gel", 1.0, "Gel shear modulus");
  params.addParam<Real>("k_n_cons", 1.0, "Nutrient consumption coefficient");
  params.addParam<Real>("D_nutrient", 1.0, "Nutrient diffusivity (reference) used in D scaling");

  params.addParam<Real>("k_exp_max", 1.0, "Max expansion/growth rate");
  params.addParam<Real>("c1", 1.0, "Tanh ramp steepness");
  params.addParam<Real>("c2", 0.0, "Tanh ramp midpoint time");

  params.addParam<Real>("k_T1_max", 1.0, "Max T1 rearrangement rate");
  params.addParam<Real>("chi_str", 1.0, "Strength of chi gate");
  params.addParam<Real>("press_str", 1.0, "Strength of pressure gate");
  params.addParam<Real>("beta_T1", 1.0, "Threshold/scale for T1 gate");
  params.addParam<Real>("m_T1", 2.0, "Hill exponent for T1 gate");
  params.addParam<Real>("q_cell", 1.0, "Exponent in cell volumetric penalty term");
  params.addParam<Real>("k_diss_0", 0.0, "Baseline matrix dissolution rate");

  params.addParam<Real>("n_c1", 2.0, "Hill exponent for nutrient activation gate");
  params.addParam<Real>("n_c2", 1.0, "Half-saturation for nutrient activation gate");

  params.addCoupledVar("n", "Nutrient scalar field (optional)");
  params.addParam<bool>("compute_dn_source_ref_dn", false, "Whether to compute dn_source_ref_dn material property");

  params.addParam<Real>("epsilon", 1e-8, "Epsilon for Miehe-style algorithmic tangent");
  params.addParam<Real>("fd_n_eps", 1e-6, "Finite-difference perturbation used to compute dP/dn (mechanics–nutrient coupling)");
  params.addParam<bool>("compute_dpk1_dn", true, "Compute and store the material property dpk1_dn = dP/dn (can be expensive)");

  return params;
}

CellGelMixtureOpt::CellGelMixtureOpt(const InputParameters & parameters)
  : ComputeLagrangianStressPK2(parameters),
    _phi_cell_0(getParam<Real>("phi_cell_0")),
    _G_cell(getParam<Real>("G_cell")),
    _G_gel(getParam<Real>("G_gel")),
    _k_n_cons(getParam<Real>("k_n_cons")),
    _D_nutrient(getParam<Real>("D_nutrient")),
    _k_exp_max(getParam<Real>("k_exp_max")),
    _c1(getParam<Real>("c1")),
    _c2(getParam<Real>("c2")),
    _k_T1_max(getParam<Real>("k_T1_max")),
    _chi_str(getParam<Real>("chi_str")),
    _press_str(getParam<Real>("press_str")),
    _beta_T1(getParam<Real>("beta_T1")),
    _m_T1(getParam<Real>("m_T1")),
    _q_cell(getParam<Real>("q_cell")),
    _k_diss_0(getParam<Real>("k_diss_0")),
    _has_n(isCoupled("n")),
    _n(_has_n ? coupledValue("n") : _zero),
    _dn_source_ref_dn(nullptr),
    _epsilon(getParam<Real>("epsilon")),
    _fd_n_eps(getParam<Real>("fd_n_eps")),
    _compute_dpk1_dn(getParam<bool>("compute_dpk1_dn")),
    _phi_cell_ref(declareProperty<Real>("phi_cell_ref")),
    _phi_cell_ref_old(getMaterialPropertyOld<Real>("phi_cell_ref")),
    _bE_cell(declareProperty<RankTwoTensor>("bE_cell")),
    _bE_cell_old(getMaterialPropertyOld<RankTwoTensor>("bE_cell")),
    _bE_pmat(declareProperty<RankTwoTensor>("bE_pmat")),
    _bE_pmat_old(getMaterialPropertyOld<RankTwoTensor>("bE_pmat")),
    _ke_old(declareProperty<Real>("ke_old")),
    _kT1_old(declareProperty<Real>("kT1_old")),
    _k_diss_old(declareProperty<Real>("k_diss_old")),
    _eta_old(declareProperty<Real>("eta_old")),
    _D_phys(declareProperty<Real>("D_phys")),
    _D_ref(declareProperty<Real>("D_ref")),
    _phi_cell(declareProperty<Real>("phi_cell")),
    _n_c1(getParam<Real>("n_c1")),
    _n_c2(getParam<Real>("n_c2")),
    _phi_ref_from_ic(declareProperty<Real>("phi_ref_from_ic")),
    _sigma_cell(declareProperty<RankTwoTensor>("sigma_cell")),
    _sigma_pmat(declareProperty<RankTwoTensor>("sigma_pmat")),
    _eta(declareProperty<Real>("eta")),
    _chi(declareProperty<Real>("chi")),
    _ke(declareProperty<Real>("ke")),
    _kT1(declareProperty<Real>("kT1")),
    _k_diss(declareProperty<Real>("k_diss")),
    _volume_ratio(declareProperty<Real>("volume_ratio")),
    _kh(declareProperty<Real>("kh")),
    _fa(declareProperty<Real>("fa")),
    _pressure(declareProperty<Real>("pressure")),
    _gp(declareProperty<Real>("gp")),
    _gate_tot(declareProperty<Real>("gate_tot")),
    _gamma_n_local(declareProperty<Real>("gamma_n_local")),
    _D_eff(declareProperty<RankTwoTensor>("D_eff")),
    _n_source_ref(declareProperty<Real>("n_source_ref")),
    _pk1_stress(declareProperty<RankTwoTensor>("pk1_stress")),
    _dpk1_dn(declareProperty<RankTwoTensor>("dpk1_dn")),
    _F_old_inv(getMaterialPropertyOldByName<RankTwoTensor>(_base_name + "deformation_gradient")),
    _S_old(getMaterialPropertyOldByName<RankTwoTensor>(_base_name + "stress_PK2")),
{
  if (getParam<bool>("compute_dn_source_ref_dn"))
    _dn_source_ref_dn = &declareProperty<Real>("dn_source_ref_dn");
}

void
CellGelMixtureOpt::initQpStatefulProperties()
{
  _phi_cell_ref[_qp] = _phi_cell_0;
  _phi_cell[_qp] = _phi_cell_0;

  _bE_cell[_qp].zero();
  _bE_cell[_qp].addIa(1.0);

  _bE_pmat[_qp].zero();
  _bE_pmat[_qp].addIa(1.0);

  _ke_old[_qp] = 0.0;
  _kT1_old[_qp] = 0.0;
  _k_diss_old[_qp] = 0.0;
  _eta_old[_qp] = 1.0;

  _phi_ref_from_ic[_qp] = -1.0;

  _pk1_stress[_qp].zero();
  _dpk1_dn[_qp].zero();
}

static inline Real
pow53(const Real x)
{
  // x^(5/3)
  return std::pow(x, 5.0 / 3.0);
}

void
CellGelMixtureOpt::computeQpProperties()
{
  const RankTwoTensor I(RankTwoTensor::initIdentity);

  // ---------------------------------------------------------------------------
  // Kinematics and deformation-rate measure
  // ---------------------------------------------------------------------------
  const Real inv_dt = (_dt > 0.0) ? (1.0 / _dt) : 0.0;
  const RankTwoTensor ell_old = inv_dt * (_F[_qp] * _F_old_inv[_qp] - I);

  const Real kdiv_old = ell_old.trace();
  _kh[_qp] = kdiv_old;

  const Real J_def = _F[_qp].det();
  _volume_ratio[_qp] = J_def;

  // ---------------------------------------------------------------------------
  // Mixture volume fraction (keep your original mapping)
  // ---------------------------------------------------------------------------
  const Real phi_ref_here = _phi_cell_ref[_qp];
  _phi_ref_from_ic[_qp] = phi_ref_here;

  const Real denom_phi = ((J_def - 1.0) * phi_ref_here + 1.0);
  _phi_cell[_qp] = (denom_phi != 0.0) ? ((J_def * phi_ref_here) / denom_phi) : phi_ref_here;

  // ---------------------------------------------------------------------------
  // Nutrient diffusivity scaling and effective diffusion tensor
  // ---------------------------------------------------------------------------
  const Real Dscale = std::max(1.0 - _phi_cell[_qp], 0.0) / (1.0 + _phi_cell[_qp] / 2.0);
  _D_phys[_qp] = Dscale * _D_nutrient;
  _D_ref[_qp]  = _D_phys[_qp] / (J_def + 1e-12);

  const RankTwoTensor Finv  = _F[_qp].inverse();
  const RankTwoTensor Dtensor = _D_phys[_qp] * I;
  _D_eff[_qp] = J_def * Finv * Dtensor * Finv.transpose();

  // ---------------------------------------------------------------------------
  // Smooth positive-part for nutrient (as in your original code)
  // ---------------------------------------------------------------------------
  const Real delta = 1e-5;

  Real n_pos = 0.0;
  Real d_npos_dn = 0.0;
  if (_has_n)
  {
    const Real n_here = _n[_qp];
    const Real root = std::sqrt(n_here * n_here + delta * delta);
    n_pos = 0.5 * (n_here + root);
    d_npos_dn = 0.5 * (1.0 + n_here / (root + 1e-12));
  }

  // Activation gate f_a(n)
  Real fa_here = 1.0;
  Real dfa_dn = 0.0;
  if (_has_n)
  {
    const Real num = std::pow(n_pos, _n_c1);
    const Real den = num + std::pow(_n_c2, _n_c1);
    fa_here = (den > 0.0) ? (num / den) : 0.0;

    if (den > 0.0)
    {
      const Real dnum_dn = (_n_c1 > 0.0) ? (_n_c1 * std::pow(n_pos, _n_c1 - 1.0) * d_npos_dn) : 0.0;
      dfa_dn = (dnum_dn * den - num * dnum_dn) / (den * den);
    }
  }
  _fa[_qp] = fa_here;

  // Local nutrient consumption factor and source term
  _gamma_n_local[_qp] = 10.0 * _k_n_cons * _phi_cell[_qp] / (_phi_cell_0 + 1e-12);
  _n_source_ref[_qp]  = -std::min(_gamma_n_local[_qp] * fa_here, 1.0);

  if (_dn_source_ref_dn && _has_n)
    (*_dn_source_ref_dn)[_qp] = -std::min(_gamma_n_local[_qp], 1.0) * dfa_dn;

  // ---------------------------------------------------------------------------
  // Stress gate g_p evaluated on OLD internal variables (bE_old) so that SAME-STEP
  // rates are well-defined without an implicit loop.
  // ---------------------------------------------------------------------------
  const Real phi_cell = _phi_cell[_qp];
  const Real k_pmat = 4.0 / 3.0 * _G_gel * std::max(1.0 - phi_cell, 0.0) / (phi_cell + 1e-12);

  const Real JE_cell_old = std::sqrt(std::max(_bE_cell_old[_qp].det(), 0.0));
  const RankTwoTensor sigma_cell_old =
      (_G_cell / (JE_cell_old + 1e-12)) * (_bE_cell_old[_qp] - std::pow(JE_cell_old, _q_cell) * I);

  const Real JE_pmat_old = std::sqrt(std::max(_bE_pmat_old[_qp].det(), 0.0));
  const RankTwoTensor sigma_pmat_old =
      (_G_gel / (pow53(JE_pmat_old) + 1e-12)) * _bE_pmat_old[_qp].deviatoric() +
      k_pmat * (JE_pmat_old - 1.0) * I;

  const RankTwoTensor sigma_old = phi_cell * sigma_cell_old + (1.0 - phi_cell) * sigma_pmat_old;
  const Real pressure_old = sigma_old.trace() / 3.0;
  const Real press_cell_old = pressure_old / (phi_cell + 1e-12);

  Real gp_here = 1.0;
  if (_press_str > 0.0)
    gp_here = 1.0 / (1.0 + std::exp(_press_str * (press_cell_old - 1.0)));
  _gp[_qp] = gp_here;

  const Real gate_tot_here = fa_here * gp_here;
  _gate_tot[_qp] = gate_tot_here;

  // ---------------------------------------------------------------------------
  // SAME-STEP rates (depend on SAME-STEP nutrient through fa_here)
  // ---------------------------------------------------------------------------
  Real k_exp0 = 0.0;
  if (_t < _c2)
    k_exp0 = 0.0;
  else
    k_exp0 = _k_exp_max * 0.5 * (std::tanh(_c1 * (_t - _c2)) + 1.0);

  _ke[_qp] = k_exp0 * gate_tot_here;

  // chi and T1 rate (keep your original logic)
  _chi[_qp] = _chi_str * (gate_tot_here - 0.5);

  if (_beta_T1 > 0.0)
  {
    const Real chi_here = std::max(_chi[_qp], 0.0);
    const Real chi_crit = _beta_T1;
    const Real m = _m_T1;

    Real gate_T1 = 0.0;
    if (chi_here > 0.0 && chi_crit > 0.0)
    {
      const Real num = std::pow(chi_here, m);
      const Real den = num + std::pow(chi_crit, m);
      gate_T1 = (den > 0.0) ? num / den : 0.0;
    }

    _kT1[_qp] = _k_T1_max * gate_T1;
  }
  else
    _kT1[_qp] = _k_T1_max / (1.0 + std::exp(-1.0 * _beta_T1 * (_chi[_qp] - _chi_str)));

  // Matrix dissolution depends on nutrient availability (fa)
  _k_diss[_qp] = _k_diss_0 * fa_here;

  // ---------------------------------------------------------------------------
  // Update internal variables bE using SAME-STEP rates (ke, kT1, k_diss)
  // ---------------------------------------------------------------------------
  {
    const RankTwoTensor bE_old = _bE_cell_old[_qp];
    const RankTwoTensor bE_old_inv = bE_old.inverse();
    const Real bE_old_inv_trace = bE_old_inv.trace();

    const RankTwoTensor bE_cell_dot =
        ell_old * bE_old + bE_old * ell_old.transpose() -
        (2.0 / 3.0) * _ke[_qp] * bE_old -
        ((4.0 / 3.0) * kdiv_old + _kT1[_qp]) *
            (bE_old - (3.0 / (bE_old_inv_trace + 1e-12)) * I);

    _bE_cell[_qp] = bE_old + bE_cell_dot * _dt;
  }

  {
    const RankTwoTensor bE_old = _bE_pmat_old[_qp];
    const RankTwoTensor bE_old_inv = bE_old.inverse();
    const Real bE_old_inv_trace = bE_old_inv.trace();
    const Real JE_pmat = std::sqrt(std::max(bE_old.det(), 0.0));

    const RankTwoTensor bE_pmat_dot =
        ell_old * bE_old + bE_old * ell_old.transpose() -
        _k_diss[_qp] * (bE_old - (3.0 / ((bE_old_inv_trace + 1e-12) * (JE_pmat + 1e-12))) * I);

    _bE_pmat[_qp] = bE_old + bE_pmat_dot * _dt;
  }

  // ---------------------------------------------------------------------------
  // Stresses from UPDATED internal variables
  // ---------------------------------------------------------------------------
  const Real JE_cell = std::sqrt(std::max(_bE_cell[_qp].det(), 0.0));
  _sigma_cell[_qp] =
      (_G_cell / (JE_cell + 1e-12)) * (_bE_cell[_qp] - std::pow(JE_cell, _q_cell) * I);

  const Real JE_pmat = std::sqrt(std::max(_bE_pmat[_qp].det(), 0.0));
  _sigma_pmat[_qp] =
      (_G_gel / (pow53(JE_pmat) + 1e-12)) * _bE_pmat[_qp].deviatoric() +
      k_pmat * (JE_pmat - 1.0) * I;

  // eta update (now uses SAME-STEP k_diss)
  _eta[_qp] = (phi_ref_here / (J_def * _phi_cell_0 + 1e-12)) * std::exp(-_k_diss[_qp] * _dt);

  // ---------------------------------------------------------------------------
  // PK1 stress and dP/dn (finite difference)
  // ---------------------------------------------------------------------------
  const RankTwoTensor sigma_cauchy = phi_cell * _sigma_cell[_qp] + (1.0 - phi_cell) * _sigma_pmat[_qp];
  const RankTwoTensor FinvT = Finv.transpose();
  const RankTwoTensor P0 = J_def * sigma_cauchy * FinvT;
  _pk1_stress[_qp] = P0;

  _dpk1_dn[_qp].zero();

  if (_compute_dpk1_dn && _has_n)
  {
    const Real n_here = _n[_qp];
    const Real eps = std::max(_fd_n_eps, 1e-12);
    const Real n_pert = n_here + eps;

    const Real rootp = std::sqrt(n_pert * n_pert + delta * delta);
    const Real n_pos_p = 0.5 * (n_pert + rootp);
    const Real nump = std::pow(n_pos_p, _n_c1);
    const Real denp = nump + std::pow(_n_c2, _n_c1);
    const Real fa_p = (denp > 0.0) ? (nump / denp) : 0.0;

    const Real gate_tot_p = fa_p * gp_here;

    const Real ke_p = k_exp0 * gate_tot_p;

    const Real chi_p = _chi_str * (gate_tot_p - 0.5);

    Real kT1_p = _k_T1_max;
    if (_beta_T1 > 0.0)
    {
      const Real chi_here_p = std::max(chi_p, 0.0);
      const Real chi_crit = _beta_T1;
      const Real m = _m_T1;

      Real gate_T1 = 0.0;
      if (chi_here_p > 0.0 && chi_crit > 0.0)
      {
        const Real num = std::pow(chi_here_p, m);
        const Real den = num + std::pow(chi_crit, m);
        gate_T1 = (den > 0.0) ? num / den : 0.0;
      }
      kT1_p = _k_T1_max * gate_T1;
    }
    else
      kT1_p = _k_T1_max / (1.0 + std::exp(-1.0 * _beta_T1 * (chi_p - _chi_str)));

    const Real k_diss_p = _k_diss_0 * fa_p;

    // cell bE update
    RankTwoTensor bE_cell_p;
    {
      const RankTwoTensor bE_old = _bE_cell_old[_qp];
      const RankTwoTensor bE_old_inv = bE_old.inverse();
      const Real tr_inv = bE_old_inv.trace();

      const RankTwoTensor bE_dot =
          ell_old * bE_old + bE_old * ell_old.transpose() -
          (2.0 / 3.0) * ke_p * bE_old -
          ((4.0 / 3.0) * kdiv_old + kT1_p) * (bE_old - (3.0 / (tr_inv + 1e-12)) * I);

      bE_cell_p = bE_old + bE_dot * _dt;
    }

    // pmat bE update
    RankTwoTensor bE_pmat_p;
    {
      const RankTwoTensor bE_old = _bE_pmat_old[_qp];
      const RankTwoTensor bE_old_inv = bE_old.inverse();
      const Real tr_inv = bE_old_inv.trace();
      const Real JE_old = std::sqrt(std::max(bE_old.det(), 0.0));

      const RankTwoTensor bE_dot =
          ell_old * bE_old + bE_old * ell_old.transpose() -
          k_diss_p * (bE_old - (3.0 / ((tr_inv + 1e-12) * (JE_old + 1e-12))) * I);

      bE_pmat_p = bE_old + bE_dot * _dt;
    }

    const Real JEc_p = std::sqrt(std::max(bE_cell_p.det(), 0.0));
    const RankTwoTensor sigma_cell_p =
        (_G_cell / (JEc_p + 1e-12)) * (bE_cell_p - std::pow(JEc_p, _q_cell) * I);

    const Real JEp_p = std::sqrt(std::max(bE_pmat_p.det(), 0.0));
    const RankTwoTensor sigma_pmat_p =
        (_G_gel / (pow53(JEp_p) + 1e-12)) * bE_pmat_p.deviatoric() +
        k_pmat * (JEp_p - 1.0) * I;

    const RankTwoTensor sigma_p = phi_cell * sigma_cell_p + (1.0 - phi_cell) * sigma_pmat_p;
    const RankTwoTensor P1 = J_def * sigma_p * FinvT;

    _dpk1_dn[_qp] = (P1 - P0) / eps;
  }

  // Keep base-class bookkeeping (if any) intact
  ComputeLagrangianStressPK2::computeQpProperties();
}

//----------Cauchy stress
void
CellGelMixtureOpt::computeQpCauchyStress()
{
  // Mix Cauchy-like stresses computed in computeQpProperties into _cauchy_stress
  _cauchy_stress[_qp] = _phi_cell[_qp] * _sigma_cell[_qp] + (1.0 - _phi_cell[_qp]) * _sigma_pmat[_qp];
  _pressure[_qp] = _cauchy_stress[_qp].trace() / 3.0;
}

//----------PK2 stress
void
CellGelMixtureOpt::computeQpPK2Stress()
{
  // (Existing FD tangent implementation; patched to use current rates rather than *_old)
  RankTwoTensor stress_old = _S_old[_qp];
  RankFourTensor jacobian;
  jacobian.zero();

  for (unsigned int k = 0; k < 9; ++k)
  {
    const unsigned int i = k / 3;
    const unsigned int j = k % 3;

    RankTwoTensor F_p = _F[_qp];
    F_p(i, j) += _epsilon;

    const Real inv_dt = (_dt > 0.0) ? (1.0 / _dt) : 0.0;
    const RankTwoTensor I(RankTwoTensor::initIdentity);
    const RankTwoTensor ell_trial = inv_dt * (F_p * _F_old_inv[_qp] - I);

    const Real kh_trial = ell_trial.trace();
    const Real J_def_trial = F_p.det();

    const Real phi_ref = _phi_cell_ref[_qp];
    const Real denom_phi = ((J_def_trial - 1.0) * phi_ref + 1.0);
    const Real phi_cell_trial = (denom_phi != 0.0) ? ((J_def_trial * phi_ref) / denom_phi) : phi_ref;

    const Real k_pmat_trial =
        4.0 / 3.0 * _G_gel * std::max(1.0 - phi_cell_trial, 0.0) / (phi_cell_trial + 1e-12);

    // cell bE trial
    RankTwoTensor bE_cell_trial;
    {
      const RankTwoTensor bE_old = _bE_cell_old[_qp];
      const RankTwoTensor bE_old_inv = bE_old.inverse();
      const Real tr_inv = bE_old_inv.trace();

      const RankTwoTensor bE_dot =
          ell_trial * bE_old + bE_old * ell_trial.transpose() -
          (2.0 / 3.0) * _ke[_qp] * bE_old -
          ((4.0 / 3.0) * kh_trial + _kT1[_qp]) * (bE_old - (3.0 / (tr_inv + 1e-12)) * I);

      bE_cell_trial = bE_old + bE_dot * _dt;
    }

    // pmat bE trial
    RankTwoTensor bE_pmat_trial;
    {
      const RankTwoTensor bE_old = _bE_pmat_old[_qp];
      const RankTwoTensor bE_old_inv = bE_old.inverse();
      const Real tr_inv = bE_old_inv.trace();
      const Real JE_old = std::sqrt(std::max(bE_old.det(), 0.0));

      const RankTwoTensor bE_dot =
          ell_trial * bE_old + bE_old * ell_trial.transpose() -
          _k_diss[_qp] * (bE_old - (3.0 / ((tr_inv + 1e-12) * (JE_old + 1e-12))) * I);

      bE_pmat_trial = bE_old + bE_dot * _dt;
    }

    // stresses
    const Real JE_cell_trial = std::sqrt(std::max(bE_cell_trial.det(), 0.0));
    const RankTwoTensor sigma_cell_trial =
        (_G_cell / (JE_cell_trial + 1e-12)) * (bE_cell_trial - std::pow(JE_cell_trial, _q_cell) * I);

    const Real JE_pmat_trial = std::sqrt(std::max(bE_pmat_trial.det(), 0.0));
    const RankTwoTensor sigma_pmat_trial =
        (_G_gel / (pow53(JE_pmat_trial) + 1e-12)) * bE_pmat_trial.deviatoric() +
        k_pmat_trial * (JE_pmat_trial - 1.0) * I;

    const RankTwoTensor sigma_cauchy_trial =
        phi_cell_trial * sigma_cell_trial + (1.0 - phi_cell_trial) * sigma_pmat_trial;

    // PK2 from sigma: S = J * F^{-1} sigma F^{-T}
    const RankTwoTensor Finv = F_p.inverse();
    const RankTwoTensor S_trial = J_def_trial * Finv * sigma_cauchy_trial * Finv.transpose();

    // FD derivative dS/dF_{ij}
    const RankTwoTensor dS = (S_trial - stress_old) / _epsilon;

    for (unsigned int m = 0; m < 3; ++m)
      for (unsigned int n = 0; n < 3; ++n)
        jacobian(m, n, i, j) = dS(m, n);
  }

  _S[_qp] = stress_old;
  _elasticity_tensor[_qp] = jacobian;
}
