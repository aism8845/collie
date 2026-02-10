#include "CellGelMixtureOpt.h"
#include "MooseMesh.h"
#include "MooseVariable.h"
#include "libmesh/utility.h"
#include "TransientInterface.h"




#include <cmath>
#include <limits>
#include <algorithm>

using namespace std;

registerMooseObject("collieApp", CellGelMixtureOpt);



namespace
{
  /// Fast exact x^(5/3) for positive x (used in porous-matrix stress).
  inline Real pow53(const Real x)
  {
    const Real t = std::cbrt(x);      // x^(1/3)
    const Real t2 = t * t;
    const Real t4 = t2 * t2;
    return t4 * t;                    // t^5 = x^(5/3)
  }
} // namespace


// ---------------- params ----------------
InputParameters
CellGelMixtureOpt::validParams()
{
  InputParameters params = DerivativeMaterialInterface<ComputeLagrangianStressPK2>::validParams();
  params.addClassDescription("CellGelMixture Material model: Implements a TNT-based biphasic multiscale model for cell-laden gel");

  // Required
  params.addRequiredParam<Real>("G_cell", "Shear modulus of cell phase");
  params.addRequiredParam<Real>("q_cell", "q = -2nu/(1-2nu) for cell phase, where nu is Poisson ratio");
  params.addRequiredParam<Real>("k_exp_max","Max isotropic expansion rate after ramping (sets time scale)");
  params.addRequiredParam<Real>("c1", "Expansion ramp paramter1 (time at which the rate ramps)");
  params.addRequiredParam<Real>("c2", "Expansion ramp paramter2 (steepness of ramp)");
  params.addParam<Real>("press_str", std::numeric_limits<Real>::quiet_NaN(),"Cell pressure threshold/scale for stress-inhibited growth gate g_p (preferred name)");
  params.addParam<Real>("press0", std::numeric_limits<Real>::quiet_NaN(),"Alias for press_str (kept for older input files)");

  // Optional (with default values shown)
  params.addParam<Real>("k_T1_max", 0.0, "T1 logistic amplitude (0 disables)");
  params.addParam<Real>("chi_str", 0.20, "shape distortion threshold for T1 yielding.");
  params.addParam<Real>("beta_T1", 100., "T1 logistic steepness parameter");
  params.addParam<Real>("epsilon", 1e-8, "Epsilon for Miehe-style algorithmic tangent");

  params.addParam<Real>("G_gel", 1.0, "shear modulus of gel matrix");
  params.addParam<Real>("k_diss_0", 0.0, "base bond dissociation rate for gel network");
  params.addParam<Real>("phi_cell_0", 1.0, "starting cell phase fraction");

  // ---- nutrient coupling + gates (optional) ----
  params.addCoupledVar("n", "Nutrient/activator field for activity gate fa (optional; if absent fa=1)");
  params.addParam<Real>("n_c1", 1.0, "Half-activation level for fa(n)");
  params.addParam<Real>("n_c2", 2.0, "Exponent in fa(n) Hill function");
  params.addParam<Real>("gamma_n0", 0.0, "Baseline nutrient consumption rate (1/time). 0 disables consumption.");
  params.addParam<Real>("D_nutrient", 0.1, "Nutrient diffusivity used in the TL pull-back (length^2/time in chosen units).");
  params.addParam<bool>("use_crowding_diffusion", true,
                     "If true, use crowding law D_phys = D_nutrient*(1-phi_cell)/(1+phi_cell/2); otherwise D_phys = D_nutrient.");
  params.addParam<Real>("D_phys_floor", 0.0, "Floor on physical diffusivity for numerical stability.");
  params.addParam<Real>("J_floor", 1e-12, "Floor on Jacobian used in diffusion pull-back.");

  params.addParam<Real>("n_eps", 1e-8,"Small regularization for dividing by n when forming an effective MatReaction rate.");

  // ---- initial microstructure / reference cell fraction field (optional) ----
  params.addCoupledVar("phi_ref_ic", "Optional field used to initialize phi_cell_ref at t=0 (seed with RandomIC for speckled 3D). ");
  params.addParam<Real>("m_T1", 0.0,"Hill exponent for T1 gating. 0 uses the existing logistic gate.");
  return params;
}

// --------------- ctor -------------------
CellGelMixtureOpt::CellGelMixtureOpt(const InputParameters & parameters)
  : DerivativeMaterialInterface<ComputeLagrangianStressPK2>(parameters),

    // Model parameters
    _G_cell(getParam<Real>("G_cell")),
    _q_cell(getParam<Real>("q_cell")),
    _k_exp_max(getParam<Real>("k_exp_max")),
    _c1(getParam<Real>("c1")),
    _c2(getParam<Real>("c2")),
    _press_str(!std::isnan(getParam<Real>("press_str")) ? getParam<Real>("press_str") : getParam<Real>("press0")),
    _k_T1_max(getParam<Real>("k_T1_max")),
    _chi_str(getParam<Real>("chi_str")),
    _beta_T1(getParam<Real>("beta_T1")),
    _G_gel(getParam<Real>("G_gel")),
    _k_diss_0(getParam<Real>("k_diss_0")),
    _phi_cell_0(getParam<Real>("phi_cell_0")),
    _n_c1(getParam<Real>("n_c1")),
    _n_c2(getParam<Real>("n_c2")),
    _m_T1(getParam<Real>("m_T1")),
    _gamma_n0(getParam<Real>("gamma_n0")),
    _D_nutrient(getParam<Real>("D_nutrient")),
    _use_crowding_diffusion(getParam<bool>("use_crowding_diffusion")),
    _D_phys_floor(getParam<Real>("D_phys_floor")),
    _J_floor(getParam<Real>("J_floor")),
    _n_eps(getParam<Real>("n_eps")),

    // nutrient coupling
    _has_n(isCoupled("n")),
    _n(_has_n ? &coupledValue("n") : nullptr),

    // IC coupling (optional)
    _has_phi_ref_ic(isCoupled("phi_ref_ic")),
    _phi_ref_ic(_has_phi_ref_ic ? &coupledValue("phi_ref_ic") : nullptr),
    
    _epsilon(getParam<Real>("epsilon")),

    // Kinematics
    _F(getMaterialPropertyByName<RankTwoTensor>(_base_name + "deformation_gradient")),
    _F_old(getMaterialPropertyOldByName<RankTwoTensor>(_base_name + "deformation_gradient")),


    // stateful (Material) properties
    _phi_cell(declareProperty<Real>("phi_cell")),
    _phi_cell_ref(declareProperty<Real>("phi_cell_ref")),

    _phi_ref_from_ic(declareProperty<Real>("phi_ref_from_ic")),

    _phi_cell_ref_old(getMaterialPropertyOld<Real>("phi_cell_ref")),


    _sigma_cell(declareProperty<RankTwoTensor>("sigma_cell")),
    _sigma_pmat(declareProperty<RankTwoTensor>("sigma_pmat")),
    _cauchy_stress(declareProperty<RankTwoTensor>("cauchy_stress")),
    _pk1_stress(declareProperty<RankTwoTensor>("pk1_stress")),
    _dcauchy_stress_dn(declarePropertyDerivative<RankTwoTensor>("cauchy_stress", "n")),
    _dpk1_stress_dn(declarePropertyDerivative<RankTwoTensor>("pk1_stress", "n")),

    _bE_cell(declareProperty<RankTwoTensor>("bE_cell")),
    _bE_cell_old(getMaterialPropertyOld<RankTwoTensor>("bE_cell")),

    _bE_pmat(declareProperty<RankTwoTensor>("bE_pmat")),
    _bE_pmat_old(getMaterialPropertyOld<RankTwoTensor>("bE_pmat")),

    _eta(declareProperty<Real>("eta")),
    _eta_old(getMaterialPropertyOld<Real>("eta")),

    _chi(declareProperty<Real>("chi")),

    _ke(declareProperty<Real>("ke")),
    _ke_old(getMaterialPropertyOld<Real>("ke")),

    _kT1(declareProperty<Real>("kT1")),
    _kT1_old(getMaterialPropertyOld<Real>("kT1")),

    _k_diss(declareProperty<Real>("k_diss")),
    _k_diss_old(getMaterialPropertyOld<Real>("k_diss")), 

    // diagnostics / couplings
    _volume_ratio(declareProperty<Real>("volume_ratio")),
    _kh(declareProperty<Real>("kh")),
    _fa(declareProperty<Real>("fa")),
    _pressure(declareProperty<Real>("pressure")),
    _gp(declareProperty<Real>("gp")),
    _gate_tot(declareProperty<Real>("gate_tot")),
    _gamma_n_local(declareProperty<Real>("gamma_n_local")),
    _D_eff(declareProperty<RealTensorValue>("D_eff")),
    _D_phys(declareProperty<Real>("D_phys")),
    _D_ref(declareProperty<Real>("D_ref")),
    _n_source_ref(declareProperty<Real>("n_source_ref")),
    // For the nutrient sink Jacobian: d(n_source_ref)/d(n)
    // (declare this pointer in the header; see patch notes below)
    _dn_source_ref_dn(nullptr)

{
  if (std::isnan(_press_str))
    paramError("press_str", "You must provide either \'press_str\' (preferred) or its alias \'press0\'.");

  // If 'n' is coupled, provide an analytic derivative so MatBodyForce (and ADMatBodyForce)
  // can form a consistent Jacobian contribution for the nonlinear gate fa(n).
  if (_has_n)
    _dn_source_ref_dn = &declarePropertyDerivative<Real>("n_source_ref", "n");
}

// ---------- init state ----------
void
CellGelMixtureOpt::initQpStatefulProperties()
{
  const RankTwoTensor I = RankTwoTensor::Identity();

  // --- initialize reference cell fraction ---
  //------for random seeding of phi_cell ----
  //----- will need to implement other methods later -----
  Real phi_ref = _phi_cell_0;
  if (_has_phi_ref_ic)
    phi_ref = (*_phi_ref_ic)[_qp];

  phi_ref = std::min(1.0 - 1e-8, std::max(1e-8, phi_ref));

  _phi_cell_ref[_qp] = phi_ref;   // freeze reference texture
  _phi_cell[_qp]     = phi_ref;   // consistent initial current value


  // --- initialize internal state variables ---
  _bE_cell[_qp] = I;
  _bE_pmat[_qp] = I;
  _sigma_cell[_qp].zero();
  _sigma_pmat[_qp].zero();

  _eta[_qp]   = 1.0;
  _chi[_qp]   = 0.0;
  _ke[_qp]    = 0.0;
  _kT1[_qp]   = 0.0;
  _k_diss[_qp]= _k_diss_0;

  // diagnostics / couplings
  _volume_ratio[_qp]   = 1.0;
  _kh[_qp]             = 0.0;
  _fa[_qp]             = 1.0;
  _pressure[_qp]       = 0.0;
  _gp[_qp]             = 1.0;
  _gate_tot[_qp]       = 1.0;
  _gamma_n_local[_qp]  = 0.0;              
  _n_source_ref[_qp]   = 0.0;
    
  
  RealTensorValue Deff0;
    Deff0.zero();
    Deff0(0,0) = _D_nutrient;
    Deff0(1,1) = _D_nutrient;
    Deff0(2,2) = _D_nutrient;
    _D_eff[_qp] = Deff0;
  _D_phys[_qp] = _D_nutrient;
  _D_ref[_qp]  = _D_nutrient;

  /******* If we want to make some property in-homogeneous at t=0, use this style or some function like it *******/
  // const Point & X = _q_point[_qp];   // X(0), X(1), X(2) are the material coordinates
}

void
CellGelMixtureOpt::computeQpProperties()
{
  const RankTwoTensor I = RankTwoTensor::Identity();
  const Real dt = std::max(_dt, 1e-30); // time step size


  _phi_cell_ref[_qp] = _phi_cell_ref_old[_qp];

  // sample IC field at QP
  Real phi_ref_here = _phi_cell_0;
  if (_has_phi_ref_ic)
    phi_ref_here = (*_phi_ref_ic)[_qp];
  phi_ref_here = std::min(1.0 - 1e-8, std::max(1e-8, phi_ref_here));

  _phi_ref_from_ic[_qp] = phi_ref_here;  // keep this for debugging

  // Freeze the reference texture ONLY at the start (avoid clobbering on restart)
  if (_t <= std::numeric_limits<Real>::epsilon())
    _phi_cell_ref[_qp] = phi_ref_here;
  else
    _phi_cell_ref[_qp] = _phi_cell_ref_old[_qp];

  /* important notations:  old: n, current : n+1 */

  /* Find velocity gradient at time n-1 */
  RankTwoTensor ell_old = (_F[_qp] - _F_old[_qp]) * _F_old[_qp].inverse() / dt;
  const Real kdiv_old = ell_old.trace();

  /* Update quantities at the current step n */

  RankTwoTensor bE_cell_dot =
      ell_old * _bE_cell_old[_qp] + _bE_cell_old[_qp] * ell_old.transpose() - (2.0 / 3.0) * _ke_old[_qp] * _bE_cell_old[_qp] 
      - ((4.0 / 3.0) * kdiv_old + _kT1_old[_qp]) * (_bE_cell_old[_qp] - (3.0 / _bE_cell_old[_qp].inverse().trace()) * I);


  _bE_cell[_qp] = _bE_cell_old[_qp] + bE_cell_dot * dt;

  const Real JE_cell = sqrt(_bE_cell[_qp].det());

  _sigma_cell[_qp] = (_G_cell / JE_cell) * (_bE_cell[_qp] - pow(JE_cell, _q_cell) * I);

  const Real JE_pmat_old = sqrt(_bE_pmat_old[_qp].det());

  RankTwoTensor bE_pmat_dot =
      ell_old * _bE_pmat_old[_qp] + _bE_pmat_old[_qp] * ell_old.transpose() -
      _k_diss_old[_qp] * (_bE_pmat_old[_qp] - (3.0 / (_bE_pmat_old[_qp].inverse().trace() * JE_pmat_old)) * I);

  _bE_pmat[_qp] = _bE_pmat_old[_qp] + bE_pmat_dot * dt;

  // --- mixture geometry & phase fractions (using per-QP reference fraction) ---
  const Real J_def = _F[_qp].det();
  _volume_ratio[_qp] = J_def;

  const Real phi_ref = _phi_cell_ref[_qp];
  _phi_cell[_qp] = (J_def * phi_ref) / ((J_def - 1.0) * phi_ref + 1.0);

  const Real JE_pmat = sqrt(_bE_pmat[_qp].det());
  const Real k_pmat = (4.0 / 3.0) * _G_gel * (1 - _phi_cell[_qp]) / _phi_cell[_qp]; // bulk modulus of porous matrix

  _sigma_pmat[_qp] = (_G_gel / pow53(JE_pmat)) * _bE_pmat[_qp].deviatoric() + k_pmat * (JE_pmat - 1) * I;

  /* ---- update internal state variables (for next step) ---- */

  // ramped "drive" (same form as before)
  //const Real k_exp0 = _k_exp_max - _k_exp_max / (1.0 + pow(_t / _c1, _c2));

  const Real T = 2.0;
  Real s = std::min(1.0, std::max(0.0, _t / T));
  Real r = 6*std::pow(s,5) - 15*std::pow(s,4) + 10*std::pow(s,3);
  const Real k_exp0 = _k_exp_max * r;


  // --- nutrient transport kinematics (Total Lagrangian pull-back of spatial diffusion) ---
  // Q = - J F^{-1} D F^{-T} Grad_X n  ==>  Div_X Q + J s = 0
  const RankTwoTensor F_inv = _F[_qp].inverse();
  const RankTwoTensor A = F_inv * F_inv.transpose();  // = F^{-1}F^{-T}

  // ---- nutrient diffusivity (physical crowding law, then TL pull-back) ----
  const Real J_use = std::max(J_def, _J_floor);

  // Physical (spatial) diffusivity in the current configuration:
  Real D_phys = _D_nutrient;
  if (_use_crowding_diffusion)
  {
    const Real phi = std::min(std::max(_phi_cell[_qp], 0.0), 0.999999);
    D_phys = _D_nutrient * (1.0 - phi) / (1.0 + 0.5 * phi);
    D_phys = std::max(D_phys, _D_phys_floor);
  }
  _D_phys[_qp] = D_phys;

  // Scalar 1D-style diagnostic (NOT used by kernels in 2D/3D):
  _D_ref[_qp]  = D_phys / J_use;

  // Total-Lagrangian pull-back for a scalar spatial diffusivity:
  //   D_eff = J * D_phys * F^{-1} F^{-T}
  const Real pref = J_use * D_phys;

  RealTensorValue Deff;
  for (unsigned i = 0; i < 3; ++i)
    for (unsigned j = 0; j < 3; ++j)
      Deff(i, j) = pref * A(i, j);

  _D_eff[_qp] = Deff;


  // ---- nutrient/activity gate fa(n) ----
  Real fa_here = 1.0;
  Real n_pos_here = 0.0;
  Real dfa_dn_here = 0.0;
  if (_has_n)
  {
    const Real n_raw = (*_n)[_qp];
    const Real eps   = std::max(_n_eps, 1e-16);
    // smooth positive-part: n_+ = 0.5*(n + sqrt(n^2 + eps^2))
    const Real s = std::sqrt(n_raw * n_raw + eps * eps);
    n_pos_here = 0.5 * (n_raw + s);

    const Real m   = _n_c2;
    const Real c   = std::pow(_n_c1, m);
    const Real a   = std::pow(n_pos_here, m);
    const Real den = a + c;

    fa_here = (den > 0.0) ? a / den : 0.0;

    // d fa / d n (needed for consistent Jacobian in MatBodyForce/ADMatBodyForce)
    // n_+ = 0.5*(n + sqrt(n^2 + eps^2))  =>  dn_+/dn = 0.5*(1 + n/sqrt(n^2+eps^2))
    const Real dnpos_dn = 0.5 * (1.0 + n_raw / s);

    // a = n_+^m, fa = a/(a+c) => dfa/dn = (da/dn)*c/(a+c)^2
    if (n_pos_here > 0.0 && m > 0.0 && den > 0.0)
    {
      const Real da_dn = m * std::pow(n_pos_here, m - 1.0) * dnpos_dn;
      dfa_dn_here = (da_dn * c) / (den * den);
    }
  }
  _fa[_qp] = fa_here;


  // ---- pressure gate gp(p_cell): suppresses expansion under compression ----
  // Convention: press_cell < 0 in compression (since sigma trace/3 is negative).
  const Real press_cell = _sigma_cell[_qp].trace() / 3.0;
  const Real x          = press_cell / _press_str;
  const Real gp_here    = (press_cell < 0.0) ? std::exp(-std::log(2.0) * x * x) : 1.0;
     
  _gp[_qp] = gp_here;

  // ---- combined gate and final ke for *next* step ----
  const Real gate_tot = fa_here * gp_here;
  _gate_tot[_qp] = gate_tot;

  _ke[_qp] = k_exp0 * gate_tot;

  // ---- cell number ratio eta ----
  // Use gated volumetric rate; forward fd
  const Real kh_here = fa_here * kdiv_old;
  _kh[_qp] = kh_here;
  _eta[_qp] = _eta_old[_qp] * std::exp(kh_here * dt);

  // ---- nutrient consumption scaling ----
  // gamma_n_local is a positive reaction-rate coefficient (consumption sink handled by MatReaction sign)
  if (_gamma_n0 == 0.0)
    _gamma_n_local[_qp] = 0.0;
  else
  {

    // nurrent consumption scaled by cell fraction, cell density, and activity gate
    const Real phi_ratio = (_phi_cell[_qp] > 0.0 && _phi_cell_0 > 0.0) ? (_phi_cell[_qp] / _phi_cell_0) : 0.0;
    _gamma_n_local[_qp] = _gamma_n0 * _phi_cell[_qp] * fa_here;
  }

  _n_source_ref[_qp] = - J_def * _gamma_n_local[_qp]; // pulled-back reaction rate for MatReaction: - (n_source_ref) * n

  // Provide derivative for the sink term Jacobian: d(n_source_ref)/d(n)
  // This is the *key* piece that prevents SNES from stalling when the gate fa(n) turns on.
  if (_dn_source_ref_dn)
  {
    if (_gamma_n0 == 0.0)
      (*_dn_source_ref_dn)[_qp] = 0.0;
    else
      (*_dn_source_ref_dn)[_qp] = J_def * (-_gamma_n0 * _phi_cell[_qp] * dfa_dn_here);
  }

  // ---- shape metric chi ----
  RankTwoTensor dev_bE_cell = _bE_cell[_qp].deviatoric();
   _chi[_qp] = sqrt(1.5 * dev_bE_cell.doubleContraction(dev_bE_cell)) / _bE_cell[_qp].trace();

  // ---- T1 gate: keep existing logistic by default; enable Hill if m_T1>0 ----
  // ---- T1 gate: skip work if disabled; otherwise logistic by default; enable Hill if m_T1>0 ----
  if (_k_T1_max == 0.0)
    _kT1[_qp] = 0.0;
  else if (_m_T1 > 0.0)
  {
    const Real chi_here  = std::max(0.0, _chi[_qp]);
    const Real chi_crit  = _chi_str;
    const Real m         = _m_T1;

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

    _k_diss[_qp] = _k_diss_0; // For future use: can make this depend on other factors

  ComputeLagrangianStressPK2::computeQpProperties();
}

//----------Cauchy stress-----------//
//-------compute pressure ----------//
void
CellGelMixtureOpt::computeQpCauchyStress()
{
  _cauchy_stress[_qp] = _phi_cell[_qp] * _sigma_cell[_qp] + (1.0 - _phi_cell[_qp]) * _sigma_pmat[_qp];
  _pressure[_qp] = _cauchy_stress[_qp].trace() / 3.0;
}

// ---------- main ----------
void
CellGelMixtureOpt::computeQpPK2Stress()
{
  // Need to explicitly compute Cauchy stress here because the _cauchy_stress variable is not
  // guaranteed to be updated before this function is called (depends on MOOSE execution order)
  const RankTwoTensor sigma_cauchy = _phi_cell[_qp] * _sigma_cell[_qp] + (1.0 - _phi_cell[_qp]) * _sigma_pmat[_qp];

  // Keep Material outputs consistent with AuxKernels (cauchy_stress, pressure)
  _cauchy_stress[_qp] = sigma_cauchy;
  _pressure[_qp]      = sigma_cauchy.trace() / 3.0;

  // Kinematics at n+1 and n
  const RankTwoTensor & F_np1 = _F[_qp];
  const RankTwoTensor & F_n   = _F_old[_qp];
  const Real J = F_np1.det();

  // Inverses used repeatedly (base state)
  const RankTwoTensor Finv  = F_np1.inverse();
  const RankTwoTensor FinvT = Finv.transpose();
  const RankTwoTensor Fn_inv = F_n.inverse();

  // PK2 stress from Cauchy
  _S[_qp] = Finv * sigma_cauchy * FinvT * J;

  // ---- compute algorithmic tangent dS/dE (finite-difference, Miehe) ----
  const Real dt = std::max(_dt, 1e-30);
  const Real inv_dt = 1.0 / dt;

  // Base velocity gradient at old step (used to form trial ell efficiently)
  const RankTwoTensor ell_old_base = (F_np1 - F_n) * Fn_inv * inv_dt;
  (void)ell_old_base;

  // Precompute constants used in trial updates (using OLD values!!!!)
  const Real JE_pmat_old        = std::sqrt(_bE_pmat_old[_qp].det());
  const Real bE_cell_inv_trace  = _bE_cell_old[_qp].inverse().trace();
  const Real bE_pmat_inv_trace  = _bE_pmat_old[_qp].inverse().trace();
  const Real phi_ref            = _phi_cell_ref[_qp];

  for (unsigned int C = 0; C < 3; ++C)
  {
    for (unsigned int D = C; D < 3; ++D)
    {
      // Symmetric perturbations in E (6 independent E directions give 6 stress comps each - 36 total)
      // I am using dE notation instead - it is 1/2 the usual dC perturbation
      RankTwoTensor dE;
      dE.zero();
      if (C == D)
        dE(C, C) = _epsilon;
      else
      {
        dE(C, D) = 0.5 * _epsilon;
        dE(D, C) = 0.5 * _epsilon;
      }

      // from Miehe: ΔF = F^{-T} ΔE (1/2 ΔC)  (F_trial = F + ΔF)
      const RankTwoTensor dF = FinvT * dE;

      // EXACT: compute inverse and det of the perturbed deformation gradient
      const RankTwoTensor F_trial    = F_np1 + dF;
      const Real          J_trial    = F_trial.det();
      const RankTwoTensor Finv_trial = F_trial.inverse();
      const RankTwoTensor FinvT_trial = Finv_trial.transpose();


      // Velocity gradient at old step for trial state
      const RankTwoTensor ell_old_trial = (F_trial - F_n) * Fn_inv * inv_dt;
      const Real kh_trial = ell_old_trial.trace();

      // --- cell phase trial update ---
      // --- using OLD ell and OLD states to form be_dot_trial
      RankTwoTensor bE_cell_dot_trial = ell_old_trial * _bE_cell_old[_qp] + _bE_cell_old[_qp] * ell_old_trial.transpose() -
          (2.0 / 3.0) * _ke_old[_qp] * _bE_cell_old[_qp] - ((4.0 / 3.0) * kh_trial + _kT1_old[_qp]) * (_bE_cell_old[_qp] - (3.0 / bE_cell_inv_trace) * RankTwoTensor::Identity());


      //--- complete cell trial be and stress update ----
      const RankTwoTensor bE_cell_trial = _bE_cell_old[_qp] + bE_cell_dot_trial * dt;
      const Real JE_cell_trial = std::sqrt(bE_cell_trial.det());
      const RankTwoTensor sigma_cell_trial = (_G_cell / JE_cell_trial) * (bE_cell_trial - std::pow(JE_cell_trial, _q_cell) * RankTwoTensor::Identity());

      // --- porous matrix trial update ---
      RankTwoTensor bE_pmat_dot_trial = ell_old_trial * _bE_pmat_old[_qp] + _bE_pmat_old[_qp] * ell_old_trial.transpose() -
          _k_diss_old[_qp] *(_bE_pmat_old[_qp] - (3.0 / (bE_pmat_inv_trace * JE_pmat_old)) * RankTwoTensor::Identity());

      // new be and Je for matrix phase at trial step
      const RankTwoTensor bE_pmat_trial = _bE_pmat_old[_qp] + bE_pmat_dot_trial * dt;
      const Real JE_pmat_trial = std::sqrt(bE_pmat_trial.det());

      // Mixture cell fraction at trial state (uses per-QP reference fraction)
      const Real phi_cell_trial = (J_trial * phi_ref) / ((J_trial - 1.0) * phi_ref + 1.0);

      // porous matrix stress at trial step
      const Real k_pmat_trial = (4.0 / 3.0) * _G_gel * (1.0 - phi_cell_trial) / phi_cell_trial;
      const RankTwoTensor sigma_pmat_trial = (_G_gel / pow53(JE_pmat_trial)) * bE_pmat_trial.deviatoric() + k_pmat_trial * (JE_pmat_trial - 1.0) * RankTwoTensor::Identity();

      // --- trial Cauchy and PK2 stress ---
      const RankTwoTensor cauchy_stress_trial = phi_cell_trial * sigma_cell_trial + (1.0 - phi_cell_trial) * sigma_pmat_trial;
      const RankTwoTensor S_trial = Finv_trial * cauchy_stress_trial * FinvT_trial * J_trial;

      // Fill Cons. tangent, minor symmetries (36 ind) and using symmetry of S and E to only run 6 times 
      for (unsigned int A = 0; A < 3; ++A)
      {
        for (unsigned int B = A; B < 3; ++B)
        {
          const Real val = (S_trial(A, B) - _S[_qp](A, B)) / _epsilon;

          _C[_qp](A, B, C, D) = val;

          if (A != B)
            _C[_qp](B, A, C, D) = val;

          if (C != D)
          {
            _C[_qp](A, B, D, C) = val;
            if (A != B)
              _C[_qp](B, A, D, C) = val;
          }
        }
      }
    }
  }

  // Provide PK1 stress for TL divergence kernels: P = F * S
  _pk1_stress[_qp] = _F[_qp] * _S[_qp];
  _dpk1_stress_dn[_qp].zero();
}
