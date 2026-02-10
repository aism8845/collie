#pragma once

#include "ComputeLagrangianStressPK2.h"
#include "DerivativeMaterialInterface.h"

// Forward Declarations
class CellGelMixtureOpt;

template <>
InputParameters validParams<CellGelMixtureOpt>();

class CellGelMixtureOpt : public DerivativeMaterialInterface<ComputeLagrangianStressPK2>
{
public:
  static InputParameters validParams();

  CellGelMixtureOpt(const InputParameters & parameters);

  void initQpStatefulProperties() override;

protected:
  virtual void computeQpProperties() override;
  virtual void computeQpPK2Stress() override;
  virtual void computeQpCauchyStress() override;

  // Parameters
  const Real _phi_cell_0;
  const Real _G_cell;
  const Real _G_gel;
  const Real _D_nutrient;

  const Real _k_exp_max;
  const Real _c1;
  const Real _c2;
  const Real _k_T1_max;
  const Real _chi_str;
  const Real _press_str;
  const Real _beta_T1;
  const Real _m_T1;
  const Real _q_cell;
  const Real _k_diss_0;
  const Real _n_c1;
  const Real _n_c2;
  const Real _gamma_n0;
  const bool _use_crowding_diffusion;
  const Real _D_phys_floor;
  const Real _J_floor;
  const Real _n_eps;

  // Flags and couplings
  const bool _has_n;
  const VariableValue * _n;
  const bool _has_phi_ref_ic;
  const VariableValue * _phi_ref_ic;

  // Optional derivative storage for nutrient source sensitivity
  MaterialProperty<Real> * _dn_source_ref_dn;

  // Epsilon for algorithmic tangent (existing)
  const Real _epsilon;

  // Kinematics
  const MaterialProperty<RankTwoTensor> & _F;
  const MaterialProperty<RankTwoTensor> & _F_old;

  // Stateful properties
  MaterialProperty<Real> & _phi_cell_ref;
  const MaterialProperty<Real> & _phi_cell_ref_old;

  MaterialProperty<RankTwoTensor> & _bE_cell;
  const MaterialProperty<RankTwoTensor> & _bE_cell_old;

  MaterialProperty<RankTwoTensor> & _bE_pmat;
  const MaterialProperty<RankTwoTensor> & _bE_pmat_old;

  const MaterialProperty<Real> & _ke_old;
  const MaterialProperty<Real> & _kT1_old;
  const MaterialProperty<Real> & _k_diss_old;
  const MaterialProperty<Real> & _eta_old;

  // Output / auxiliary properties
  MaterialProperty<Real> & _D_phys;
  MaterialProperty<Real> & _D_ref;
  MaterialProperty<Real> & _phi_cell;
  MaterialProperty<Real> & _phi_ref_from_ic;
  MaterialProperty<RankTwoTensor> & _sigma_cell;
  MaterialProperty<RankTwoTensor> & _sigma_pmat;
  MaterialProperty<RankTwoTensor> & _cauchy_stress;
  MaterialProperty<Real> & _eta;
  MaterialProperty<Real> & _chi;
  MaterialProperty<Real> & _ke;
  MaterialProperty<Real> & _kT1;
  MaterialProperty<Real> & _k_diss;
  MaterialProperty<Real> & _volume_ratio;
  MaterialProperty<Real> & _kh;
  MaterialProperty<Real> & _fa;
  MaterialProperty<Real> & _pressure;
  MaterialProperty<Real> & _gp;
  MaterialProperty<Real> & _gate_tot;
  MaterialProperty<Real> & _gamma_n_local;
  MaterialProperty<RealTensorValue> & _D_eff;
  MaterialProperty<Real> & _n_source_ref;
  MaterialProperty<RankTwoTensor> & _dcauchy_stress_dn;

  // Mechanics–nutrient coupling: PK1 stress and its derivative w.r.t nutrient
  MaterialProperty<RankTwoTensor> & _pk1_stress;
  MaterialProperty<RankTwoTensor> & _dpk1_stress_dn;
};
