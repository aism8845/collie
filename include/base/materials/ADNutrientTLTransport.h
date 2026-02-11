#pragma once

#include "ADMaterial.h"
#include "MooseEnum.h"

class ADNutrientTLTransport : public ADMaterial
{
public:
  static InputParameters validParams();
  ADNutrientTLTransport(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;

  // coupled nutrient
  const ADVariableValue & _n;

  // coupled displacements (RZ: ux=u_r, uy=u_z)
  const ADVariableValue & _ur;
  const ADVariableValue & _uz;
  const ADVariableGradient & _grad_ur;
  const ADVariableGradient & _grad_uz;

  // old/older (Real) states for Jdot
  const VariableValue & _ur_old;
  const VariableValue & _uz_old;
  const VariableGradient & _grad_ur_old;
  const VariableGradient & _grad_uz_old;

  const VariableValue & _ur_older;
  const VariableValue & _uz_older;
  const VariableGradient & _grad_ur_older;
  const VariableGradient & _grad_uz_older;

  // axisym controls
  const bool _axisymmetric;
  const unsigned int _radial_coord;
  const Real _r_eps;
  const unsigned int _axial_coord;

  // referential cell fraction (stateful, non-AD) provided by CellGelMixtureOpt
  const MaterialProperty<Real> & _phi_cell_ref;

  // parameters
  const Real _D0;
  const Real _D_floor;
  const Real _gamma_n0;
  const Real _phi_max;
  const MooseEnum _crowding_model;
  const Real _crowd_exp;
  const Real _n_c1;
  const unsigned int _n_c2;
  const Real _n_c1_pow;
  const Real _smooth_eps_c;
  const Real _smooth_eps_D;

  // outputs used by AD kernels
  ADMaterialProperty<Real> & _J_nutr;
  ADMaterialProperty<Real> & _Jdot_nutr;
  ADMaterialProperty<RealTensorValue> & _D_eff_nutr;
  ADMaterialProperty<Real> & _n_source_ref_nutr;
  ADMaterialProperty<Real> & _jdot_rate_nutr;

  // diagnostics (unique names to avoid colliding with CellGelMixtureOpt)
  ADMaterialProperty<Real> & _D_phys_nutr;
  ADMaterialProperty<Real> & _D_ref_nutr;
  ADMaterialProperty<Real> & _D_iso_nutr;
};
