#pragma once

#include "ADMaterial.h"
#include "RankTwoTensor.h"
#include "ADRankTwoTensorForward.h"
#include "MathUtils.h"

/**
 * TL (reference-form) nutrient transport coefficients for scalar n(X,t).
 *
 * This material is designed to produce FULL AD coupling of the nutrient residual wrt
 * displacements by building F = I + Grad(u) directly from AD displacement gradients.
 *
 * Outputs (AD material properties):
 *   - J_nutr              : det(F)
 *   - Jdot_nutr           : time derivative of J (BDF1/BDF2-consistent)
 *   - D_eff_nutr          : pulled-back diffusivity tensor K = J * D_phys(phi) * F^{-1} * F^{-T}
 *   - n_source_ref_nutr   : "rate" for ADMatReaction to add + J * gamma_local * n
 *                           (since ADMatReaction adds -rate*u, we set rate = -J*gamma_local)
 *   - jdot_rate_nutr      : "rate" for ADMatReaction to add + Jdot * n (set rate = -Jdot)
 *
 * Notes:
 * - Avoids requesting the framework 'deformation_gradient' property (prevents AD/non-AD collision).
 * - Requires coupling to disp_r, disp_z (ux, uy in RZ).
 */
class ADNutrientTLTransport : public ADMaterial
{
public:
  static InputParameters validParams();
  ADNutrientTLTransport(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;

  // Coupled variables
  const ADVariableValue & _n;

  const ADVariableValue & _ur;
  const ADVariableValue & _uz;

  const ADVariableGradient & _grad_ur;
  const ADVariableGradient & _grad_uz;

  // Old/older for Jdot
  const VariableValue & _ur_old;
  const VariableValue & _uz_old;
  const VariableGradient & _grad_ur_old;
  const VariableGradient & _grad_uz_old;

  const VariableValue & _ur_older;
  const VariableValue & _uz_older;
  const VariableGradient & _grad_ur_older;
  const VariableGradient & _grad_uz_older;

  // Axisymmetry handling (RZ)
  const bool _axisymmetric;
  const unsigned int _radial_coord;
  const Real _r_eps;

  // Convenience: in 2D, the "other" coordinate index (axial if radial_coord is radial)
  const unsigned int _axial_coord;

  // Referenital cell volume fraction (stateful, from your mechanics material)
  const MaterialProperty<Real> & _phi_cell_ref;

  // Parameters
  const Real _D0;
  const Real _D_floor;
  const Real _gamma_n0;
  const Real _phi_max;
  const Real _n_c1;
  const unsigned int _n_c2;     // integer Hill exponent
  const Real _n_c1_pow;         // n_c1^n_c2 (precomputed)

  // Output AD properties
  ADMaterialProperty<Real> & _J_nutr;
  ADMaterialProperty<Real> & _Jdot_nutr;
  ADMaterialProperty<RealTensorValue> & _D_eff_nutr;
  ADMaterialProperty<Real> & _n_source_ref_nutr;
  ADMaterialProperty<Real> & _jdot_rate_nutr;
};
