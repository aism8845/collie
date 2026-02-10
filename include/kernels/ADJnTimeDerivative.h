#pragma once

#include "ADTimeKernelValue.h"

class ADJnTimeDerivative : public ADTimeKernelValue
{
public:
  static InputParameters validParams();
  ADJnTimeDerivative(const InputParameters & parameters);

protected:
  virtual ADReal precomputeQpResidual() override;

  const ADMaterialProperty<Real> & _coef;
  const ADMaterialProperty<Real> & _coef_dot;
  const bool _include_jdot;
};
