#pragma once
#include "ADTimeKernelValue.h"

class ADMatCoefTimeDerivative : public ADTimeKernelValue
{
public:
  static InputParameters validParams();
  ADMatCoefTimeDerivative(const InputParameters & parameters);

protected:
  ADReal precomputeQpResidual() override;

  const ADMaterialProperty<Real> & _coef;   // AD property
};
