#pragma once

#include "ADKernelGrad.h"

class ADMatTensorDiffusion : public ADKernelGrad
{
public:
  static InputParameters validParams();
  ADMatTensorDiffusion(const InputParameters & parameters);

protected:
  virtual ADRealVectorValue precomputeQpResidual() override;

  const ADMaterialProperty<RealTensorValue> & _diffusivity;
};
