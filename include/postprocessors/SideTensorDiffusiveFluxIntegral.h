#pragma once

#include "SideIntegralPostprocessor.h"

class SideTensorDiffusiveFluxIntegral : public SideIntegralPostprocessor
{
public:
  static InputParameters validParams();
  SideTensorDiffusiveFluxIntegral(const InputParameters & parameters);

protected:
  virtual Real computeQpIntegral() override;

  const VariableGradient & _grad_n;
  const MaterialProperty<RealTensorValue> & _diffusivity;
};
