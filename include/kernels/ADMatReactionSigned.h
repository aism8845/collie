#pragma once

#include "ADKernelValue.h"

class ADMatReactionSigned : public ADKernelValue
{
public:
  static InputParameters validParams();
  ADMatReactionSigned(const InputParameters & parameters);

protected:
  virtual ADReal precomputeQpResidual() override;

  const ADMaterialProperty<Real> & _reaction_rate;
};
