#pragma once

#include "Kernel.h"

/**
 * PK1StressDivergenceNOffDiag
 *
 * Add ONLY the off-diagonal Jacobian entry dRu/dn for a displacement component,
 * using a material-provided tensor dpk1_dn = dP/dn.
 *
 * Residual and diagonal Jacobian contributions are zero to avoid double-counting.
 */
class PK1StressDivergenceNOffDiag : public Kernel
{
public:
  static InputParameters validParams();

  PK1StressDivergenceNOffDiag(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override { return 0.0; }
  virtual Real computeQpJacobian() override { return 0.0; }
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

  const unsigned int _component;
  const unsigned int _n_var;

  const MaterialProperty<RankTwoTensor> & _dpk1_dn;
};
