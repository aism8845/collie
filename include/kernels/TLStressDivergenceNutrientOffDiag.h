#pragma once

#include "TotalLagrangianStressDivergenceAxisymmetricCylindrical.h"

class TLStressDivergenceNutrientOffDiag
  : public TotalLagrangianStressDivergenceAxisymmetricCylindrical
{
public:
  static InputParameters validParams();
  TLStressDivergenceNutrientOffDiag(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override { return 0.0; }
  virtual Real computeQpJacobian() override { return 0.0; }
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

  const MooseVariable * _n_var;
  const unsigned int _n_var_num;
  const MaterialProperty<RankTwoTensor> & _dpk1_dn;
};
