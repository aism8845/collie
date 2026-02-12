#pragma once

#include "NodalBC.h"

class ScalarCoupledDirichletBC : public NodalBC
{
public:
  static InputParameters validParams();
  ScalarCoupledDirichletBC(const InputParameters & parameters);

  virtual void computeOffDiagJacobianScalar(unsigned int jvar) override;

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;

  const VariableValue & _bath;
  const unsigned int _bath_var_number;
};
