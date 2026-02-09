#include "ADMatCoefTimeDerivative.h"

registerMooseObject("collieApp", ADMatCoefTimeDerivative);

InputParameters
ADMatCoefTimeDerivative::validParams()
{
  InputParameters params = ADTimeKernelValue::validParams();
  params.addRequiredParam<MaterialPropertyName>(
      "coef", "AD material property multiplying u_dot (i.e. coef * du/dt)");
  return params;
}

ADMatCoefTimeDerivative::ADMatCoefTimeDerivative(const InputParameters & parameters)
  : ADTimeKernelValue(parameters),
    _coef(getADMaterialProperty<Real>(getParam<MaterialPropertyName>("coef")))
{
}

ADReal
ADMatCoefTimeDerivative::precomputeQpResidual()
{
  return _coef[_qp] * _u_dot[_qp];
}
