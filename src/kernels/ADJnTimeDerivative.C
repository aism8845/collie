#include "ADJnTimeDerivative.h"

registerMooseObject("collieApp", ADJnTimeDerivative);

InputParameters
ADJnTimeDerivative::validParams()
{
  InputParameters params = ADTimeKernelValue::validParams();
  params.addRequiredParam<MaterialPropertyName>("coef", "J material property name.");
  params.addRequiredParam<MaterialPropertyName>("coef_dot", "Jdot material property name.");
  params.addParam<bool>("include_jdot", true, "If true, include Jdot * n term.");
  params.addClassDescription("AD time derivative for J*n with optional Jdot*n.");
  return params;
}

ADJnTimeDerivative::ADJnTimeDerivative(const InputParameters & parameters)
  : ADTimeKernelValue(parameters),
    _coef(getADMaterialProperty<Real>(getParam<MaterialPropertyName>("coef"))),
    _coef_dot(getADMaterialProperty<Real>(getParam<MaterialPropertyName>("coef_dot"))),
    _include_jdot(getParam<bool>("include_jdot"))
{
}

ADReal
ADJnTimeDerivative::precomputeQpResidual()
{
  ADReal res = _coef[_qp] * _u_dot[_qp];
  if (_include_jdot)
    res += _coef_dot[_qp] * _u[_qp];
  return res;
}
