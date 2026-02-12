#include "SideTensorDiffusiveFluxIntegral.h"

registerMooseObject("collieApp", SideTensorDiffusiveFluxIntegral);

InputParameters
SideTensorDiffusiveFluxIntegral::validParams()
{
  InputParameters params = SideIntegralPostprocessor::validParams();
  params.addRequiredCoupledVar("variable", "Nutrient variable used to compute diffusive flux.");
  params.addParam<MaterialPropertyName>(
      "diffusivity",
      "D_eff_nutr",
      "Tensor diffusivity material property name used to compute boundary flux.");
  params.addClassDescription("Computes boundary-integrated tensor diffusive nutrient flux.");
  return params;
}

SideTensorDiffusiveFluxIntegral::SideTensorDiffusiveFluxIntegral(const InputParameters & parameters)
  : SideIntegralPostprocessor(parameters),
    _grad_n(coupledGradient("variable")),
    _diffusivity(getMaterialProperty<RealTensorValue>(getParam<MaterialPropertyName>("diffusivity")))
{
}

Real
SideTensorDiffusiveFluxIntegral::computeQpIntegral()
{
  const RealVectorValue flux = -(_diffusivity[_qp] * _grad_n[_qp]);
  return flux * _normals[_qp];
}
