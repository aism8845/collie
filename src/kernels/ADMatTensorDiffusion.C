#include "ADMatTensorDiffusion.h"

registerMooseObject("collieApp", ADMatTensorDiffusion);

InputParameters
ADMatTensorDiffusion::validParams()
{
  InputParameters params = ADKernelGrad::validParams();
  params.addRequiredParam<MaterialPropertyName>("diffusivity", "Tensor diffusivity material property name.");
  params.addClassDescription("AD tensor diffusion using a material property.");
  return params;
}

ADMatTensorDiffusion::ADMatTensorDiffusion(const InputParameters & parameters)
  : ADKernelGrad(parameters),
    _diffusivity(getADMaterialProperty<RealTensorValue>(
        getParam<MaterialPropertyName>("diffusivity")))
{
}

ADRealVectorValue
ADMatTensorDiffusion::precomputeQpResidual()
{
  return _diffusivity[_qp] * _grad_u[_qp];
}
