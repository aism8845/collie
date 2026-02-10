#include "ADMatReactionSigned.h"

registerMooseObject("collieApp", ADMatReactionSigned);

InputParameters
ADMatReactionSigned::validParams()
{
  InputParameters params = ADKernelValue::validParams();
  params.addRequiredParam<MaterialPropertyName>("reaction_rate", "Reaction rate material property name.");
  params.addClassDescription("AD reaction kernel with MatReaction sign convention.");
  return params;
}

ADMatReactionSigned::ADMatReactionSigned(const InputParameters & parameters)
  : ADKernelValue(parameters),
    _reaction_rate(getADMaterialProperty<Real>(
        getParam<MaterialPropertyName>("reaction_rate")))
{
}

ADReal
ADMatReactionSigned::precomputeQpResidual()
{
  return -_reaction_rate[_qp] * _u[_qp];
}
