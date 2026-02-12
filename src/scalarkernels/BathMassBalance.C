#include "BathMassBalance.h"

#include "Function.h"

registerMooseObject("collieApp", BathMassBalance);

InputParameters
BathMassBalance::validParams()
{
  InputParameters params = ScalarKernel::validParams();
  params.addRequiredParam<PostprocessorName>(
      "flux_pp",
      "Postprocessor that provides outward diffusive nutrient flux integrated on the bath boundary.");
  params.addRequiredParam<Real>("V_bath", "Bath volume for converting flux to concentration rate.");
  params.addParam<Real>("n_feed", 1.0, "Feed concentration used by the refill term.");
  params.addParam<FunctionName>(
      "k_refill",
      "Optional refill-rate function of time. If omitted, refill is disabled.");
  params.addClassDescription("Applies scalar bath nutrient mass balance with optional refill.");
  return params;
}

BathMassBalance::BathMassBalance(const InputParameters & parameters)
  : ScalarKernel(parameters),
    _flux_pp(getPostprocessorValue("flux_pp")),
    _v_bath(getParam<Real>("V_bath")),
    _n_feed(getParam<Real>("n_feed")),
    _k_refill(isParamValid("k_refill") ? &getFunction("k_refill") : nullptr)
{
  if (_v_bath <= 0.0)
    paramError("V_bath", "must be positive.");
}

void
BathMassBalance::reinit()
{
}

Real
BathMassBalance::computeQpResidual()
{
  const Real k_refill = _k_refill ? _k_refill->value(_t) : 0.0;
  return -(_flux_pp / _v_bath + k_refill * (_n_feed - _u[_i]));
}

Real
BathMassBalance::computeQpJacobian()
{
  const Real k_refill = _k_refill ? _k_refill->value(_t) : 0.0;
  return k_refill;
}
