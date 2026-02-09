#include "PK1StressDivergenceNOffDiag.h"

registerMooseObject("collieApp", PK1StressDivergenceNOffDiag);

InputParameters
PK1StressDivergenceNOffDiag::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addClassDescription(
      "Adds only the off-diagonal Jacobian contribution ∂Ru/∂n using a material "
      "property dpk1_dn = dP/dn (First Piola-Kirchhoff stress).");

  params.addRequiredCoupledVar("n", "Coupled nutrient (scalar) variable.");
  params.addRequiredParam<unsigned int>("component",
                                        "Displacement component index i for this kernel (0=x, 1=y, 2=z).");
  params.addParam<MaterialPropertyName>(
      "dpk1_dn", "dpk1_dn", "Material property name for the tensor dP/dn (First Piola-Kirchhoff).");

  return params;
}

PK1StressDivergenceNOffDiag::PK1StressDivergenceNOffDiag(const InputParameters & parameters)
  : Kernel(parameters),
    _component(getParam<unsigned int>("component")),
    _n_var(coupled("n")),
    _dpk1_dn(getMaterialProperty<RankTwoTensor>(getParam<MaterialPropertyName>("dpk1_dn")))
{
}

Real
PK1StressDivergenceNOffDiag::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (jvar != _n_var)
    return 0.0;

  const RealVectorValue dP_row = _dpk1_dn[_qp].row(_component);
  return (dP_row * _grad_test[_i][_qp]) * _phi[_j][_qp];
}
