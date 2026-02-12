#include "ScalarCoupledDirichletBC.h"

#include "FEProblemBase.h"
#include "MooseVariableScalar.h"

registerMooseObject("collieApp", ScalarCoupledDirichletBC);

InputParameters
ScalarCoupledDirichletBC::validParams()
{
  InputParameters params = NodalBC::validParams();
  params.addRequiredCoupledVar("bath", "Coupled scalar bath concentration variable.");
  params.addClassDescription("Pins a nodal variable to a coupled scalar bath variable.");
  return params;
}

ScalarCoupledDirichletBC::ScalarCoupledDirichletBC(const InputParameters & parameters)
  : NodalBC(parameters), _bath(coupledScalarValue("bath")), _bath_var_number(coupledScalar("bath"))
{
  if (!isCoupledScalar("bath"))
    paramError("bath", "must couple to a scalar variable.");
}

Real
ScalarCoupledDirichletBC::computeQpResidual()
{
  return _u[_qp] - _bath[0];
}

Real
ScalarCoupledDirichletBC::computeQpJacobian()
{
  return 1.0;
}

void
ScalarCoupledDirichletBC::computeOffDiagJacobianScalar(unsigned int jvar)
{
  if (!_var.isNodalDefined() || jvar != _bath_var_number)
    return;

  MooseVariableScalar & bath_var = _sys.getScalarVariable(_tid, jvar);
  const auto & bath_dofs = bath_var.dofIndices();
  if (bath_dofs.empty())
    return;

  const dof_id_type row = _var.nodalDofIndex();
  addJacobianElement(_fe_problem.assembly(0, _sys.number()), -1.0, row, bath_dofs[0], 1.0);
}
