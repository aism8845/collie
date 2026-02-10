#include "ADMaterialPropertyMinLocation.h"

#include "MooseMesh.h"
#include "MooseUtils.h"

#include <limits>

registerMooseObject("collieApp", ADMaterialPropertyMinLocation);

InputParameters
ADMaterialPropertyMinLocation::validParams()
{
  InputParameters params = ElementVectorPostprocessor::validParams();
  params.addRequiredParam<MaterialPropertyName>("ad_material_property",
                                                "AD material property to scan for a minimum.");
  params.set<ExecFlagEnum>("execute_on") = "nonlinear timestep_end";
  params.addClassDescription("Find minimum AD material property value and its location.");
  return params;
}

ADMaterialPropertyMinLocation::ADMaterialPropertyMinLocation(const InputParameters & parameters)
  : ElementVectorPostprocessor(parameters),
    _prop(getMaterialProperty<Real>(getParam<MaterialPropertyName>("ad_material_property"))),
    _local_min(std::numeric_limits<Real>::max()),
    _local_elem_id(std::numeric_limits<dof_id_type>::max()),
    _local_x(0.0),
    _local_y(0.0),
    _local_z(0.0),
    _min_value(declareVector("min_value")),
    _elem_id(declareVector("elem_id")),
    _x(declareVector("x")),
    _y(declareVector("y")),
    _z(declareVector("z"))
{
}

void
ADMaterialPropertyMinLocation::initialize()
{
  _local_min = std::numeric_limits<Real>::max();
  _local_elem_id = std::numeric_limits<dof_id_type>::max();
  _local_x = 0.0;
  _local_y = 0.0;
  _local_z = 0.0;

  _min_value.clear();
  _elem_id.clear();
  _x.clear();
  _y.clear();
  _z.clear();
}

void
ADMaterialPropertyMinLocation::execute()
{
  const unsigned int nqp = _qrule->n_points();
  for (unsigned int qp = 0; qp < nqp; ++qp)
  {
    const Real val = MetaPhysicL::raw_value(_prop[qp]);
    if (val < _local_min)
    {
      _local_min = val;
      _local_elem_id = _current_elem->id();
      _local_x = _q_point[qp](0);
      _local_y = _q_point[qp](1);
      _local_z = _q_point[qp](2);
    }
  }
}

void
ADMaterialPropertyMinLocation::threadJoin(const UserObject & uo)
{
  const auto & other = static_cast<const ADMaterialPropertyMinLocation &>(uo);
  if (other._local_min < _local_min)
  {
    _local_min = other._local_min;
    _local_elem_id = other._local_elem_id;
    _local_x = other._local_x;
    _local_y = other._local_y;
    _local_z = other._local_z;
  }
}

void
ADMaterialPropertyMinLocation::finalize()
{
  Real min_copy = _local_min;
  unsigned int rank = 0;
  comm().minloc(min_copy, rank);

  Real elem_id_out = 0.0;
  Real x_out = 0.0;
  Real y_out = 0.0;
  Real z_out = 0.0;

  if (rank == processor_id())
  {
    elem_id_out = static_cast<Real>(_local_elem_id);
    x_out = _local_x;
    y_out = _local_y;
    z_out = _local_z;
  }

  comm().sum(elem_id_out);
  comm().sum(x_out);
  comm().sum(y_out);
  comm().sum(z_out);

  _min_value.push_back(min_copy);
  _elem_id.push_back(elem_id_out);
  _x.push_back(x_out);
  _y.push_back(y_out);
  _z.push_back(z_out);
}
