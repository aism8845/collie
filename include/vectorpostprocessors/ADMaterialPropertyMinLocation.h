#pragma once

#include "ElementVectorPostprocessor.h"

class ADMaterialPropertyMinLocation : public ElementVectorPostprocessor
{
public:
  static InputParameters validParams();
  ADMaterialPropertyMinLocation(const InputParameters & parameters);

  virtual void initialize() override;
  virtual void execute() override;
  virtual void finalize() override;
  virtual void threadJoin(const UserObject & uo) override;

protected:
  const MaterialProperty<Real> & _prop;

  Real _local_min;
  dof_id_type _local_elem_id;
  Real _local_x;
  Real _local_y;
  Real _local_z;

  std::vector<Real> & _min_value;
  std::vector<Real> & _elem_id;
  std::vector<Real> & _x;
  std::vector<Real> & _y;
  std::vector<Real> & _z;
};
