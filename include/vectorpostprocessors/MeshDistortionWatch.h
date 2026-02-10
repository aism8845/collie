#pragma once

#include "ElementVectorPostprocessor.h"

class MeshDistortionWatch : public ElementVectorPostprocessor
{
public:
  static InputParameters validParams();
  MeshDistortionWatch(const InputParameters & parameters);

  virtual void initialize() override;
  virtual void execute() override;
  virtual void finalize() override;
  virtual void threadJoin(const UserObject & y) override;

protected:
  struct Candidate
  {
    Real val;
    Real x;
    Real y;
    Real z;
    Real h;
    dof_id_type elem_id;
  };

  const MaterialProperty<RankTwoTensor> & _F;
  const Real _exclude_radius_factor;

  std::vector<Candidate> _j_candidates;
  std::vector<Candidate> _m2_candidates;

  std::vector<Real> & _j_min_1;
  std::vector<Real> & _j_x_1;
  std::vector<Real> & _j_y_1;
  std::vector<Real> & _j_min_2;
  std::vector<Real> & _j_x_2;
  std::vector<Real> & _j_y_2;

  std::vector<Real> & _m2_min_1;
  std::vector<Real> & _m2_x_1;
  std::vector<Real> & _m2_y_1;
  std::vector<Real> & _m2_min_2;
  std::vector<Real> & _m2_x_2;
  std::vector<Real> & _m2_y_2;

  void computeFirstSecond(const std::vector<Candidate> & candidates,
                          Real & min1,
                          Real & x1,
                          Real & y1,
                          Real & min2,
                          Real & x2,
                          Real & y2) const;
};
