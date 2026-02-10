#include "MeshDistortionWatch.h"

#include "MooseMesh.h"

#include <cmath>
#include <limits>

registerMooseObject("collieApp", MeshDistortionWatch);

InputParameters
MeshDistortionWatch::validParams()
{
  InputParameters params = ElementVectorPostprocessor::validParams();
  params.addParam<MaterialPropertyName>("deformation_gradient_property",
                                         "deformation_gradient",
                                         "Mechanics deformation gradient material property name.");
  params.addParam<Real>("exclude_radius_factor",
                        2.0,
                        "Second-worst exclusion radius = factor * h(first-worst element).");
  params.set<ExecFlagEnum>("execute_on") = "timestep_end";
  params.addClassDescription(
      "Tracks two worst locations per step for J=det(F) and a secondary distortion metric "
      "metric2=F_theta_theta in RZ.");
  return params;
}

MeshDistortionWatch::MeshDistortionWatch(const InputParameters & parameters)
  : ElementVectorPostprocessor(parameters),
    _F(getMaterialProperty<RankTwoTensor>(
        getParam<MaterialPropertyName>("deformation_gradient_property"))),
    _exclude_radius_factor(getParam<Real>("exclude_radius_factor")),
    _j_min_1(declareVector("J_min_1")),
    _j_x_1(declareVector("J_min_1_x")),
    _j_y_1(declareVector("J_min_1_y")),
    _j_min_2(declareVector("J_min_2")),
    _j_x_2(declareVector("J_min_2_x")),
    _j_y_2(declareVector("J_min_2_y")),
    _m2_min_1(declareVector("metric2_min_1")),
    _m2_x_1(declareVector("metric2_min_1_x")),
    _m2_y_1(declareVector("metric2_min_1_y")),
    _m2_min_2(declareVector("metric2_min_2")),
    _m2_x_2(declareVector("metric2_min_2_x")),
    _m2_y_2(declareVector("metric2_min_2_y"))
{
}

void
MeshDistortionWatch::initialize()
{
  _j_candidates.clear();
  _m2_candidates.clear();

  _j_min_1.clear();
  _j_x_1.clear();
  _j_y_1.clear();
  _j_min_2.clear();
  _j_x_2.clear();
  _j_y_2.clear();

  _m2_min_1.clear();
  _m2_x_1.clear();
  _m2_y_1.clear();
  _m2_min_2.clear();
  _m2_x_2.clear();
  _m2_y_2.clear();
}

void
MeshDistortionWatch::execute()
{
  const unsigned int nqp = _qrule->n_points();
  if (!nqp)
    return;

  Real j_elem_min = std::numeric_limits<Real>::max();
  Real m2_elem_min = std::numeric_limits<Real>::max();

  for (unsigned int qp = 0; qp < nqp; ++qp)
  {
    const Real j = _F[qp].det();
    const Real m2 = _F[qp](2, 2); // RZ hoop stretch proxy

    if (j < j_elem_min)
      j_elem_min = j;
    if (m2 < m2_elem_min)
      m2_elem_min = m2;
  }

  const Point c = _current_elem->vertex_average();
  const Real h = _current_elem->hmax();

  _j_candidates.push_back({j_elem_min, c(0), c(1), c(2), h, _current_elem->id()});
  _m2_candidates.push_back({m2_elem_min, c(0), c(1), c(2), h, _current_elem->id()});
}

void
MeshDistortionWatch::threadJoin(const UserObject & y)
{
  const auto & other = static_cast<const MeshDistortionWatch &>(y);
  _j_candidates.insert(_j_candidates.end(), other._j_candidates.begin(), other._j_candidates.end());
  _m2_candidates.insert(_m2_candidates.end(), other._m2_candidates.begin(), other._m2_candidates.end());
}

void
MeshDistortionWatch::computeFirstSecond(const std::vector<Candidate> & candidates,
                                        Real & min1,
                                        Real & x1,
                                        Real & y1,
                                        Real & min2,
                                        Real & x2,
                                        Real & y2) const
{
  const Real inf = std::numeric_limits<Real>::max();
  const Real nan = std::numeric_limits<Real>::quiet_NaN();

  Candidate local_first{inf, 0.0, 0.0, 0.0, 0.0, DofObject::invalid_id};
  for (const auto & c : candidates)
    if (c.val < local_first.val)
      local_first = c;

  Real min1_copy = local_first.val;
  unsigned int rank1 = 0;
  comm().minloc(min1_copy, rank1);

  Real x1_out = 0.0;
  Real y1_out = 0.0;
  Real h1_out = 0.0;
  if (rank1 == processor_id() && local_first.elem_id != DofObject::invalid_id)
  {
    x1_out = local_first.x;
    y1_out = local_first.y;
    h1_out = local_first.h;
  }

  comm().sum(x1_out);
  comm().sum(y1_out);
  comm().sum(h1_out);

  min1 = min1_copy;
  x1 = x1_out;
  y1 = y1_out;

  const Real r_excl = _exclude_radius_factor * h1_out;
  const Real r2_excl = r_excl * r_excl;

  Candidate local_second{inf, 0.0, 0.0, 0.0, 0.0, DofObject::invalid_id};
  for (const auto & c : candidates)
  {
    if (c.elem_id == local_first.elem_id)
      continue;

    const Real dx = c.x - x1_out;
    const Real dy = c.y - y1_out;
    const Real d2 = dx * dx + dy * dy;
    if (d2 <= r2_excl)
      continue;

    if (c.val < local_second.val)
      local_second = c;
  }

  Real min2_copy = local_second.val;
  unsigned int rank2 = 0;
  comm().minloc(min2_copy, rank2);

  Real x2_out = 0.0;
  Real y2_out = 0.0;
  if (rank2 == processor_id() && local_second.elem_id != DofObject::invalid_id)
  {
    x2_out = local_second.x;
    y2_out = local_second.y;
  }

  comm().sum(x2_out);
  comm().sum(y2_out);

  if (min2_copy >= 0.5 * inf)
  {
    min2 = nan;
    x2 = nan;
    y2 = nan;
  }
  else
  {
    min2 = min2_copy;
    x2 = x2_out;
    y2 = y2_out;
  }
}

void
MeshDistortionWatch::finalize()
{
  Real j1 = 0.0, jx1 = 0.0, jy1 = 0.0, j2 = 0.0, jx2 = 0.0, jy2 = 0.0;
  computeFirstSecond(_j_candidates, j1, jx1, jy1, j2, jx2, jy2);

  Real m21 = 0.0, m2x1 = 0.0, m2y1 = 0.0, m22 = 0.0, m2x2 = 0.0, m2y2 = 0.0;
  computeFirstSecond(_m2_candidates, m21, m2x1, m2y1, m22, m2x2, m2y2);

  _j_min_1.push_back(j1);
  _j_x_1.push_back(jx1);
  _j_y_1.push_back(jy1);
  _j_min_2.push_back(j2);
  _j_x_2.push_back(jx2);
  _j_y_2.push_back(jy2);

  _m2_min_1.push_back(m21);
  _m2_x_1.push_back(m2x1);
  _m2_y_1.push_back(m2y1);
  _m2_min_2.push_back(m22);
  _m2_x_2.push_back(m2x2);
  _m2_y_2.push_back(m2y2);
}
