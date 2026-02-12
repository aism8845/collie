#pragma once

#include "ScalarKernel.h"

class Function;

class BathMassBalance : public ScalarKernel
{
public:
  static InputParameters validParams();
  BathMassBalance(const InputParameters & parameters);
  virtual void reinit() override;

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;

  const PostprocessorValue & _flux_pp;
  const Real _v_bath;
  const Real _n_feed;
  const Function * const _k_refill;
};
