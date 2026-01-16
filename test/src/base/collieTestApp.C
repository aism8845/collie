//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
#include "collieTestApp.h"
#include "collieApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "MooseSyntax.h"

InputParameters
collieTestApp::validParams()
{
  InputParameters params = collieApp::validParams();
  params.set<bool>("use_legacy_material_output") = false;
  params.set<bool>("use_legacy_initial_residual_evaluation_behavior") = false;
  return params;
}

collieTestApp::collieTestApp(const InputParameters & parameters) : MooseApp(parameters)
{
  collieTestApp::registerAll(
      _factory, _action_factory, _syntax, getParam<bool>("allow_test_objects"));
}

collieTestApp::~collieTestApp() {}

void
collieTestApp::registerAll(Factory & f, ActionFactory & af, Syntax & s, bool use_test_objs)
{
  collieApp::registerAll(f, af, s);
  if (use_test_objs)
  {
    Registry::registerObjectsTo(f, {"collieTestApp"});
    Registry::registerActionsTo(af, {"collieTestApp"});
  }
}

void
collieTestApp::registerApps()
{
  registerApp(collieApp);
  registerApp(collieTestApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
// External entry point for dynamic application loading
extern "C" void
collieTestApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  collieTestApp::registerAll(f, af, s);
}
extern "C" void
collieTestApp__registerApps()
{
  collieTestApp::registerApps();
}
