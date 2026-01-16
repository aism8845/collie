#include "collieApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"
#include "MooseSyntax.h"

InputParameters
collieApp::validParams()
{
  InputParameters params = MooseApp::validParams();
  params.set<bool>("use_legacy_material_output") = false;
  params.set<bool>("use_legacy_initial_residual_evaluation_behavior") = false;
  return params;
}

collieApp::collieApp(const InputParameters & parameters) : MooseApp(parameters)
{
  collieApp::registerAll(_factory, _action_factory, _syntax);
}

collieApp::~collieApp() {}

void
collieApp::registerAll(Factory & f, ActionFactory & af, Syntax & syntax)
{
  ModulesApp::registerAllObjects<collieApp>(f, af, syntax);
  Registry::registerObjectsTo(f, {"collieApp"});
  Registry::registerActionsTo(af, {"collieApp"});

  /* register custom execute flags, action syntax, etc. here */
}

void
collieApp::registerApps()
{
  registerApp(collieApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
extern "C" void
collieApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  collieApp::registerAll(f, af, s);
}
extern "C" void
collieApp__registerApps()
{
  collieApp::registerApps();
}
