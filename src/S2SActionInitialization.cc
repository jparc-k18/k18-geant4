// -*- C++ -*-

#include "S2SActionInitialization.hh"
#include "S2SPrimaryGeneratorAction.hh"
#include "S2SRunAction.hh"
#include "S2SEventAction.hh"
#include "S2SSteppingAction.hh"

//_____________________________________________________________________________
S2SActionInitialization::S2SActionInitialization()
{
}

//_____________________________________________________________________________
S2SActionInitialization::~S2SActionInitialization()
{}

//_____________________________________________________________________________
void
S2SActionInitialization::BuildForMaster() const
{
  SetUserAction(new S2SRunAction);
}

//_____________________________________________________________________________
void
S2SActionInitialization::Build() const
{
  SetUserAction(new S2SPrimaryGeneratorAction);
  SetUserAction(new S2SRunAction);
  SetUserAction(new S2SEventAction);
  SetUserAction(new S2SSteppingAction);
}
