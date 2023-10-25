// -*- C++ -*-

#include "S2SRunAction.hh"

#include <G4Run.hh>
#include <G4RunManager.hh>
#include <G4Timer.hh>
#include <Randomize.hh>

#include "S2SAnaManager.hh"
#include "FuncName.hh"

namespace
{
auto& anaMan = S2SAnaManager::GetInstance();
G4Timer timer;
}

//_____________________________________________________________________________
S2SRunAction::S2SRunAction()
  : G4UserRunAction()
{
}

//_____________________________________________________________________________
S2SRunAction::~S2SRunAction()
{
}

//_____________________________________________________________________________
void
S2SRunAction::BeginOfRunAction(const G4Run* aRun)
{
  G4cout << FUNC_NAME << G4endl
	 << "   Run# = " << aRun->GetRunID() << G4endl;
  G4Random::setTheSeed(std::time(nullptr));
  // G4int initSeed = GetIntFromKernelEntropyPool()&0x7FFFFFFF;
  // G4Random::setTheSeed(initSeed);

#ifdef DEBUG
  G4cout << "   Initial Seed = " << G4Random::getTheSeed() << G4endl;
  G4Random::showEngineStatus();
#endif

  anaMan.BeginOfRun( aRun );
  timer.Start();
}

//_____________________________________________________________________________
void
S2SRunAction::EndOfRunAction(const G4Run* aRun)
{
  timer.Stop();
  anaMan.EndOfRun(aRun);
  G4cout << FUNC_NAME << G4endl
	 << "   Process end  = " << timer.GetClockTime()
	 << "   Event number = " << aRun->GetNumberOfEvent() << G4endl
	 << "   Elapsed time = " << timer << G4endl << G4endl;
}
