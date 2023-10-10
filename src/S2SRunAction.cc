/*
  "S2SRunAction.cc"

  T.Gogami, 10 Nov 2014
*/

#include "S2SRunAction.hh"
//#include "S2SParamManager.hh"
#include "S2SAnaManager.hh"

#include "Randomize.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"

#include <TFile.h>
#include <TTree.h>
#include <TObjArray.h>
#include <TH1F.h>
#include <stdlib.h>

namespace
{
auto& anaMan = S2SAnaManager::GetInstance();

}

S2SRunAction::S2SRunAction()
{
}

S2SRunAction::~S2SRunAction()
{
}

void S2SRunAction::BeginOfRunAction(const G4Run* aRun)
{
  G4Random::setTheSeed(std::time(nullptr));
  // int initSeed=GetIntFromKernelEntropyPool()&0x7FFFFFFF;
  // G4Random::setTheSeed(initSeed);
#ifdef DEBUG
  int startSeed=CLHEP::HepRandom::getTheSeed();
  G4cout << "*** Initial Seed = " << startSeed << G4endl;
  G4Random::showEngineStatus();
#endif

  // G4Random::setTheEngine(new RanecuEngine);

  anaMan.BeginOfRun( aRun );
}

void S2SRunAction::EndOfRunAction(const G4Run* aRun)
{
  anaMan.EndOfRun( aRun );
  G4cout << "  END... bye ! " << G4endl;
}
