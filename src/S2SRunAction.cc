// -*- C++ -*-

#include "S2SRunAction.hh"

#include <G4Run.hh>
#include <G4RunManager.hh>
#include <G4Timer.hh>
#include <Randomize.hh>

#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <unistd.h>

#include "K18RunControl.hh"
#include "S2SAnaManager.hh"
#include "FuncName.hh"

namespace
{
auto& anaMan = S2SAnaManager::GetInstance();
G4Timer timer;

G4long
MakeRunSeed()
{
  const auto now =
    std::chrono::high_resolution_clock::now().time_since_epoch().count();
  const auto pid = static_cast<std::uint64_t>(getpid());
  std::uint64_t seed = static_cast<std::uint64_t>(now);
  seed ^= pid+0x9e3779b97f4a7c15ULL+(seed << 6)+(seed >> 2);
  seed &= 0x7fffffffUL;
  return static_cast<G4long>(seed ? seed : 1UL);
}
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
  // A fixed seed makes regression runs reproducible.  Production runs use a
  // high-resolution clock combined with the process id.
  G4long seed = 0;
  if(K18RunControl::ReadEnvironmentLong("K18G4_SEED", seed)){
    if(seed <= 0){
      G4cerr << "K18G4_SEED must be positive; got " << seed << G4endl;
      std::exit(EXIT_FAILURE);
    }
    G4Random::setTheSeed(seed);
    G4cout << "   Fixed seed (K18G4_SEED) = " << G4Random::getTheSeed()
           << G4endl;
  } else {
    G4Random::setTheSeed(MakeRunSeed());
  }
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
