// -*- C++ -*-

#include "S2SPhysicsList.hh"

#include <iomanip>
#include <CLHEP/Units/SystemOfUnits.h>

#include <globals.hh>
#include <G4ios.hh>

#include <G4DecayPhysics.hh>
#include <G4EmStandardPhysics.hh>
#include <G4EmExtraPhysics.hh>
#include <G4IonPhysics.hh>
#include <G4StoppingPhysics.hh>
#include <G4HadronElasticPhysics.hh>
#include <G4NeutronTrackingCut.hh>

#include <G4HadronPhysicsQGSP_BERT.hh>

#include "ConfMan.hh"

namespace
{
const auto& confMan = ConfMan::GetInstance();
}

//_____________________________________________________________________________
S2SPhysicsList::S2SPhysicsList(G4int verbose_level)
  : G4VModularPhysicsList()
{
  if(verbose_level > 0){
    G4cout << "<<< Geant4 Physics List simulation engine: S2SPhysicsList"<<G4endl;
    G4cout <<G4endl;
  }
  defaultCutValue = 0.7*CLHEP::mm;
  SetVerboseLevel(verbose_level);

  if(confMan.Get<G4bool>("EM")){
    // EM Physics
    RegisterPhysics(new G4EmStandardPhysics(verbose_level));
    // Synchroton Radiation & GN Physics
    RegisterPhysics(new G4EmExtraPhysics(verbose_level));
  }

  if(confMan.Get<G4bool>("DECAY")){
    // Decays
    RegisterPhysics(new G4DecayPhysics(verbose_level));
  }

  if(confMan.Get<G4bool>("HADRON")){
    // Hadron Elastic scattering
    RegisterPhysics(new G4HadronElasticPhysics(verbose_level));
    // Hadron Physics
    RegisterPhysics(new G4HadronPhysicsQGSP_BERT(verbose_level));
  }

  ///// Others
  // Stopping Physics
  RegisterPhysics(new G4StoppingPhysics(verbose_level));
  // Ion Physics
  RegisterPhysics(new G4IonPhysics(verbose_level));
  // Neutron tracking cut
  RegisterPhysics(new G4NeutronTrackingCut(verbose_level));
}
