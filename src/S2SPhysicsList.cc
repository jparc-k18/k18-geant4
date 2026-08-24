// -*- C++ -*-

#include "S2SPhysicsList.hh"

#include <cstdlib>
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
#include <G4ParticleDefinition.hh>
#include <G4ParticleTable.hh>
#include <G4ProcessManager.hh>
#include <G4VPhysicsConstructor.hh>

#include "ConfMan.hh"
#include "FuncName.hh"
#include "K18MissingMassReaction.hh"

#define G4MT_physicsVector                                                    \
  ((G4VMPLsubInstanceManager.offset[g4vmplInstanceID]).physicsVector)

namespace
{
const auto& confMan = ConfMan::GetInstance();
}

//_____________________________________________________________________________
S2SPhysicsList::S2SPhysicsList(G4int ver)
  : G4VModularPhysicsList()
{
  defaultCutValue = 0.7*CLHEP::mm;
  verboseLevel = ver;

  // EM Physics
  RegisterPhysics(new G4EmStandardPhysics(ver));
  // Synchrotron Radiation & GN Physics
  RegisterPhysics(new G4EmExtraPhysics(ver));

  // Decays
  RegisterPhysics(new G4DecayPhysics(ver));

  // Hadron Elastic scattering
  RegisterPhysics(new G4HadronElasticPhysics(ver));
  // Hadron Physics
  RegisterPhysics(new G4HadronPhysicsQGSP_BERT(ver));

  ///// Others
  // Stopping Physics
  RegisterPhysics(new G4StoppingPhysics(ver));
  // Ion Physics
  RegisterPhysics(new G4IonPhysics(ver));
  // Neutron tracking cut
  RegisterPhysics(new G4NeutronTrackingCut(ver));
}

//_____________________________________________________________________________
void
S2SPhysicsList::ConstructParticle()
{
  for(auto itr = G4MT_physicsVector->cbegin();
      itr != G4MT_physicsVector->cend(); ++itr)
  {
    (*itr)->ConstructParticle();
  }
}

//_____________________________________________________________________________
void
S2SPhysicsList::ConstructProcess()
{
  // G4AutoLock l(&constructProcessMutex);
  AddTransportation();

  for(auto itr = G4MT_physicsVector->cbegin();
      itr != G4MT_physicsVector->cend(); ++itr)
  {
    auto name = (*itr)->GetPhysicsName();
    if(name == "G4EmStandard" || name == "G4GammaLeptoNuclearPhys"){
      if(!confMan.Get<G4bool>("EM")) continue;
    }else if(name == "Decay"){
      if(!confMan.Get<G4bool>("DECAY")) continue;
    }else{
      if(!confMan.Get<G4bool>("HADRON")) continue;
    }

    if(verboseLevel > 0)
      G4cout << FUNC_NAME << " Construct " << name << G4endl;

    (*itr)->ConstructProcess();
  }

  // The missing-mass reaction is attached to the transported incident beam.
  // Its final state is generated at the target plane selected by the primary
  // action, after ordinary field transport and material interactions.
  if(confMan.Get<G4int>("Generator") == 6381){
    G4String beam_name = confMan.Get<G4String>("ReactionBeamParticle");
    if(beam_name.empty())
      beam_name = "kaon-";
    auto* particle = G4ParticleTable::GetParticleTable()->FindParticle(
      beam_name);
    if(!particle || !particle->GetProcessManager()){
      G4cerr << FUNC_NAME << " cannot attach K18 missing-mass reaction to "
             << beam_name << G4endl;
      std::exit(EXIT_FAILURE);
    }
    particle->GetProcessManager()->AddDiscreteProcess(
      new K18MissingMassReactionProcess());
  }
}

//_____________________________________________________________________________
void
S2SPhysicsList::SetCuts()
{
  SetCutsWithDefault();
}
