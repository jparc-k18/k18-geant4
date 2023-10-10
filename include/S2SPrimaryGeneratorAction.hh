/*
  "S2SPrimaryGeneratorAction.hh"

  Modified by Toshi Gogami , 10Nov2014
*/


#ifndef PRIMARY_GENERATOR_ACTION_H
#define PRIMARY_GENERATOR_ACTION_H

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"

class G4ParticleGun;
class G4ParticleDefinition;
class G4Event;


// ====================================================================
//
// class definition
//
// ====================================================================

class S2SPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {
public:
  static G4String ClassName();
  S2SPrimaryGeneratorAction();
  ~S2SPrimaryGeneratorAction();

  void GeneratePrimaries(G4Event* anEvent);
  void GenerateUniform0(G4Event* anEvent);
  void GenerateFocusCheck(G4Event* anEvent);
  void GenerateMonoBeam(G4Event* anEvent);

private:
  G4ParticleGun* particleGun;
  G4ParticleDefinition* m_particle;
  G4double T2Distance;
  //void SetMom(G4ThreeVector gMom, G4ThreeVector gPos, G4Event* anEvent, G4int evID);
  //void GenParUni(G4Event* anEvent);
  G4double momcent, mombite;
  G4double thetamax;
  G4int GenPID;
  G4int generator;
  G4double beamx, beamy, beamz;
};

inline G4String
S2SPrimaryGeneratorAction::ClassName()
{
  static G4String s_name("S2SPrimaryGeneratorAction");
  return s_name;
}

#endif
