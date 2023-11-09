// -*- C++ -*-

#ifndef PRIMARY_GENERATOR_ACTION_H
#define PRIMARY_GENERATOR_ACTION_H

#include <G4VUserPrimaryGeneratorAction.hh>
#include <globals.hh>

class G4ParticleGun;
class G4ParticleDefinition;
class G4Event;

//_____________________________________________________________________________
class S2SPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {
public:
  static G4String ClassName();
  S2SPrimaryGeneratorAction();
  ~S2SPrimaryGeneratorAction();

private:
  G4ParticleGun* m_particleGun;
  G4int          m_generator;

public:
  void GeneratePrimaries(G4Event* anEvent);

private:
  void GenerateDemo(G4Event* anEvent);
  void GenerateUniform0(G4Event* anEvent);
  void GenerateFocusCheck(G4Event* anEvent);
  void GenerateMonochromeBeam(G4Event* anEvent);
  void GenerateUniformSpherical(G4Event* anEvent);
};

inline G4String
S2SPrimaryGeneratorAction::ClassName()
{
  static G4String s_name("S2SPrimaryGeneratorAction");
  return s_name;
}

#endif
