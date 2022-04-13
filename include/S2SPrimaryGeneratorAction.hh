/*
  "S2SPrimaryGeneratorAction.hh"
  
  Modified by Toshi Gogami , 10Nov2014
*/


#ifndef PRIMARY_GENERATOR_ACTION_H
#define PRIMARY_GENERATOR_ACTION_H

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"
#include "ConfMan.hh"

class G4ParticleGun;
class G4ParticleDefinition;
class G4Event;
class S2SAnalysis;


// ====================================================================
//
// class definition
//
// ====================================================================

class S2SPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {
public:
  S2SPrimaryGeneratorAction(S2SAnalysis *ana,
			    G4double parMomPre=1300,
			    G4double T2Distance=600);
//  S2SPrimaryGeneratorAction(S2SAnalysis*,
//			    ConfMan*,
//			    G4double T2Distance=600);
  S2SPrimaryGeneratorAction(S2SAnalysis*, ConfMan*);
  ~S2SPrimaryGeneratorAction();

  void GeneratePrimaries(G4Event* anEvent);

private:
  // use G4 particle gun
  S2SAnalysis *anaMan_;
  G4ParticleGun* particleGun;
  G4ParticleDefinition* particle;
  G4double parMomPre;
  G4double T2Distance;
  //void SetMom(G4ThreeVector gMom, G4ThreeVector gPos, G4Event* anEvent, G4int evID);
  //void GenParUni(G4Event* anEvent);
  G4double momcent, mombite;
  G4double thetamax;
  ConfMan* confman;
  G4int GenPID;
  G4double beamx, beamy, beamz;
};

#endif
