/*
  "S2SPhysicsList.hh"
  
  2007/4  K.Shirotori
  Modified by Toshi Gogami , 10Nov2014
*/

#ifndef PhysicsList_h
#define PhysicsList_h 1

#include "G4VUserPhysicsList.hh"
//#include "globals.hh"
#include "G4Decay.hh"
#include "ConfMan.hh"


class S2SPhysicsList: public G4VUserPhysicsList
{
public:
  S2SPhysicsList();
  S2SPhysicsList(ConfMan*);
  ~S2SPhysicsList();
  
protected:
  // Construct particle and physics
  void ConstructParticle();
  void ConstructProcess();
  void SetCuts();
  
protected:
  // these methods Construct particles 
  void ConstructBosons();
  void ConstructLeptons();
  void ConstructMesons();
  void ConstructBaryons();
  void ConstructHeavyIon();
  
private:
  // these methods Construct physics processes and register them
  void AddTransportation();
  void ConstructDecay();
  void ConstructEM();
  void ConstructHadronic();
  
private:
  //S2SParamManager* paramMan;
  G4int emflag;
  G4int decayflag;
  G4int hadronflag;
  G4double defaultCutValue;
  
};

#endif

 
