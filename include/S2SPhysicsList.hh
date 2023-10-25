// -*- C++ -*-

#ifndef S2S_PHYSICS_LIST_HH
#define S2S_PHYSICS_LIST_HH 1

#include <globals.hh>
#include <G4VModularPhysicsList.hh>

//_____________________________________________________________________________
class S2SPhysicsList: public G4VModularPhysicsList
{
public:
  static G4String ClassName();
  S2SPhysicsList(G4int verbose_level=1);
  virtual ~S2SPhysicsList() = default;

  S2SPhysicsList(const S2SPhysicsList&) = delete;
  S2SPhysicsList& operator =(const S2SPhysicsList&) = delete;

public:
  virtual void ConstructParticle();
  virtual void ConstructProcess();
  virtual void SetCuts();
};

//_____________________________________________________________________________
inline G4String
S2SPhysicsList::ClassName()
{
  static G4String s_name("S2SPhysicsList");
  return s_name;
}

#endif
