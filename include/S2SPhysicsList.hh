// -*- C++ -*-

#ifndef S2S_PHYSICS_LIST_HH
#define S2S_PHYSICS_LIST_HH 1

#include <globals.hh>
#include <G4VModularPhysicsList.hh>

//_____________________________________________________________________________
class S2SPhysicsList: public G4VModularPhysicsList
{
public:
  S2SPhysicsList(G4int verbose_level=1);
  virtual ~S2SPhysicsList() = default;

  S2SPhysicsList(const S2SPhysicsList&) = delete;
  S2SPhysicsList& operator =(const S2SPhysicsList&) = delete;
};

#endif
