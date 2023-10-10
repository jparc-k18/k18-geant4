// -*- C++ -*-

#ifndef ACTION_INITIALIZATION_HH
#define ACTION_INITIALIZATION_HH

#include "G4VUserActionInitialization.hh"

//_____________________________________________________________________________
class S2SActionInitialization : public G4VUserActionInitialization
{
public:
  S2SActionInitialization();
  ~S2SActionInitialization() override;

  void BuildForMaster() const override;
  void Build() const override;
};

#endif
