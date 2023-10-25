// -*- C++ -*-

#ifndef S2S_EVENT_ACTION_HH
#define S2S_EVENT_ACTION_HH

#include <G4UserEventAction.hh>
#include <globals.hh>

//_____________________________________________________________________________
class S2SEventAction : public G4UserEventAction
{
public:
  static G4String ClassName();
  S2SEventAction();
  ~S2SEventAction();

public:
  virtual void BeginOfEventAction(const G4Event*);
  virtual void EndOfEventAction(const G4Event*);
};

//_____________________________________________________________________________
inline G4String
S2SEventAction::ClassName()
{
  static G4String s_name("S2SEventAction");
  return s_name;
}

#endif
