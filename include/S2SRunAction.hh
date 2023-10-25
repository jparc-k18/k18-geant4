// -*- C++ -*-

#ifndef S2S_RUN_ACTION_HH
#define S2S_RUN_ACTION_HH 1

#include <G4UserRunAction.hh>

#include <G4String.hh>
#include <globals.hh>

class G4Run;

//_____________________________________________________________________________
class S2SRunAction : public G4UserRunAction
{
public:
  static G4String ClassName();
  S2SRunAction();
  virtual ~S2SRunAction();

public:
  virtual void BeginOfRunAction(const G4Run*);
  virtual void EndOfRunAction(const G4Run*);
};

//_____________________________________________________________________________
inline G4String
S2SRunAction::ClassName()
{
  static G4String s_name("S2SRunAction");
  return s_name;
}

#endif
