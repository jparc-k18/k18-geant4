/*
  "S2SSteppingAction.hh"
  
  Modified by Toshiyuki Gogami on 20Nov2014
*/

#ifndef S2SSteppingAction_H
#define S2SSteppingAction_H 1

#include "globals.hh"
#include "G4UserSteppingAction.hh"

class S2SSteppingAction : public G4UserSteppingAction
{
  public:
    S2SSteppingAction();
    virtual ~S2SSteppingAction();

    virtual void UserSteppingAction(const G4Step*);
};

#endif

