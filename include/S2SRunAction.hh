/*
  "S2SRunAction.hh"

  T.Gogami, 10 Nov 2014
*/

#ifndef S2SRunAction_h
#define S2SRunAction_h 1

//#include "S2SParamManager.hh"

#include "G4UserRunAction.hh"
#include "ConfMan.hh"
#include "S2SAnaManager.hh"
#include "globals.hh"

#include <TFile.h>
#include <TTree.h>
#include <TObjArray.h>
#include <TH1F.h>

class G4Run;

class S2SRunAction : public G4UserRunAction
{
public:
  S2SRunAction();
  virtual ~S2SRunAction();

  void BeginOfRunAction(const G4Run*);
  void EndOfRunAction(const G4Run*);

};

#endif
