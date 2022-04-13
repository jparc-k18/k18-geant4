/*
  "S2SEventAction.hh"
  
  Modified by Toshiyuki Gogami on 20Nov2014.
*/

#ifndef S2SEventAction_h
#define S2SEventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"
class S2SAnalysis;

class S2SEventAction : public G4UserEventAction
{
public:
  S2SEventAction(S2SAnalysis* ana);
  ~S2SEventAction();
  
private:
  S2SAnalysis* anaMan_;
  
public:
  void  BeginOfEventAction(const G4Event*);
  void  EndOfEventAction(const G4Event*);
  void SetDrawFlag   (G4String val)  {drawFlag = val;};
  void SetPrintModulo(G4int    val)  {printModulo = val;};
  
private:
  G4String  drawFlag;
  G4int printModulo;
  G4int DCcolID;
  time_t start,end;
};

#endif
