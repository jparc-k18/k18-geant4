// -*- C++ -*-

#ifndef WC_SD_HH
#define WC_SD_HH 1

#include <G4VSensitiveDetector.hh>

#include "WCHit.hh"
#include "TGraph2D.h"
#include "TF1.h"

class G4Step;
class G4HCofThisEvent;
class G4TouchableHistory;

//_____________________________________________________________________________
class WCSD : public G4VSensitiveDetector
{
public:
  WCSD(const G4String& name);
  ~WCSD();

private:
  int EMFlag;
  WCHitsCollection *WCCollection;
  TGraph2D* xyTable;
  TF1* yTable;

public:
  void EndOfEvent(G4HCofThisEvent* HCE);
  void Initialize(G4HCofThisEvent* HCE);
  G4bool ProcessHits(G4Step* aStep, G4TouchableHistory* ROhist);

  void DefineXYTable();
  void DefineYTable();
  TGraph2D* GetXYTable(){ return xyTable; }
  TF1*      GetYTable(){ return yTable;  }
  G4double  GetXYdepFactor(G4double x, G4double y)
  { return xyTable->Interpolate(x,y); }
  G4double  GetYdepFactor(G4double y){ return yTable->Eval(y); }
};

#endif
