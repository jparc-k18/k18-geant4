/*
  "WCSD.hh"
 
  Toshiyuki Gogami (2015)
*/

#ifndef WCSD_h
#define WCSD_h 1

#include "G4VSensitiveDetector.hh"
#include "WCHit.hh"
#include "TGraph2D.h"
#include "TF1.h"

class G4Step;
class G4HCofThisEvent;
class G4TouchableHistory;

class WCSD : public G4VSensitiveDetector
{
public:
  WCSD( G4String name );
  ~WCSD();
  
  void Initialize( G4HCofThisEvent *HCE );
  G4bool ProcessHits( G4Step *aStep, G4TouchableHistory *ROhist );
  void EndOfEvent( G4HCofThisEvent *HCE );
  
  //   void DrawAll() const;
  //   void PrintAll() const;
  void clear();
  void DefineXYTable();
  void DefineYTable();
  TGraph2D* GetXYTable() {return xyTable; };
  TF1* GetYTable()       {return yTable;  };
  double GetXYdepFactor(double x, double y){
    return xyTable->Interpolate(x,y);
  };
  double GetYdepFactor(double y){
    return yTable->Eval(y);
  };
  
private:
  int EMFlag;
  WCHitsCollection *WCCollection;
  TGraph2D* xyTable;
  TF1* yTable;
};

#endif
