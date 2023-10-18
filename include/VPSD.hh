// -*- C++ -*-

#ifndef VP_SD_HH
#define VP_SD_HH

#include "G4VSensitiveDetector.hh"
#include "VPHit.hh"

class G4Step;
class G4HCofThisEvent;
class G4TouchableHistory;

//_____________________________________________________________________________
class VPSD : public G4VSensitiveDetector
{
public:
  VPSD(G4String name);
  ~VPSD();

  void Initialize( G4HCofThisEvent *HCE );
  G4bool ProcessHits( G4Step *aStep, G4TouchableHistory *ROhist );
  void EndOfEvent( G4HCofThisEvent *HCE );

//   void DrawAll() const;
//   void PrintAll() const;
  void clear();

private:
  int EMFlag;
  VPHitsCollection *VPCollection;

};

#endif
