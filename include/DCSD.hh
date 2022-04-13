/*
  DCSD.hh
  2007/4  K.Shirotori
*/

#ifndef DCSD_h
#define DCSD_h 1

#include "G4VSensitiveDetector.hh"
#include "DCHit.hh"

class G4Step;
class G4HCofThisEvent;
class G4TouchableHistory;

class DCSD : public G4VSensitiveDetector
{
public:
  DCSD( G4String name );
  ~DCSD();

  void Initialize( G4HCofThisEvent *HCE );
  G4bool ProcessHits( G4Step *aStep, G4TouchableHistory *ROhist );
  void EndOfEvent( G4HCofThisEvent *HCE );

//   void DrawAll() const;
//   void PrintAll() const;
  void clear();

private:
  int EMFlag;
  DCHitsCollection *DCCollection;
};

#endif
