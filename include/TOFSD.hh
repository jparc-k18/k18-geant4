/*
  TOFSD.hh
  2007/4  K.Shirotori
*/

#ifndef TOFSD_h
#define TOFSD_h 1

#include "G4VSensitiveDetector.hh"
#include "TOFHit.hh"

class G4Step;
class G4HCofThisEvent;
class G4TouchableHistory;

class TOFSD : public G4VSensitiveDetector
{
public:
  TOFSD( G4String name );
  ~TOFSD();

  void Initialize( G4HCofThisEvent *HCE );
  G4bool ProcessHits( G4Step *aStep, G4TouchableHistory *ROhist );
  void EndOfEvent( G4HCofThisEvent *HCE );

//   void DrawAll() const;
//   void PrintAll() const;
  void clear();

private:
  int EMFlag;
  TOFHitsCollection *TOFCollection;
};

#endif
