/*
  ACSD.hh
  2007/4  K.Shirotori
*/

#ifndef ACSD_h
#define ACSD_h 1

#include "G4VSensitiveDetector.hh"
#include "ACHit.hh"

class G4Step;
class G4HCofThisEvent;
class G4TouchableHistory;

class ACSD : public G4VSensitiveDetector
{
public:
  ACSD( G4String name );
  ~ACSD();

  void Initialize( G4HCofThisEvent *HCE );
  G4bool ProcessHits( G4Step *aStep, G4TouchableHistory *ROhist );
  void EndOfEvent( G4HCofThisEvent *HCE );

//   void DrawAll() const;
//   void PrintAll() const;
  void clear();

private:
  ACHitsCollection *ACCollection;
};

#endif
