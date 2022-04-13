/*
  SlSD.hh
  2007/4  K.Shirotori
*/

#ifndef SlSD_h
#define SlSD_h 1

#include "G4VSensitiveDetector.hh"
#include "SlHit.hh"

class G4Step;
class G4HCofThisEvent;
class G4TouchableHistory;

class SlSD : public G4VSensitiveDetector
{
public:
  SlSD( G4String name );
  ~SlSD();

  void Initialize( G4HCofThisEvent *HCE );
  G4bool ProcessHits( G4Step *aStep, G4TouchableHistory *ROhist );
  void EndOfEvent( G4HCofThisEvent *HCE );

//   void DrawAll() const;
//   void PrintAll() const;
  void clear();

private:
  int EMFlag;
  SlHitsCollection *SlCollection;

};

#endif
