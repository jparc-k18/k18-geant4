// -*- C++ -*-

#ifndef DC_SD_HH
#define DC_SD_HH 1

#include <G4VSensitiveDetector.hh>

#include "DCHit.hh"

class G4Step;
class G4HCofThisEvent;
class G4TouchableHistory;

//_____________________________________________________________________________
class DCSD : public G4VSensitiveDetector
{
public:
  DCSD(const G4String& name);
  ~DCSD();

private:
  DCHitsCollection *DCCollection;

public:
  void EndOfEvent(G4HCofThisEvent* HCE);
  void Initialize(G4HCofThisEvent* HCE);
  G4bool ProcessHits(G4Step* aStep, G4TouchableHistory* ROhist);
};

#endif
