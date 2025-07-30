// -*- C++ -*-

#ifndef PDSD_h
#define PDSD_h 1

#include <G4VSensitiveDetector.hh>

#include "PDHit.hh"

class G4Step;
class G4HCofThisEvent;
class G4TouchableHistory;

//_____________________________________________________________________________
class PDSD : public G4VSensitiveDetector
{
public:
  PDSD(const G4String& name);
  ~PDSD();

private:
  PDHitsCollection* PDCollection;

public:
  void   EndOfEvent(G4HCofThisEvent *HCE);
  void   Initialize(G4HCofThisEvent *HCE);
  G4bool ProcessHits(G4Step *aStep, G4TouchableHistory *ROhist);
};

#endif
