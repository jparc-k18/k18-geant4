// -*- C++ -*-

#ifndef TOFSD_h
#define TOFSD_h 1

#include <G4VSensitiveDetector.hh>

#include "TOFHit.hh"

class G4Step;
class G4HCofThisEvent;
class G4TouchableHistory;

//_____________________________________________________________________________
class TOFSD : public G4VSensitiveDetector
{
public:
  TOFSD(const G4String& name);
  ~TOFSD();

private:
  TOFHitsCollection* TOFCollection;
  int EMFlag;

public:
  void   EndOfEvent(G4HCofThisEvent *HCE);
  void   Initialize(G4HCofThisEvent *HCE);
  G4bool ProcessHits(G4Step *aStep, G4TouchableHistory *ROhist);
};

#endif
