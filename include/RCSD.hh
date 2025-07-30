// -*- C++ -*-

#ifndef RCSD_h
#define RCSD_h 1

#include <G4VSensitiveDetector.hh>

#include "RCHit.hh"

class G4Step;
class G4HCofThisEvent;
class G4TouchableHistory;

//_____________________________________________________________________________
class RCSD : public G4VSensitiveDetector
{
public:
  RCSD(const G4String& name);
  ~RCSD();

private:
  RCHitsCollection* RCCollection;

public:
  void   EndOfEvent(G4HCofThisEvent *HCE);
  void   Initialize(G4HCofThisEvent *HCE);
  G4bool ProcessHits(G4Step *aStep, G4TouchableHistory *ROhist);
};

#endif
