// -*- C++ -*-

#ifndef GE_SD_HH
#define GE_SD_HH 1

#include <G4VSensitiveDetector.hh>

#include "GeHit.hh"

class G4HCofThisEvent;
class G4Step;
class G4TouchableHistory;

class GeSD : public G4VSensitiveDetector
{
public:
  GeSD(const G4String& name);
  virtual ~GeSD();

  virtual void Initialize(G4HCofThisEvent* hce);
  virtual G4bool ProcessHits(G4Step* step, G4TouchableHistory* history);
  virtual void EndOfEvent(G4HCofThisEvent* hce);

private:
  GeHitsCollection* m_collection;
};

#endif
