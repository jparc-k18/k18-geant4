// -*- C++ -*-

#include "GeSD.hh"

#include <G4HCofThisEvent.hh>
#include <G4Step.hh>
#include <G4SystemOfUnits.hh>

GeSD::GeSD(const G4String& name)
  : G4VSensitiveDetector(name),
    m_collection()
{
  collectionName.insert(name);
}

GeSD::~GeSD()
{
}

void
GeSD::Initialize(G4HCofThisEvent* hce)
{
  m_collection = new GeHitsCollection(SensitiveDetectorName,
                                      collectionName[0]);
  hce->AddHitsCollection(GetCollectionID(0), m_collection);
}

G4bool
GeSD::ProcessHits(G4Step* step, G4TouchableHistory* /* history */)
{
  if(step->GetTotalEnergyDeposit() <= 0.*keV)
    return false;
  m_collection->insert(new GeHit(SensitiveDetectorName, step));
  return true;
}

void
GeSD::EndOfEvent(G4HCofThisEvent* /* hce */)
{
}
