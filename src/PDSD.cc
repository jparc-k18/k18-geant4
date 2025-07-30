// -*- C++ -*-

#include "PDSD.hh"

#include <G4HCofThisEvent.hh>
#include <G4VPhysicalVolume.hh>
#include <G4TouchableHistory.hh>
#include <G4Step.hh>
#include <G4Track.hh>
#include <G4VVisManager.hh>
#include <G4TouchableHandle.hh>
#include <G4SystemOfUnits.hh>

#include "ConfMan.hh"

namespace
{
const auto& confMan = ConfMan::GetInstance();
}

//_____________________________________________________________________________
PDSD::PDSD(const G4String& name)
  : G4VSensitiveDetector(name)
{
  collectionName.insert(name);
}

//_____________________________________________________________________________
PDSD::~PDSD()
{
}

//_____________________________________________________________________________
void
PDSD::Initialize(G4HCofThisEvent *HCE)
{
  PDCollection = new PDHitsCollection(SensitiveDetectorName,
                                        collectionName[0]);
  HCE->AddHitsCollection(GetCollectionID(0), PDCollection);
}

//_____________________________________________________________________________
G4bool
PDSD::ProcessHits(G4Step *aStep, G4TouchableHistory* /* ROhist*/)
{
  const auto preStepPoint = aStep->GetPreStepPoint();
  const auto aTrack = aStep->GetTrack();
  const auto Definition = aTrack->GetDefinition();
  const G4String particleName = Definition->GetParticleName();
  const G4String particleType = Definition->GetParticleType();
  if(preStepPoint->GetStepStatus() != fGeomBoundary)
    return false;
  if(Definition->GetPDGCharge() == 0.)
    return false;
  PDCollection->insert(new PDHit(SensitiveDetectorName, aStep));
  return true;
}

//_____________________________________________________________________________
void
PDSD::EndOfEvent(G4HCofThisEvent* /* HCE */)
{
}
