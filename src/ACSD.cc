// -*- C++ -*-

#include "ACSD.hh"

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
ACSD::ACSD(const G4String& name)
  : G4VSensitiveDetector(name)
{
  collectionName.insert(name);
}

//_____________________________________________________________________________
ACSD::~ACSD()
{
}

//_____________________________________________________________________________
void
ACSD::Initialize(G4HCofThisEvent* HCE)
{
  ACCollection = new ACHitsCollection(SensitiveDetectorName,
                                      collectionName[0]);
  HCE->AddHitsCollection(GetCollectionID(0), ACCollection);
}

//_____________________________________________________________________________
G4bool
ACSD::ProcessHits(G4Step* aStep, G4TouchableHistory *ROhist)
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
  ACCollection->insert(new ACHit(SensitiveDetectorName, aStep));
  return true;
}

//_____________________________________________________________________________
void
ACSD::EndOfEvent(G4HCofThisEvent* /* HCE */)
{
}
