// -*- C++ -*-

#include "DCSD.hh"

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
DCSD::DCSD(const G4String& name)
  : G4VSensitiveDetector(name),
    EMFlag(0)
{
  collectionName.insert(name);
  EMFlag = confMan.GetEMFlag();
}

//_____________________________________________________________________________
DCSD::~DCSD()
{
}

//_____________________________________________________________________________
void
DCSD::Initialize(G4HCofThisEvent* HCE)
{
  DCCollection = new DCHitsCollection(SensitiveDetectorName,
                                      collectionName[0]);
  HCE->AddHitsCollection(GetCollectionID(0), DCCollection);
}

//_____________________________________________________________________________
G4bool
DCSD::ProcessHits(G4Step* aStep, G4TouchableHistory* /* ROhist */)
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
  DCCollection->insert(new DCHit(SensitiveDetectorName, aStep));
  return true;
}

//_____________________________________________________________________________
void
DCSD::EndOfEvent(G4HCofThisEvent* /* HCE */)
{
}
