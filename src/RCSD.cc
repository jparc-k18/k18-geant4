// -*- C++ -*-

#include "RCSD.hh"

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
RCSD::RCSD(const G4String& name)
  : G4VSensitiveDetector(name)
{
  collectionName.insert(name);
}

//_____________________________________________________________________________
RCSD::~RCSD()
{
}

//_____________________________________________________________________________
void
RCSD::Initialize(G4HCofThisEvent *HCE)
{
  RCCollection = new RCHitsCollection(SensitiveDetectorName,
                                        collectionName[0]);
  HCE->AddHitsCollection(GetCollectionID(0), RCCollection);
}

//_____________________________________________________________________________
G4bool
RCSD::ProcessHits(G4Step *aStep, G4TouchableHistory* /* ROhist*/)
{
  const auto preStepPoint = aStep->GetPreStepPoint();
  const auto aTrack = aStep->GetTrack();
  const auto Definition = aTrack->GetDefinition();
  const G4String particleName = Definition->GetParticleName();
  const G4String particleType = Definition->GetParticleType();
  // if(preStepPoint->GetStepStatus() != fGeomBoundary)
  //   return false;
  // if(Definition->GetPDGCharge() == 0.)
  //   return false;
  RCCollection->insert(new RCHit(SensitiveDetectorName, aStep));
  return true;
}

//_____________________________________________________________________________
void
RCSD::EndOfEvent(G4HCofThisEvent* /* HCE */)
{
}
