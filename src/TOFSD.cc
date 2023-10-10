// -*- C++ -*-

#include "TOFSD.hh"

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
//const G4double DEThreshold = 0.001*keV;
//const G4double DEThreshold = 0.00001*keV;
//const G4double DEThreshold = 0.1*keV;
const G4double DEThreshold = 10.0*keV;
// const double PositionSeparationThreshold = 2.0*cm;
// const double TimeSeparationThreshold     = 5.0*ns;
}

//_____________________________________________________________________________
TOFSD::TOFSD(const G4String& name)
  : G4VSensitiveDetector(name),
    EMFlag(0)
{
  collectionName.insert(name);
  EMFlag = confMan.GetEMFlag();
}

//_____________________________________________________________________________
TOFSD::~TOFSD()
{
}

//_____________________________________________________________________________
void
TOFSD::Initialize(G4HCofThisEvent *HCE)
{
  TOFCollection = new TOFHitsCollection(SensitiveDetectorName,
                                        collectionName[0]);
  HCE->AddHitsCollection(GetCollectionID(0), TOFCollection);
}

//_____________________________________________________________________________
G4bool
TOFSD::ProcessHits(G4Step *aStep, G4TouchableHistory* /* ROhist*/)
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
  TOFCollection->insert(new TOFHit(SensitiveDetectorName, aStep));
  return true;
}

//_____________________________________________________________________________
void
TOFSD::EndOfEvent(G4HCofThisEvent* /* HCE */)
{
}
