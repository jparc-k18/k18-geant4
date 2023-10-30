// -*- C++ -*-

#include "S2SSteppingAction.hh"

#include <G4SteppingManager.hh>
#include <G4Track.hh>
#include <G4Step.hh>
#include <G4StepPoint.hh>
#include <G4TrackStatus.hh>
#include <G4VPhysicalVolume.hh>
#include <G4ParticleDefinition.hh>
#include <G4ParticleTypes.hh>

//_____________________________________________________________________________
S2SSteppingAction::S2SSteppingAction()
{
}

//_____________________________________________________________________________
S2SSteppingAction::~S2SSteppingAction()
{
}

//_____________________________________________________________________________
void
S2SSteppingAction::UserSteppingAction(const G4Step* theStep)
{
  auto theTrack = theStep->GetTrack();
  auto thePrePoint = theStep->GetPreStepPoint();
  auto thePrePVname = thePrePoint->GetPhysicalVolume()->GetName();
  auto partName = theTrack->GetDefinition()->GetParticleName();
  auto thePreMaterial = thePrePoint->GetMaterial();
  if(thePreMaterial->GetName() == "Fe"){
    theTrack->SetTrackStatus(fKillTrackAndSecondaries);
  }
}
