/*
  "S2SSteppingAction.cc"

  Modified by Toshiyuki Gogami on 20Nov2014
*/

#include "S2SSteppingAction.hh"
#include "G4SteppingManager.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4TrackStatus.hh"
#include "G4VPhysicalVolume.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"

S2SSteppingAction::S2SSteppingAction()
{}

S2SSteppingAction::~S2SSteppingAction()
{}

void S2SSteppingAction::UserSteppingAction(const G4Step * theStep)
{
  G4Track * theTrack = theStep->GetTrack();

  G4StepPoint * thePrePoint = theStep->GetPreStepPoint();
  G4String thePrePVname = thePrePoint->GetPhysicalVolume()->GetName();
  G4String partName = theTrack->GetDefinition()->GetParticleName();

  //  std::cout<<"physVolume="<<thePrePVname<<"  particleName="<<partName<<std::endl;

  if(thePrePVname=="physQ1Pole" || thePrePVname=="physQ2Pole" || thePrePVname=="physD1Pole" || thePrePVname=="physEnd"){
    theTrack->SetTrackStatus(fKillTrackAndSecondaries);
  }
}
