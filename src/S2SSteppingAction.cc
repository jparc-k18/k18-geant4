// -*- C++ -*-

#include "S2SSteppingAction.hh"

#include <algorithm>
#include <cmath>

#include <G4SteppingManager.hh>
#include <G4Track.hh>
#include <G4Step.hh>
#include <G4StepPoint.hh>
#include <G4TrackStatus.hh>
#include <G4VPhysicalVolume.hh>
#include <G4ParticleDefinition.hh>
#include <G4ParticleTypes.hh>
#include <G4SystemOfUnits.hh>

#include <ConfMan.hh>
#include "DCGeomMan.hh"
#include "S2SAnaManager.hh"

namespace
{
  const auto& confMan = ConfMan::GetInstance();
  const auto& geomMan = DCGeomMan::GetInstance();

  // Record the transported primary state at the nominal target plane.  Linear
  // interpolation within the crossing step avoids dependence on step length.
  void RecordTargetTruth(const G4Step* step)
  {
    const G4int generator = confMan.Get<G4int>("Generator");
    if(!confMan.Get<G4bool>("UseK18Beamline")
       || (generator != 6380 && generator != 6381))
      return;
    auto* track = step ? step->GetTrack() : nullptr;
    auto& ana = S2SAnaManager::GetInstance();
    if(!track || track->GetParentID() != 0 || ana.HasTargetTruth())
      return;
    const auto* particle = track->GetDefinition();
    if(!particle)
      return;
    const G4int abs_pdg = std::abs(particle->GetPDGEncoding());
    if(abs_pdg != 211 && abs_pdg != 321)
      return;

    const auto* pre = step->GetPreStepPoint();
    const auto* post = step->GetPostStepPoint();
    if(!pre || !post)
      return;
    const auto pre_position = pre->GetPosition();
    const auto post_position = post->GetPosition();
    const G4double delta_z = post_position.z()-pre_position.z();
    if(std::abs(delta_z) < 1.e-12*mm)
      return;

    const G4double target_z = geomMan.GetGlobalPosition("Target").z()*mm;
    const G4double fraction = (target_z-pre_position.z())/delta_z;
    if(fraction < -1.e-9 || fraction > 1.+1.e-9)
      return;
    const G4double clamped = std::max(0., std::min(1., fraction));
    const auto position = pre_position+
      clamped*(post_position-pre_position);
    const auto momentum = pre->GetMomentum()+
      clamped*(post->GetMomentum()-pre->GetMomentum());
    if(momentum.mag2() <= 0. || std::abs(momentum.z()) < 1.e-12*MeV)
      return;
    const G4double time = pre->GetGlobalTime()+
      clamped*(post->GetGlobalTime()-pre->GetGlobalTime());
    const auto relative =
      position-geomMan.GetGlobalPosition("Target")*mm;
    ana.SetTargetTruthData(
      relative.x(), relative.y(), relative.z(),
      momentum.x()/momentum.z(), momentum.y()/momentum.z(),
      momentum.phi(), momentum.theta(), momentum.mag(), time);
  }
}

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
  RecordTargetTruth(theStep);
  auto theTrack = theStep->GetTrack();
  auto thePrePoint = theStep->GetPreStepPoint();
  auto thePrePVname = thePrePoint->GetPhysicalVolume()->GetName();
  auto partName = theTrack->GetDefinition()->GetParticleName();
  auto thePreMaterial = thePrePoint->GetMaterial();
  if(confMan.Get<G4bool>("StopAtIron")
     && thePreMaterial->GetName() == "Fe"){
    theTrack->SetTrackStatus(fKillTrackAndSecondaries);
  }
}
