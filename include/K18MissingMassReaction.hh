// -*- C++ -*-

#ifndef K18_MISSING_MASS_REACTION_HH
#define K18_MISSING_MASS_REACTION_HH 1

#include <G4LorentzVector.hh>
#include <G4ParticleChange.hh>
#include <G4VDiscreteProcess.hh>
#include <globals.hh>

class G4Event;
class G4Track;

namespace K18MissingMassReaction
{
  // The primary action supplies the beam source and requested reaction plane.
  // Kinematics are created only after that beam reaches the plane in Geant4.
  void Arm(const G4Event* event,
           const G4ThreeVector& reaction_vertex,
           const G4LorentzVector& beam_source_p4,
           const G4LorentzVector& beam_source_vertex,
           const G4LorentzVector& beam_pretarget_p4);

  G4bool IsArmedFor(const G4Track& track);
  G4double DistanceToVertex(const G4Track& track);
  G4bool GenerateAtVertex(const G4Track& track,
                          G4ParticleChange& particle_change);
}

// Invokes the configured reaction when the armed beam reaches its vertex.
class K18MissingMassReactionProcess : public G4VDiscreteProcess
{
public:
  K18MissingMassReactionProcess();
  ~K18MissingMassReactionProcess() override = default;

  G4bool IsApplicable(const G4ParticleDefinition&) override;
  G4double PostStepGetPhysicalInteractionLength(
    const G4Track&, G4double, G4ForceCondition*) override;
  G4VParticleChange* PostStepDoIt(const G4Track&, const G4Step&) override;

protected:
  G4double GetMeanFreePath(const G4Track&, G4double,
                           G4ForceCondition*) override;

private:
  G4ParticleChange m_particle_change;
};

#endif
