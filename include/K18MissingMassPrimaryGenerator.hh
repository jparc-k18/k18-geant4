// -*- C++ -*-

#ifndef K18_MISSING_MASS_PRIMARY_GENERATOR_HH
#define K18_MISSING_MASS_PRIMARY_GENERATOR_HH 1

#include <memory>

#include <globals.hh>

class G4Event;
class G4ParticleGun;
class K18BeamlineField;
class K18PhaseSpaceProfile;

// Generates the incident K1.8 beam used by missing-mass studies.  The beam is
// launched upstream and the reaction process creates final-state particles
// only after Geant4 transports it to the sampled target position.
class K18MissingMassPrimaryGenerator
{
public:
  K18MissingMassPrimaryGenerator();
  ~K18MissingMassPrimaryGenerator();

  K18MissingMassPrimaryGenerator(
    const K18MissingMassPrimaryGenerator&) = delete;
  K18MissingMassPrimaryGenerator& operator=(
    const K18MissingMassPrimaryGenerator&) = delete;

  void GeneratePhaseSpaceBeam(G4Event* event, G4ParticleGun* particle_gun);
  void GenerateMissingMassBeam(G4Event* event, G4ParticleGun* particle_gun);

private:
  void GenerateBeam(G4Event* event, G4ParticleGun* particle_gun,
                    G4bool attach_reaction);

  std::unique_ptr<K18PhaseSpaceProfile> m_profile;
  std::unique_ptr<K18BeamlineField> m_backpropagation_field;
};

#endif
