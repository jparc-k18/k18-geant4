// -*- C++ -*-

#ifndef K18_REACTION_MODEL_SAMPLER_HH
#define K18_REACTION_MODEL_SAMPLER_HH 1

#include <G4ThreeVector.hh>
#include <globals.hh>

// Sampling helpers shared by K18 missing-mass reaction generators.
namespace K18ReactionModelSampler
{
G4double SampleFlatteCuspMass(G4double mass_min,
                              G4double mass_max,
                              G4double scattering_length_real,
                              G4double scattering_length_imag,
                              G4double momentum_transfer_mev);

G4ThreeVector SampleFermiMomentum(const G4String& model,
                                  G4double sigma_gev);
}

#endif
