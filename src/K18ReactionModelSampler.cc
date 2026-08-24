// -*- C++ -*-

#include "K18ReactionModelSampler.hh"

#include <cmath>
#include <cstdlib>

#include <G4SystemOfUnits.hh>
#include <G4ios.hh>
#include <Randomize.hh>

#include "FdComplex.hh"
#include "FermiMotion.hh"

namespace K18ReactionModelSampler
{
G4double
SampleFlatteCuspMass(G4double mass_min,
                     G4double mass_max,
                     G4double scattering_length_real,
                     G4double scattering_length_imag,
                     G4double momentum_transfer_mev)
{
  if(!(mass_max > mass_min)){
    G4cerr << "SampleFlatteCuspMass requires mass_max > mass_min"
           << G4endl;
    std::exit(EXIT_FAILURE);
  }

  // Scan the configured interval once to obtain a conservative rejection
  // envelope for the configured Flatte-like line shape.
  G4double density_max = 0.;
  for(G4int i=0; i<200; ++i){
    const G4double mass = mass_min
      + (mass_max-mass_min)*static_cast<G4double>(i)/199.;
    const G4double density = FdComplex::f_single(
      mass/CLHEP::MeV, scattering_length_real,
      scattering_length_imag, momentum_transfer_mev);
    if(std::isfinite(density) && density > density_max)
      density_max = density;
  }
  if(!(density_max > 0.)){
    G4cerr << "SampleFlatteCuspMass found no positive density in "
           << mass_min/CLHEP::MeV << ":" << mass_max/CLHEP::MeV
           << " MeV" << G4endl;
    std::exit(EXIT_FAILURE);
  }

  for(G4int trial=0; trial<1000000; ++trial){
    const G4double mass = G4RandFlat::shoot(mass_min, mass_max);
    const G4double density = FdComplex::f_single(
      mass/CLHEP::MeV, scattering_length_real,
      scattering_length_imag, momentum_transfer_mev);
    if(std::isfinite(density) && density > 0.
       && G4RandFlat::shoot(0., 1.3*density_max) <= density)
      return mass;
  }

  G4cerr << "SampleFlatteCuspMass exceeded rejection-sampling trials"
         << G4endl;
  std::exit(EXIT_FAILURE);
}

G4ThreeVector
SampleFermiMomentum(const G4String& model,
                    G4double sigma_gev)
{
  G4ThreeVector momentum;
  if(model == "gaussian" || model == "Gaussian"){
    momentum.set(G4RandGauss::shoot(0., sigma_gev),
                 G4RandGauss::shoot(0., sigma_gev),
                 G4RandGauss::shoot(0., sigma_gev));
  } else if(model == "deuteron" || model == "Deuteron"){
    const TVector3 sampled = FermiMotion::GetMomentum();
    momentum.set(sampled.x(), sampled.y(), sampled.z());
  } else {
    G4cerr << "K18ReactionModelSampler: unknown ReactionFermiModel="
           << model << G4endl;
    std::exit(EXIT_FAILURE);
  }
  return momentum;
}
}
