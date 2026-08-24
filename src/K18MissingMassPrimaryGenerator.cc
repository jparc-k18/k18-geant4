// -*- C++ -*-

#include "K18MissingMassPrimaryGenerator.hh"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <exception>
#include <limits>

#include <G4Event.hh>
#include <G4LorentzVector.hh>
#include <G4ParticleDefinition.hh>
#include <G4ParticleGun.hh>
#include <G4ParticleTable.hh>
#include <G4SystemOfUnits.hh>
#include <Randomize.hh>

#include "ConfMan.hh"
#include "DCGeomMan.hh"
#include "GeneratorParticleBranches.hh"
#include "K18BeamlineField.hh"
#include "K18MissingMassReaction.hh"
#include "K18PhaseSpaceProfile.hh"
#include "K18RunControl.hh"
#include "S2SAnaManager.hh"

namespace
{
const auto& confMan = ConfMan::GetInstance();
const auto& geomMan = DCGeomMan::GetInstance();
const auto particleTable = G4ParticleTable::GetParticleTable();

[[noreturn]] void
Fail(const G4String& message)
{
  G4cerr << "[K18MissingMassPrimaryGenerator] " << message << G4endl;
  std::exit(EXIT_FAILURE);
}

G4double
ConfDoubleOr(const G4String& key, G4double fallback)
{
  return confMan.Get<G4String>(key).empty()
    ? fallback : confMan.Get<G4double>(key);
}

G4double
ConfDoubleOr(const G4String& key, const G4String& alias, G4double fallback)
{
  if(!confMan.Get<G4String>(key).empty())
    return confMan.Get<G4double>(key);
  return confMan.Get<G4String>(alias).empty()
    ? fallback : confMan.Get<G4double>(alias);
}

G4bool
ConfBoolOr(const G4String& key, G4bool fallback)
{
  return confMan.Get<G4String>(key).empty()
    ? fallback : confMan.Get<G4bool>(key);
}

G4bool
ConfBoolOr(const G4String& key, const G4String& alias, G4bool fallback)
{
  if(!confMan.Get<G4String>(key).empty())
    return confMan.Get<G4bool>(key);
  return confMan.Get<G4String>(alias).empty()
    ? fallback : confMan.Get<G4bool>(alias);
}

G4String
ConfStringOr(const G4String& key, const G4String& fallback)
{
  const G4String value = confMan.Get<G4String>(key);
  return value.empty() ? fallback : value;
}

G4String
ConfStringOr(const G4String& key, const G4String& alias,
             const G4String& fallback)
{
  const G4String value = confMan.Get<G4String>(key);
  if(!value.empty())
    return value;
  const G4String alias_value = confMan.Get<G4String>(alias);
  return alias_value.empty() ? fallback : alias_value;
}

G4double
BeamDouble(G4bool reaction, const G4String& reaction_key,
           const G4String& phase_space_key, G4double fallback)
{
  return reaction
    ? ConfDoubleOr(reaction_key, phase_space_key, fallback)
    : ConfDoubleOr(phase_space_key, fallback);
}

G4bool
BeamBool(G4bool reaction, const G4String& reaction_key,
         const G4String& phase_space_key, G4bool fallback)
{
  return reaction
    ? ConfBoolOr(reaction_key, phase_space_key, fallback)
    : ConfBoolOr(phase_space_key, fallback);
}

G4String
BeamString(G4bool reaction, const G4String& reaction_key,
           const G4String& phase_space_key, const G4String& fallback)
{
  return reaction
    ? ConfStringOr(reaction_key, phase_space_key, fallback)
    : ConfStringOr(phase_space_key, fallback);
}

G4String
BackpropagationFieldPrefix()
{
  if(confMan.Get<G4int>("Generator") != 6381)
    return "K18PhaseSpaceBackProp";
  static const char* suffixes[] = {
    "GlobalScale", "Q10Scale", "Q11Scale", "D4Scale", "Q12Scale",
    "Q13Scale", "K18Q10Scale", "K18Q11Scale", "K18D4Scale",
    "K18Q12Scale", "K18Q13Scale"
  };
  for(const char* suffix : suffixes){
    if(!confMan.Get<G4String>("ReactionBeamBackProp" +
                              G4String(suffix)).empty())
      return "ReactionBeamBackProp";
  }
  return "K18PhaseSpaceBackProp";
}

G4double
Shoot(G4double mean, G4double sigma, G4double half_width)
{
  if(half_width > 0.)
    return G4RandFlat::shoot(mean-half_width, mean+half_width);
  if(sigma > 0.)
    return G4RandGauss::shoot(mean, sigma);
  return mean;
}

struct BeamSample
{
  G4double momentum;
  G4double x;
  G4double y;
  G4double z;
  G4double u;
  G4double v;
};

struct TransportState
{
  G4ThreeVector position;
  G4ThreeVector momentum; // GeV/c in the K1.8 internal frame.
};

G4double
UpstreamCoordinate(const G4ThreeVector& position)
{
  return (position-K18BeamlineFrame::Vin()).dot(K18BeamlineFrame::Axis(0.));
}

TransportState
Derivative(const TransportState& state, const K18BeamlineField& field,
           G4double charge)
{
  TransportState derivative;
  if(state.momentum.mag2() <= 0.)
    return derivative;

  const G4ThreeVector direction = state.momentum.unit();
  G4double point[4] = {state.position.x(), state.position.y(),
                       state.position.z(), 0.};
  G4double value[3] = {0., 0., 0.};
  field.GetFieldValue(point, value);
  const G4ThreeVector field_tesla(value[0]/tesla, value[1]/tesla,
                                  value[2]/tesla);
  derivative.position = direction;
  derivative.momentum =
    charge*0.000299792458*direction.cross(field_tesla);
  return derivative;
}

void
StepRungeKutta(TransportState& state, G4double step,
               const K18BeamlineField& field, G4double charge)
{
  const G4double momentum = state.momentum.mag();
  const auto k1 = Derivative(state, field, charge);
  const TransportState s2{state.position+0.5*step*k1.position,
                          state.momentum+0.5*step*k1.momentum};
  const auto k2 = Derivative(s2, field, charge);
  const TransportState s3{state.position+0.5*step*k2.position,
                          state.momentum+0.5*step*k2.momentum};
  const auto k3 = Derivative(s3, field, charge);
  const TransportState s4{state.position+step*k3.position,
                          state.momentum+step*k3.momentum};
  const auto k4 = Derivative(s4, field, charge);

  state.position += step/6.*(k1.position+2.*k2.position+
                             2.*k3.position+k4.position);
  state.momentum += step/6.*(k1.momentum+2.*k2.momentum+
                             2.*k3.momentum+k4.momentum);
  if(state.momentum.mag2() > 0.)
    state.momentum = momentum*state.momentum.unit();
}

G4bool
BackPropagate(TransportState& state, G4double source_coordinate,
              G4double max_path, G4double step,
              const K18BeamlineField& field, G4double charge)
{
  step = std::abs(step);
  if(step <= 0.)
    Fail("back-propagation step must be positive");
  if(max_path <= 0.)
    Fail("back-propagation maximum path must be positive");

  const G4double momentum = state.momentum.mag();
  TransportState previous = state;
  G4double previous_delta =
    UpstreamCoordinate(previous.position)-source_coordinate;
  if(previous_delta <= 0.)
    return true;

  const G4int max_steps = static_cast<G4int>(std::ceil(max_path/step));
  for(G4int i=0; i<max_steps; ++i){
    TransportState current = previous;
    StepRungeKutta(current, -step, field, charge);
    const G4double current_delta =
      UpstreamCoordinate(current.position)-source_coordinate;
    if(current_delta <= 0.){
      const G4double denominator = previous_delta-current_delta;
      const G4double fraction = denominator != 0.
        ? previous_delta/denominator : 0.;
      state.position = previous.position+
        fraction*(current.position-previous.position);
      state.momentum = previous.momentum+
        fraction*(current.momentum-previous.momentum);
      if(state.momentum.mag2() > 0.)
        state.momentum = momentum*state.momentum.unit();
      return true;
    }
    previous = current;
    previous_delta = current_delta;
  }
  return false;
}
}

K18MissingMassPrimaryGenerator::K18MissingMassPrimaryGenerator()
  : m_profile(),
    m_backpropagation_field(
      new K18BeamlineField(BackpropagationFieldPrefix()))
{
  const G4bool reaction = confMan.Get<G4int>("Generator") == 6381;
  const G4String profile_name = BeamString(
    reaction, "ReactionBeamProfile", "K18PhaseSpaceProfile", "");
  if(profile_name.empty())
    return;

  const G4String profile_tree = BeamString(
    reaction, "ReactionBeamProfileTree", "K18PhaseSpaceProfileTree", "");
  const G4double default_profile_p = BeamDouble(
    reaction, "ReactionBeamP", "K18PhaseSpaceP",
    ConfDoubleOr("PK18", 1.4));
  const G4double profile_p = BeamDouble(
    reaction, "ReactionBeamProfileP", "K18PhaseSpaceProfileP",
    default_profile_p);
  const G4String particle_name = BeamString(
    reaction, "ReactionBeamParticle", "K18PhaseSpaceParticle", "kaon-");
  const auto* profile_particle = particleTable->FindParticle(particle_name);
  if(!profile_particle)
    Fail("unknown profile beam particle: " + particle_name);
  const G4double profile_particle_charge =
    profile_particle->GetPDGCharge()/eplus;
  if(profile_particle_charge == 0.)
    Fail("profile beam particle must be charged: " + particle_name);
  const G4int profile_charge = profile_particle_charge > 0. ? 1 : -1;
  const G4double profile_z = BeamDouble(
    reaction, "ReactionBeamProfileTargetZ", "K18PhaseSpaceProfileTargetZ",
    1599.6);
  const G4double profile_chisqr_max = BeamDouble(
    reaction, "ReactionBeamProfileChiSqMax",
    "K18PhaseSpaceProfileChiSqMax", 5.);
  const G4double profile_max_tracks = BeamDouble(
    reaction, "ReactionBeamProfileMaxTracks",
    "K18PhaseSpaceProfileMaxTracks", 200000.);
  if(!std::isfinite(profile_max_tracks) || profile_max_tracks < 1.
     || std::floor(profile_max_tracks) != profile_max_tracks
     || profile_max_tracks > std::numeric_limits<G4long>::max())
    Fail("K18PhaseSpaceProfileMaxTracks must be a positive integer");
  try {
    m_profile.reset(new K18PhaseSpaceProfile(
      profile_name, profile_tree, profile_z, profile_chisqr_max,
      profile_charge, profile_p,
      static_cast<G4long>(std::lround(profile_max_tracks))));
  } catch(const std::exception& error) {
    Fail(error.what());
  }
}

K18MissingMassPrimaryGenerator::~K18MissingMassPrimaryGenerator() = default;

void
K18MissingMassPrimaryGenerator::GeneratePhaseSpaceBeam(
  G4Event* event, G4ParticleGun* particle_gun)
{
  GenerateBeam(event, particle_gun, false);
}

void
K18MissingMassPrimaryGenerator::GenerateMissingMassBeam(
  G4Event* event, G4ParticleGun* particle_gun)
{
  GenerateBeam(event, particle_gun, true);
}

void
K18MissingMassPrimaryGenerator::GenerateBeam(G4Event* event,
                                             G4ParticleGun* particle_gun,
                                             G4bool attach_reaction)
{
  if(!event || !particle_gun)
    Fail("event and particle gun must be valid");

  const G4String particle_name = BeamString(
    attach_reaction, "ReactionBeamParticle", "K18PhaseSpaceParticle",
    "kaon-");
  auto* particle = particleTable->FindParticle(particle_name);
  if(!particle)
    Fail("unknown incident particle: " + particle_name);

  G4double momentum_mean = BeamDouble(
    attach_reaction,
    "ReactionBeamP", "K18PhaseSpaceP",
    ConfDoubleOr("PK18", 1.4))*GeV;
  G4double momentum_sigma = BeamDouble(
    attach_reaction,
    "ReactionBeamPSigma", "K18PhaseSpacePSigma", 0.)*GeV;
  const G4double relative_sigma = BeamDouble(
    attach_reaction,
    "ReactionBeamDPOverP", "K18PhaseSpaceDPOverP", 0.);
  if(momentum_sigma <= 0. && relative_sigma > 0.)
    momentum_sigma = momentum_mean*relative_sigma;
  const G4double momentum_half_width = BeamDouble(
    attach_reaction,
    "ReactionBeamPHalfWidth", "K18PhaseSpacePHalfWidth", 0.)*GeV;
  const G4double momentum_min = BeamDouble(
    attach_reaction,
    "ReactionBeamPMin", "K18PhaseSpacePMin", 0.)*GeV;
  const G4double momentum_max = BeamDouble(
    attach_reaction,
    "ReactionBeamPMax", "K18PhaseSpacePMax", 0.)*GeV;
  const G4bool reject_outside = BeamBool(
    attach_reaction,
    "ReactionBeamRejectOutsideRange",
    "K18PhaseSpaceRejectOutsideRange", false);
  const G4bool use_profile_momentum = BeamBool(
    attach_reaction,
    "ReactionBeamUseProfileMomentum",
    "K18PhaseSpaceUseProfileMomentum", false);
  if(momentum_mean <= 0.)
    Fail("incident momentum must be positive");
  if(use_profile_momentum && !m_profile)
    Fail("profile momentum requested without a beam profile");

  auto shoot_momentum = [&]()
  {
    for(G4int i=0; i<1000000; ++i){
      G4double momentum = Shoot(momentum_mean, momentum_sigma,
                                momentum_half_width);
      if(reject_outside){
        if(momentum_min > 0. && momentum < momentum_min)
          continue;
        if(momentum_max > momentum_min && momentum > momentum_max)
          continue;
      } else {
        if(momentum_min > 0.)
          momentum = std::max(momentum, momentum_min);
        if(momentum_max > momentum_min)
          momentum = std::min(momentum, momentum_max);
      }
      return momentum;
    }
    Fail("cannot sample momentum inside the configured range");
  };

  auto shoot_sample = [&](G4long sample_index)
  {
    G4double momentum = use_profile_momentum
      ? momentum_mean : shoot_momentum();
    if(m_profile){
      const G4double profile_momentum = BeamDouble(
        attach_reaction,
        "ReactionBeamProfileP", "K18PhaseSpaceProfileP",
        momentum/GeV);
      const auto profile = m_profile->Shoot(
        particle->GetPDGCharge()/eplus, profile_momentum, sample_index);
      if(use_profile_momentum && !profile.has_p)
        Fail("selected beam profile does not provide p_gev");
      if(use_profile_momentum)
        momentum = profile.p_gev*GeV;
      return BeamSample{
        momentum, profile.x_mm*mm, profile.y_mm*mm,
        Shoot(BeamDouble(attach_reaction, "ReactionVertexZMean",
                         "K18PhaseSpaceZMean", 0.)*mm,
              BeamDouble(attach_reaction, "ReactionVertexZSigma",
                         "K18PhaseSpaceZSigma", 0.)*mm,
              BeamDouble(attach_reaction, "ReactionVertexZHalfWidth",
                         "K18PhaseSpaceZHalfWidth", 0.)*mm),
        profile.u, profile.v};
    }
    return BeamSample{
      momentum,
      Shoot(BeamDouble(attach_reaction, "ReactionBeamXMean",
                       "K18PhaseSpaceXMean", 0.)*mm,
            BeamDouble(attach_reaction, "ReactionBeamXSigma",
                       "K18PhaseSpaceXSigma", 30.)*mm,
            BeamDouble(attach_reaction, "ReactionBeamXHalfWidth",
                       "K18PhaseSpaceXHalfWidth", 0.)*mm),
      Shoot(BeamDouble(attach_reaction, "ReactionBeamYMean",
                       "K18PhaseSpaceYMean", 0.)*mm,
            BeamDouble(attach_reaction, "ReactionBeamYSigma",
                       "K18PhaseSpaceYSigma", 5.)*mm,
            BeamDouble(attach_reaction, "ReactionBeamYHalfWidth",
                       "K18PhaseSpaceYHalfWidth", 0.)*mm),
      Shoot(BeamDouble(attach_reaction, "ReactionVertexZMean",
                       "K18PhaseSpaceZMean", 0.)*mm,
            BeamDouble(attach_reaction, "ReactionVertexZSigma",
                       "K18PhaseSpaceZSigma", 0.)*mm,
            BeamDouble(attach_reaction, "ReactionVertexZHalfWidth",
                       "K18PhaseSpaceZHalfWidth", 0.)*mm),
      Shoot(BeamDouble(attach_reaction, "ReactionBeamUMean",
                       "K18PhaseSpaceUMean", 0.),
            BeamDouble(attach_reaction, "ReactionBeamUSigma",
                       "K18PhaseSpaceUSigma", 0.),
            BeamDouble(attach_reaction, "ReactionBeamUHalfWidth",
                       "K18PhaseSpaceUHalfWidth", 0.)),
      Shoot(BeamDouble(attach_reaction, "ReactionBeamVMean",
                       "K18PhaseSpaceVMean", 0.),
            BeamDouble(attach_reaction, "ReactionBeamVSigma",
                       "K18PhaseSpaceVSigma", 0.),
            BeamDouble(attach_reaction, "ReactionBeamVHalfWidth",
                       "K18PhaseSpaceVHalfWidth", 0.))};
  };

  const G4ThreeVector target = geomMan.GetGlobalPosition("Target")*mm;
  const G4double target_coordinate = geomMan.GetLocalZ("K18Target")*mm;
  const K18BeamlineFrame frame(target, target_coordinate);
  G4double source_coordinate = BeamDouble(
    attach_reaction,
    "ReactionBeamSourceOffset", "K18PhaseSpaceSourceOffset", -50.)*mm;
  const G4String source_plane = BeamString(
    attach_reaction,
    "ReactionBeamSourcePlane", "K18PhaseSpaceSourcePlane", "");
  if(source_plane == "VI"){
    source_coordinate = 0.;
  } else if(source_plane == "Q10Upstream" || source_plane == "Q10Front" ||
            source_plane == "Q10Entrance"){
    source_coordinate = K18BeamlineFrame::Q10UpstreamOffset();
  } else if(!source_plane.empty()) {
    Fail("unknown beam source plane: " + source_plane);
  }

  const G4double step = BeamDouble(
    attach_reaction,
    "ReactionBeamBackPropStep", "K18PhaseSpaceBackPropStep", 1.)*mm;
  const G4double max_path = BeamDouble(
    attach_reaction,
    "ReactionBeamBackPropMaxPath", "K18PhaseSpaceBackPropMaxPath", 12000.)*mm;
  const G4int max_trials = std::max<G4int>(1, std::lround(BeamDouble(
    attach_reaction,
    "ReactionBeamBackPropMaxTrial", "K18PhaseSpaceBackPropMaxTrial", 1000.)));
  const G4double charge = particle->GetPDGCharge()/eplus;
  const G4long global_event = K18RunControl::GlobalEventIndex(event);
  BeamSample sample = shoot_sample(global_event);
  TransportState source_state;
  G4bool accepted = false;
  const G4bool use_target_plane = attach_reaction
    || ConfBoolOr("K18PhaseSpaceUseTargetPlane", true);
  if(use_target_plane){
    for(G4int trial=0; trial<max_trials; ++trial){
      const G4ThreeVector vertex =
        target+G4ThreeVector(sample.x, sample.y, sample.z);
      const G4ThreeVector internal_direction = frame.S2SToInternalVector(
        G4ThreeVector(sample.u, sample.v, 1.).unit()).unit();
      source_state = TransportState{
        frame.S2SToInternalPoint(vertex),
        sample.momentum/GeV*internal_direction};
      if(BackPropagate(source_state, source_coordinate, max_path, step,
                       *m_backpropagation_field, charge)){
        accepted = true;
        break;
      }
      sample = shoot_sample(global_event+trial+1);
    }
  } else {
    // Source-plane mode interprets x/y/u/v in the incoming K1.8 frame.
    source_state.position = K18BeamlineFrame::Vin()
      +source_coordinate*K18BeamlineFrame::Axis(0.)
      +sample.x*K18BeamlineFrame::XAxis(0.)
      +sample.y*K18BeamlineFrame::YAxis();
    source_state.momentum = sample.momentum/GeV
      *(K18BeamlineFrame::Axis(0.)
        +sample.u*K18BeamlineFrame::XAxis(0.)
        +sample.v*K18BeamlineFrame::YAxis()).unit();
    accepted = true;
  }
  if(!accepted)
    Fail("cannot back-propagate a sampled target state to the source plane");

  const G4ThreeVector source_position =
    frame.InternalToS2SPoint(source_state.position);
  const G4ThreeVector source_direction =
    frame.InternalToS2SVector(source_state.momentum.unit()).unit();
  const G4ThreeVector target_direction =
    G4ThreeVector(sample.u, sample.v, 1.).unit();
  const G4double mass = particle->GetPDGMass();
  const G4double energy =
    std::sqrt(sample.momentum*sample.momentum+mass*mass);
  const G4LorentzVector source_momentum(
    sample.momentum*source_direction, energy);
  const G4LorentzVector source_vertex(source_position, 0.);
  const G4LorentzVector pretarget_momentum(
    sample.momentum*target_direction, energy);
  const G4ThreeVector reaction_position =
    target+G4ThreeVector(sample.x, sample.y, sample.z);

  // Launch the incident beam upstream.  The reaction process creates the
  // final state when the transported track reaches the target.
  particle_gun->SetParticleDefinition(particle);
  particle_gun->SetParticleMomentumDirection(source_direction);
  particle_gun->SetParticleEnergy(energy-mass);
  particle_gun->SetParticlePosition(source_position);
  particle_gun->GeneratePrimaryVertex(event);
  if(attach_reaction){
    K18MissingMassReaction::Arm(
      event, reaction_position, source_momentum, source_vertex,
      pretarget_momentum);
    return;
  }

  auto& analysis = S2SAnaManager::GetInstance();
  const G4ThreeVector relative = source_position-target;
  const G4double u = std::abs(source_direction.z()) > 1.e-12
    ? source_direction.x()/source_direction.z() : 0.;
  const G4double v = std::abs(source_direction.z()) > 1.e-12
    ? source_direction.y()/source_direction.z() : 0.;
  analysis.SetPrimaryParticle(
    0, particle->GetPDGEncoding(), source_momentum, source_vertex);
  analysis.SetPrimaryData(
    relative.x(), relative.y(), relative.z(), u, v,
    source_direction.phi(), source_direction.theta(), sample.momentum,
    sample.momentum, particle->GetPDGEncoding());
  analysis.SetGeneratedParticle(
    GeneratorParticleBranches::kPrimPi,
    static_cast<G4int>(GeneratorParticleBranches::ParticleId::None),
    particle->GetPDGEncoding(), source_momentum, source_vertex);
}
