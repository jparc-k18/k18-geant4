// -*- C++ -*-

#include "K18MissingMassReaction.hh"

#include <algorithm>
#include <cerrno>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <vector>

#include <G4DynamicParticle.hh>
#include <G4Event.hh>
#include <G4EventManager.hh>
#include <G4ParticleDefinition.hh>
#include <G4ParticleTable.hh>
#include <G4RandomDirection.hh>
#include <G4Step.hh>
#include <G4SystemOfUnits.hh>
#include <G4Track.hh>
#include <G4VPhysicalVolume.hh>
#include <Randomize.hh>

#include "ConfMan.hh"
#include "DCGeomMan.hh"
#include "K18ReactionModelSampler.hh"
#include "GeneratorParticleBranches.hh"
#include "K18RunControl.hh"
#include "S2SAnaManager.hh"

namespace
{
  const auto& confMan = ConfMan::GetInstance();
  const auto& geomMan = DCGeomMan::GetInstance();
  const auto particleTable = G4ParticleTable::GetParticleTable();
  namespace GenBranch = GeneratorParticleBranches;

  struct PendingReaction
  {
    G4bool armed = false;
    G4int event_id = -1;
    G4int approach_steps = 0;
    G4double vertex_z = 0.;
    G4ThreeVector reaction_vertex;
    G4LorentzVector beam_source_p4;
    G4LorentzVector beam_source_vertex;
    G4LorentzVector beam_pretarget_p4;
  };

  // G4S2S currently uses the serial run manager.  thread_local also keeps the
  // event contract correct if the action/process is later used by worker
  // threads.
  thread_local PendingReaction pending;

  G4double
  ConfDoubleOr(const G4String& key, G4double fallback)
  {
    const G4String raw = confMan.Get<G4String>(key);
    return raw.empty() ? fallback : confMan.Get<G4double>(key);
  }

  G4bool
  ConfBoolOr(const G4String& key, G4bool fallback)
  {
    const G4String raw = confMan.Get<G4String>(key);
    return raw.empty() ? fallback : confMan.Get<G4bool>(key);
  }

  G4String
  ConfStringOr(const G4String& key, const G4String& fallback)
  {
    const G4String raw = confMan.Get<G4String>(key);
    return raw.empty() ? fallback : raw;
  }

  [[noreturn]] void
  Fail(const G4String& message)
  {
    G4cerr << "K18MissingMassReaction: " << message << G4endl;
    std::exit(EXIT_FAILURE);
  }

  void
  ApplyReactionKinematicsSeed()
  {
    const G4String raw = confMan.Get<G4String>(
      "ReactionKinematicsSeedBase");
    if(raw.empty())
      return;
    errno = 0;
    char* end = nullptr;
    const G4long base = std::strtol(raw.c_str(), &end, 10);
    if(errno != 0 || end == raw.c_str() || *end != '\0' || base <= 0)
      Fail("ReactionKinematicsSeedBase must be a positive integer");
    const auto* event = G4EventManager::GetEventManager()->GetConstCurrentEvent();
    if(!event)
      Fail("cannot seed reaction kinematics without a current G4 event");
    const G4long global_event = K18RunControl::GlobalEventIndex(event);
    if(global_event < 0 || base > std::numeric_limits<G4long>::max()-global_event)
      Fail("ReactionKinematicsSeedBase + global event is out of range");
    const G4long seed = base + global_event;
    G4Random::setTheSeed(seed);
    if(event->GetEventID() == 0)
      G4cout << "[K18MissingMassReaction] reaction-kinematics seed base=" << base
             << " global_event=" << global_event
             << " first_seed=" << seed << G4endl;
  }

  G4ParticleDefinition*
  FindParticleOrDie(const G4String& name)
  {
    auto* particle = particleTable->FindParticle(name);
    if(!particle)
      Fail("unknown particle: " + name);
    return particle;
  }

  G4double
  ParticleMassOrConfig(const G4String& particle_key,
                       const G4String& mass_key,
                       const G4String& fallback_particle)
  {
    const G4double configured = ConfDoubleOr(mass_key, 0.)*CLHEP::GeV;
    if(configured > 0.)
      return configured;
    const auto name = ConfStringOr(particle_key, fallback_particle);
    return FindParticleOrDie(name)->GetPDGMass();
  }

  G4double
  TwoBodyMomentum(G4double parent_mass,
                  G4double child1_mass,
                  G4double child2_mass)
  {
    const G4double term1 = parent_mass*parent_mass
      - (child1_mass + child2_mass)*(child1_mass + child2_mass);
    const G4double term2 = parent_mass*parent_mass
      - (child1_mass - child2_mass)*(child1_mass - child2_mass);
    if(term1 <= 0. || term2 <= 0.)
      return 0.;
    return 0.5*std::sqrt(term1*term2)/parent_mass;
  }

  G4double
  SampleFlatteCuspMass(G4double mass_min,
                       G4double mass_max,
                       G4double scattering_length_real,
                       G4double scattering_length_imag,
                       G4double momentum_transfer_mev)
  {
    return K18ReactionModelSampler::SampleFlatteCuspMass(
      mass_min, mass_max, scattering_length_real,
      scattering_length_imag, momentum_transfer_mev);
  }

  G4ThreeVector
  DirectionWithCosTheta(const G4ThreeVector& axis,
                        G4double cos_theta,
                        G4double phi)
  {
    const G4ThreeVector z_axis = axis.unit();
    const G4ThreeVector ref = std::abs(z_axis.z()) < 0.9
      ? G4ThreeVector(0., 0., 1.) : G4ThreeVector(0., 1., 0.);
    const G4ThreeVector x_axis = (ref.cross(z_axis)).unit();
    const G4ThreeVector y_axis = (z_axis.cross(x_axis)).unit();
    const G4double sin_theta =
      std::sqrt(std::max(0., 1. - cos_theta*cos_theta));
    return (cos_theta*z_axis
            + sin_theta*std::cos(phi)*x_axis
            + sin_theta*std::sin(phi)*y_axis).unit();
  }

  G4bool
  SolveTwoBodyLabMomentum(const G4LorentzVector& beam_p4,
                          const G4LorentzVector& target_p4,
                          G4double scattered_mass,
                          G4double recoil_mass,
                          const G4ThreeVector& scattered_dir,
                          G4double scan_pmax,
                          G4int scan_bins,
                          G4bool high_solution,
                          G4double& scattered_p)
  {
    const G4LorentzVector total_p4 = beam_p4 + target_p4;
    const G4ThreeVector dir = scattered_dir.unit();
    auto balance = [&](G4double p)
    {
      const G4double e = std::sqrt(p*p + scattered_mass*scattered_mass);
      const G4LorentzVector scat(dir*p, e);
      const G4LorentzVector recoil = total_p4 - scat;
      return recoil.mag2() - recoil_mass*recoil_mass;
    };

    scan_pmax = std::max(scan_pmax, total_p4.vect().mag() + 3.*CLHEP::GeV);
    scan_bins = std::max(100, scan_bins);
    const G4double root_eps = 1.e-6*CLHEP::GeV*CLHEP::GeV;
    std::vector<G4double> roots;
    auto push_root = [&](G4double root)
    {
      if(roots.empty() || std::abs(roots.back()-root) > 1.e-6*CLHEP::GeV)
        roots.push_back(root);
    };
    G4double prev_p = 0.;
    G4double prev_f = balance(prev_p);
    for(G4int i=1; i<=scan_bins; ++i){
      const G4double p = scan_pmax*static_cast<G4double>(i)/scan_bins;
      const G4double f = balance(p);
      if(std::abs(prev_f) < root_eps){
        push_root(prev_p);
      } else if(prev_f*f < 0.){
        G4double lo = prev_p;
        G4double hi = p;
        G4double flo = prev_f;
        for(G4int iter=0; iter<80; ++iter){
          const G4double mid = 0.5*(lo + hi);
          const G4double fm = balance(mid);
          if(flo*fm <= 0.){
            hi = mid;
          } else {
            lo = mid;
            flo = fm;
          }
        }
        push_root(0.5*(lo + hi));
      }
      prev_p = p;
      prev_f = f;
    }
    if(std::abs(prev_f) < root_eps)
      push_root(prev_p);
    if(roots.empty())
      return false;
    scattered_p = high_solution ? roots.back() : roots.front();
    return scattered_p > 0.;
  }

  struct TargetSample
  {
    G4LorentzVector target;
    G4LorentzVector residual;
  };

  TargetSample
  SampleTarget(G4double target_mass,
               G4bool use_bound_state,
               G4bool use_fermi,
               const G4String& fermi_model,
               G4double fermi_sigma,
               G4double target_binding,
               G4double target_nucleus_mass,
               G4double residual_mass)
  {
    if(use_bound_state)
      return {G4LorentzVector(0., 0., 0., target_nucleus_mass),
              G4LorentzVector()};
    if(!use_fermi)
      return {G4LorentzVector(0., 0., 0., target_mass),
              G4LorentzVector()};

    const G4ThreeVector fermi_gev =
      K18ReactionModelSampler::SampleFermiMomentum(
        fermi_model, fermi_sigma/CLHEP::GeV);
    const G4ThreeVector fermi = fermi_gev*CLHEP::GeV;
    const G4bool residual_spectator =
      target_nucleus_mass > 0. && residual_mass > 0.;
    if(residual_spectator){
      const G4ThreeVector residual_mom = -fermi;
      const G4double residual_e =
        std::sqrt(residual_mass*residual_mass + residual_mom.mag2());
      const G4double target_e =
        target_nucleus_mass - residual_e - target_binding;
      if(target_e <= 0.)
        Fail("invalid off-shell target energy in residual-spectator model");
      return {G4LorentzVector(fermi, target_e),
              G4LorentzVector(residual_mom, residual_e)};
    }
    const G4double target_e =
      std::sqrt(target_mass*target_mass + fermi.mag2()) - target_binding;
    if(target_e <= 0.)
      Fail("invalid target energy after ReactionTargetBindingEnergy");
    return {G4LorentzVector(fermi, target_e), G4LorentzVector()};
  }

  G4bool
  InPhiAcceptance(const G4ThreeVector& direction,
                  const G4ThreeVector& beam_direction,
                  G4double phi_min,
                  G4double phi_max)
  {
    if(phi_max - phi_min >= 2.*CLHEP::pi)
      return true;
    const G4ThreeVector z_axis = beam_direction.unit();
    const G4ThreeVector ref = std::abs(z_axis.z()) < 0.9
      ? G4ThreeVector(0., 0., 1.) : G4ThreeVector(0., 1., 0.);
    const G4ThreeVector x_axis = (ref.cross(z_axis)).unit();
    const G4ThreeVector y_axis = (z_axis.cross(x_axis)).unit();
    G4double phi = std::atan2(direction.dot(y_axis), direction.dot(x_axis));
    if(phi < 0.) phi += 2.*CLHEP::pi;
    if(phi_min <= phi_max)
      return phi >= phi_min && phi <= phi_max;
    return phi >= phi_min || phi <= phi_max;
  }
}

namespace K18MissingMassReaction
{
void
Arm(const G4Event* event,
    const G4ThreeVector& reaction_vertex,
    const G4LorentzVector& beam_source_p4,
    const G4LorentzVector& beam_source_vertex,
    const G4LorentzVector& beam_pretarget_p4)
{
  if(!event)
    Fail("cannot arm without a G4Event");
  pending.armed = true;
  pending.event_id = event->GetEventID();
  pending.approach_steps = 0;
  pending.vertex_z = reaction_vertex.z();
  pending.reaction_vertex = reaction_vertex;
  pending.beam_source_p4 = beam_source_p4;
  pending.beam_source_vertex = beam_source_vertex;
  pending.beam_pretarget_p4 = beam_pretarget_p4;
}

G4bool
IsArmedFor(const G4Track& track)
{
  if(!pending.armed || track.GetParentID() != 0)
    return false;
  const auto* event = G4EventManager::GetEventManager()->GetConstCurrentEvent();
  if(!event || event->GetEventID() != pending.event_id)
    return false;
  const G4String beam_name = ConfStringOr("ReactionBeamParticle", "kaon-");
  return track.GetDefinition()->GetParticleName() == beam_name;
}

G4double
DistanceToVertex(const G4Track& track)
{
  if(!IsArmedFor(track))
    return std::numeric_limits<G4double>::max();
  const G4double dz = pending.vertex_z - track.GetPosition().z();
  const G4double tolerance =
    ConfDoubleOr("ReactionVertexPositionTolerance", 1.e-3)*mm;
  if(dz < -tolerance)
    Fail("live beam passed the armed reaction-z plane without firing");
  if(std::abs(dz) <= tolerance)
    return 0.;
  const G4double dir_z = track.GetMomentumDirection().z();
  if(dir_z <= 0.)
    return std::numeric_limits<G4double>::max();
  // dz/dir_z is only a straight-line estimate.  In the K18 magnetic field a
  // single such step can overshoot the requested z plane by O(10 um).  Stop
  // short, let Geant4 update the live direction, and converge over a few
  // no-op process calls instead of weakening the vertex-position tolerance.
  const G4double approach_fraction =
    ConfDoubleOr("ReactionVertexApproachFraction", 0.9);
  if(!std::isfinite(approach_fraction)
     || !(approach_fraction > 0.) || !(approach_fraction < 1.))
    Fail("ReactionVertexApproachFraction must be in (0,1)");
  return approach_fraction*dz/dir_z;
}

G4bool
GenerateAtVertex(const G4Track& track, G4ParticleChange& change)
{
  if(!IsArmedFor(track))
    return false;
  const G4double position_tolerance =
    ConfDoubleOr("ReactionVertexPositionTolerance", 1.e-3)*mm;
  const G4double containment_tolerance =
    ConfDoubleOr("ReactionTargetContainmentTolerance", 0.)*mm;
  if(!std::isfinite(containment_tolerance) || containment_tolerance < 0.)
    Fail("ReactionTargetContainmentTolerance must be finite and non-negative");
  const G4double dz = pending.vertex_z - track.GetPosition().z();
  if(dz < -position_tolerance){
    Fail("forced process did not stop at the armed reaction-z plane; dz="
         + std::to_string((track.GetPosition().z()-pending.vertex_z)/mm)
         + " mm");
  }
  if(dz > position_tolerance){
    ++pending.approach_steps;
    const G4int max_approach_steps = std::max<G4int>(
      1, std::lround(ConfDoubleOr(
        "ReactionMaxVertexApproachSteps", 16.)));
    if(pending.approach_steps > max_approach_steps)
      Fail("live beam did not converge to the armed reaction-z plane after "
           + std::to_string(max_approach_steps) + " approach steps; dz="
           + std::to_string(dz/mm) + " mm");
    return false;
  }

  // Transport consumes a geometry-dependent number of random values.  Reset
  // only when the live beam has reached the reaction plane so that the target
  // and material variants of a paired study draw identical reaction
  // kinematics for the same global event number.  Without the optional key,
  // sampling continues from the transport random stream.
  ApplyReactionKinematicsSeed();

  const auto* beam_particle = track.GetDefinition();
  auto* scattered_particle = FindParticleOrDie(
    ConfStringOr("ReactionScatteredParticle", "kaon+"));
  const G4String target_name =
    ConfStringOr("ReactionTargetParticle", "proton");
  const G4String recoil_name =
    ConfStringOr("ReactionRecoilParticle", "xi-");
  auto* target_particle = particleTable->FindParticle(target_name);
  auto* recoil_particle = particleTable->FindParticle(recoil_name);
  const G4double recoil_mass_config =
    ConfDoubleOr("ReactionRecoilMass", 0.)*CLHEP::GeV;
  if(!recoil_particle && recoil_mass_config <= 0.)
    Fail("unknown recoil particle and ReactionRecoilMass is not set");

  const G4double scattered_mass = scattered_particle->GetPDGMass();
  const G4double target_mass = ParticleMassOrConfig(
    "ReactionTargetParticle", "ReactionTargetMass", target_name);
  const G4double recoil_mass = recoil_mass_config > 0.
    ? recoil_mass_config : recoil_particle->GetPDGMass();
  const G4int target_pdg_config = static_cast<G4int>(std::lround(
    ConfDoubleOr("ReactionTargetPDG", 0.)));
  const G4int target_pdg = target_pdg_config != 0 ? target_pdg_config
    : (target_particle ? target_particle->GetPDGEncoding() : 0);
  const G4int recoil_pdg = recoil_particle
    ? recoil_particle->GetPDGEncoding() : 0;

  const G4ThreeVector beam_momentum = track.GetMomentum();
  const G4LorentzVector beam_vertex_p4(
    beam_momentum, track.GetTotalEnergy());
  const G4ThreeVector reaction_vertex = track.GetPosition();
  const G4LorentzVector reaction_vtx(
    reaction_vertex, track.GetGlobalTime());
  const G4ThreeVector beam_direction = beam_momentum.unit();
  // A generated reaction is valid only when the transported beam reaches the
  // configured target.  This check is part of the generator contract.
  const auto center = geomMan.GetGlobalPosition("Target")*mm;
  const auto relative = reaction_vertex - center;
  const G4String expected_volume =
    ConfStringOr("ReactionTargetPhysicalVolume", "");
  const G4String containment_mode = ConfStringOr(
    "ReactionTargetContainmentMode",
    expected_volume.empty() ? "box" : "physical_volume");
  G4bool outside = false;
  if(containment_mode == "physical_volume"){
    if(expected_volume.empty())
      Fail("ReactionTargetPhysicalVolume is required for physical_volume mode");
    const auto* live_volume = track.GetVolume();
    outside = !live_volume || live_volume->GetName() != expected_volume;
  } else if(containment_mode == "box"){
    const G4double half_x = 0.5*ConfDoubleOr("TargetSizeX", 0.)*mm;
    const G4double half_y = 0.5*ConfDoubleOr("TargetSizeY", 0.)*mm;
    const G4double half_z = 0.5*ConfDoubleOr("TargetSizeZ", 0.)*mm;
    if(half_x <= 0. || half_y <= 0. || half_z <= 0.)
      Fail("TargetSizeX/Y/Z must be positive for box containment mode");
    outside = std::abs(relative.x()) > half_x + containment_tolerance
      || std::abs(relative.y()) > half_y + containment_tolerance
      || std::abs(relative.z()) > half_z + containment_tolerance;
  } else {
    Fail("unknown ReactionTargetContainmentMode=" + containment_mode);
  }
  if(outside){
    if(pending.event_id == 0)
      G4cout << "[K18MissingMassReaction] beam missed target; mode="
             << containment_mode << " volume="
             << (track.GetVolume() ? track.GetVolume()->GetName() : "none")
             << " relative_position=" << relative/mm << " mm" << G4endl;
    change.SetNumberOfSecondaries(0);
    change.ProposeTrackStatus(fStopAndKill);
    change.ProposeEnergy(0.);
    pending.armed = false;
    return true;
  }

  const G4String fermi_model =
    ConfStringOr("ReactionFermiModel", "none");
  const G4bool fermi_enabled =
    fermi_model != "none" && fermi_model != "None"
    && fermi_model != "NONE" && fermi_model != "0"
    && fermi_model != "false" && fermi_model != "False"
    && fermi_model != "FALSE";
  const G4bool use_fermi =
    ConfBoolOr("ReactionUseFermiMotion", fermi_enabled);
  const G4double fermi_sigma =
    ConfDoubleOr("ReactionFermiSigmaP", 0.090)*CLHEP::GeV;
  const G4double target_binding = ConfDoubleOr(
    "ReactionTargetBindingEnergy",
    ConfDoubleOr("ReactionSeparationEnergy", 0.))*CLHEP::GeV;
  const G4double target_nucleus_mass =
    ConfDoubleOr("ReactionTargetNucleusMass", 0.)*CLHEP::GeV;
  const G4double residual_mass =
    ConfDoubleOr("ReactionResidualMass", 0.)*CLHEP::GeV;
  const G4double bound_state_mass =
    ConfDoubleOr("ReactionBoundStateMass", 0.)*CLHEP::GeV;
  const G4double bound_state_dm = ConfDoubleOr(
    "ReactionBoundDMMeV",
    ConfDoubleOr("ReactionBoundExcessEnergyMeV", 0.))*CLHEP::MeV;
  const G4bool use_bound_state =
    ConfBoolOr("ReactionUseBoundState", bound_state_mass > 0.);

  G4double reaction_recoil_mass = recoil_mass;
  if(use_bound_state){
    if(target_nucleus_mass <= 0.)
      Fail("ReactionUseBoundState requires ReactionTargetNucleusMass");
    if(bound_state_mass > 0.){
      reaction_recoil_mass = bound_state_mass;
    } else {
      if(residual_mass <= 0.)
        Fail("bound state requires ReactionResidualMass or mass override");
      reaction_recoil_mass = residual_mass + recoil_mass + bound_state_dm;
    }
  }

  const G4String recoil_model =
    ConfStringOr("ReactionRecoilMassModel", "fixed");
  if(recoil_model == "flatte_cusp"){
    if(use_bound_state)
      Fail("Flatte cusp recoil model cannot be combined with bound state");
    const G4double threshold = FindParticleOrDie("sigma+")->GetPDGMass()
      + FindParticleOrDie("neutron")->GetPDGMass();
    const G4double mass_min =
      ConfDoubleOr("ReactionRecoilMassMin",
                   threshold/CLHEP::GeV - 0.100)*CLHEP::GeV;
    const G4double mass_max =
      ConfDoubleOr("ReactionRecoilMassMax",
                   threshold/CLHEP::GeV + 0.100)*CLHEP::GeV;
    reaction_recoil_mass = SampleFlatteCuspMass(
      mass_min, mass_max,
      ConfDoubleOr("ReactionCuspScatteringLengthReal", 2.06),
      ConfDoubleOr("ReactionCuspScatteringLengthImag", 4.64),
      ConfDoubleOr("ReactionCuspMomentumTransferMeV", 186.));
  } else if(recoil_model != "fixed"){
    Fail("unknown ReactionRecoilMassModel=" + recoil_model);
  }

  const G4double theta_min =
    ConfDoubleOr("ReactionThetaMin", 0.)*CLHEP::deg;
  const G4double theta_max =
    ConfDoubleOr("ReactionThetaMax", 16.)*CLHEP::deg;
  const G4double phi_min =
    ConfDoubleOr("ReactionPhiMin", 0.)*CLHEP::deg;
  const G4double phi_max =
    ConfDoubleOr("ReactionPhiMax", 360.)*CLHEP::deg;
  const G4bool flat_costheta =
    ConfBoolOr("ReactionFlatCosTheta", false);
  const G4bool sample_cm =
    ConfBoolOr("ReactionSampleCM", use_fermi);
  const G4bool high_solution =
    ConfBoolOr("ReactionHighMomentumSolution", true);
  const G4double solve_pmax =
    ConfDoubleOr("ReactionSolvePMax", 5.)*CLHEP::GeV;
  const G4int solve_bins = std::max<G4int>(
    100, std::lround(ConfDoubleOr("ReactionSolveBins", 1000.)));
  const G4int max_trials = std::max<G4int>(
    1, std::lround(ConfDoubleOr("ReactionMaxTrial", 1000.)));

  G4LorentzVector target_p4;
  G4LorentzVector residual_p4;
  G4LorentzVector scattered_p4;
  G4bool accepted = false;
  for(G4int trial=0; trial<max_trials; ++trial){
    const auto target = SampleTarget(
      target_mass, use_bound_state, use_fermi, fermi_model,
      fermi_sigma, target_binding, target_nucleus_mass, residual_mass);
    target_p4 = target.target;
    residual_p4 = target.residual;
    if(sample_cm){
      const G4LorentzVector total = beam_vertex_p4 + target_p4;
      const G4double sqrt_s = total.mag();
      const G4double p_cm = TwoBodyMomentum(
        sqrt_s, scattered_mass, reaction_recoil_mass);
      if(!(p_cm > 0.))
        continue;
      const G4ThreeVector cm_dir = G4RandomDirection();
      scattered_p4 = G4LorentzVector(
        cm_dir*p_cm,
        std::sqrt(p_cm*p_cm + scattered_mass*scattered_mass));
      scattered_p4.boost(total.vect()/total.e());
      const G4ThreeVector direction = scattered_p4.vect().unit();
      const G4double theta = std::acos(std::max(
        -1., std::min(1., beam_direction.dot(direction))));
      if(theta < theta_min || theta > theta_max
         || !InPhiAcceptance(direction, beam_direction, phi_min, phi_max))
        continue;
      accepted = true;
      break;
    }

    G4double theta;
    if(flat_costheta){
      theta = std::acos(G4RandFlat::shoot(
        std::cos(theta_max), std::cos(theta_min)));
    } else {
      theta = G4RandFlat::shoot(theta_min, theta_max);
    }
    const G4double phi = G4RandFlat::shoot(phi_min, phi_max);
    const G4ThreeVector direction =
      DirectionWithCosTheta(beam_direction, std::cos(theta), phi);
    G4double momentum = 0.;
    if(!SolveTwoBodyLabMomentum(
         beam_vertex_p4, target_p4, scattered_mass,
         reaction_recoil_mass, direction,
         solve_pmax, solve_bins, high_solution, momentum))
      continue;
    scattered_p4 = G4LorentzVector(
      direction*momentum,
      std::sqrt(momentum*momentum + scattered_mass*scattered_mass));
    accepted = true;
    break;
  }
  if(!accepted)
    Fail("no two-body solution at the live transported vertex after "
         + std::to_string(max_trials) + " trials");

  const G4LorentzVector recoil_p4 =
    beam_vertex_p4 + target_p4 - scattered_p4;
  auto* secondary = new G4DynamicParticle(
    scattered_particle, scattered_p4.vect().unit(),
    scattered_p4.e() - scattered_mass);
  change.SetNumberOfSecondaries(1);
  change.AddSecondary(secondary);
  change.ProposeTrackStatus(fStopAndKill);
  change.ProposeEnergy(0.);

  auto& ana = S2SAnaManager::GetInstance();
  const G4int none = static_cast<G4int>(GenBranch::ParticleId::None);
  const G4int scattered_pdg = scattered_particle->GetPDGEncoding();
  ana.SetGeneratedParticle(GenBranch::kPrimPi, none, scattered_pdg,
                           scattered_p4, reaction_vtx);
  ana.SetGeneratedParticle(GenBranch::kReactionBeam, none,
                           beam_particle->GetPDGEncoding(),
                           pending.beam_pretarget_p4,
                           G4LorentzVector(pending.reaction_vertex, 0.));
  // Store the live transported state twice under explicit contracts.  The
  // analyzer compares this transport audit branch to ReactionBeamVertex and
  // aborts if later changes ever substitute a different beam state.
  ana.SetGeneratedParticle(GenBranch::kReactionBeamVertexTransport, none,
                           beam_particle->GetPDGEncoding(),
                           beam_vertex_p4, reaction_vtx);
  ana.SetGeneratedParticle(GenBranch::kReactionBeamVertex, none,
                           beam_particle->GetPDGEncoding(),
                           beam_vertex_p4, reaction_vtx);
  ana.SetGeneratedParticle(GenBranch::kReactionTarget, none, target_pdg,
                           target_p4, reaction_vtx);
  ana.SetGeneratedParticle(GenBranch::kReactionResidual, none, 0,
                           residual_p4, reaction_vtx);
  ana.SetGeneratedParticle(GenBranch::kReactionScat, none, scattered_pdg,
                           scattered_p4, reaction_vtx);
  ana.SetGeneratedParticle(GenBranch::kReactionRecoil, none,
                           use_bound_state ? 0 : recoil_pdg,
                           recoil_p4, reaction_vtx);

  const auto target_center = geomMan.GetGlobalPosition("Target")*mm;
  const auto target_relative = reaction_vertex - target_center;
  const auto scattered_dir = scattered_p4.vect().unit();
  const G4double u = std::abs(scattered_dir.z()) > 1.e-12
    ? scattered_dir.x()/scattered_dir.z() : 0.;
  const G4double v = std::abs(scattered_dir.z()) > 1.e-12
    ? scattered_dir.y()/scattered_dir.z() : 0.;
  ana.SetPrimaryParticle(0, scattered_pdg, scattered_p4, reaction_vtx);
  ana.SetPrimaryParticle(1, beam_particle->GetPDGEncoding(),
                         pending.beam_source_p4,
                         pending.beam_source_vertex);
  ana.SetPrimaryData(target_relative.x(), target_relative.y(),
                     target_relative.z(), u, v,
                     scattered_dir.phi(), scattered_dir.theta(),
                     scattered_p4.vect().mag(),
                     pending.beam_pretarget_p4.vect().mag(), scattered_pdg);

  if(pending.event_id == 0){
    G4cout << "[K18MissingMassReaction] generated at transported vertex; z="
           << reaction_vertex.z()/mm << " mm, p_beam="
           << beam_vertex_p4.vect().mag()/CLHEP::GeV << " GeV/c"
           << G4endl;
  }
  pending.armed = false;
  return true;
}
}

K18MissingMassReactionProcess::K18MissingMassReactionProcess()
  : G4VDiscreteProcess("K18MissingMassReaction", fUserDefined)
{
}

G4bool
K18MissingMassReactionProcess::IsApplicable(const G4ParticleDefinition& particle)
{
  return particle.GetParticleName()
    == ConfStringOr("ReactionBeamParticle", "kaon-");
}

G4double
K18MissingMassReactionProcess::PostStepGetPhysicalInteractionLength(
  const G4Track& track, G4double, G4ForceCondition* condition)
{
  *condition = NotForced;
  return K18MissingMassReaction::DistanceToVertex(track);
}

G4VParticleChange*
K18MissingMassReactionProcess::PostStepDoIt(const G4Track& track, const G4Step&)
{
  m_particle_change.Initialize(track);
  if(!K18MissingMassReaction::GenerateAtVertex(track, m_particle_change))
    return &m_particle_change;
  return &m_particle_change;
}

G4double
K18MissingMassReactionProcess::GetMeanFreePath(
  const G4Track&, G4double, G4ForceCondition* condition)
{
  *condition = NotForced;
  return std::numeric_limits<G4double>::max();
}
