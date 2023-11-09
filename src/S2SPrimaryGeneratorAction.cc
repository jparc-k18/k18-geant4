// -*- C++ -*-

#include "S2SPrimaryGeneratorAction.hh"

#include <G4Event.hh>
#include <G4ParticleGun.hh>
#include <G4ParticleTable.hh>
#include <G4ParticleDefinition.hh>
#include <G4LorentzVector.hh>
#include <G4ThreeVector.hh>
#include <Randomize.hh>

#include <TMath.h>

#include "ConfMan.hh"
#include "DCGeomMan.hh"
#include "FuncName.hh"
#include "S2SAnaManager.hh"
#include "DetSizeMan.hh"

namespace
{
using CLHEP::mm;
using CLHEP::deg;
using CLHEP::radian;
using CLHEP::GeV;
const auto& confMan = ConfMan::GetInstance();
const auto& geomMan = DCGeomMan::GetInstance();
const auto& sizeMan = DetSizeMan::GetInstance();
auto& anaMan = S2SAnaManager::GetInstance();
const auto particleTable = G4ParticleTable::GetParticleTable();
}

//_____________________________________________________________________________
S2SPrimaryGeneratorAction::S2SPrimaryGeneratorAction()
  : G4VUserPrimaryGeneratorAction(),
    m_generator(confMan.Get<G4int>("Generator"))
{
}

//_____________________________________________________________________________
S2SPrimaryGeneratorAction::~S2SPrimaryGeneratorAction()
{
  if(m_particleGun) delete m_particleGun;
}

//_____________________________________________________________________________
void
S2SPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  if(m_particleGun) delete m_particleGun;
  switch(m_generator){
  case 0: GenerateDemo(anEvent); break;
  case 1: GenerateMonochromeBeam(anEvent); break;
  case 2: GenerateUniformSpherical(anEvent); break;
  default:
    G4cerr << " * Generator number error : " << m_generator << G4endl;
    break;
  }
}

//_____________________________________________________________________________
void
S2SPrimaryGeneratorAction::GenerateDemo(G4Event* anEvent)
{
  static const G4int n_particle = 1;
  m_particleGun = new G4ParticleGun(n_particle);
  static const G4String name = "kaon+";
  // static const G4String name = "proton";
  static const auto particle = particleTable->FindParticle(name);
  static const auto pdg = particle->GetPDGEncoding();
  static const auto& target_pos = geomMan.GetGlobalPosition("Target")*mm;
  static const auto& target_size = sizeMan.GetSize("Target")*mm/2;
  static const G4double m0 = particle->GetPDGMass();
  static const G4int experiment = confMan.Get<G4int>("Experiment");
  const auto evnum = anEvent->GetEventID();
  G4double p0 = (experiment == 10) ? 0.9*GeV : 1.3*GeV;
  p0 = (evnum%3 == 0) ? p0
    : (evnum%3 == 1) ? p0*1.075
    : p0*0.925;
  auto n = G4RandFlat::shootInt(5);
  G4double theta = (n%5 == 0) ? 0*deg
    : (n%5 == 1) ? 2*deg
    : (n%5 == 2) ? 4*deg
    : (n%5 == 3) ? -2*deg
    : -4*deg;
  G4double phi = 0*deg; // G4RandFlat::shoot(0., 360.)*deg;
  G4LorentzVector p(0, 0, 0, TMath::Sqrt(p0*p0 + m0*m0));
  p.setRThetaPhi(p0, theta, phi);
  G4double x0 =  G4RandFlat::shoot(-target_size.x(), target_size.x());
  G4double y0 =  G4RandFlat::shoot(-target_size.y(), target_size.y());
  G4double z0 =  G4RandFlat::shoot(-target_size.z(), target_size.z());
  G4LorentzVector v(target_pos, 0);
  // G4LorentzVector v(target_pos + G4ThreeVector(0, 0, 300*CLHEP::mm), 0);
#if 0
  G4cout << FUNC_NAME << G4endl
         << " " << p0 << " " << theta/deg << " " << phi/deg << G4endl
         << " " << p << " " << p.theta()/deg << " " << v << G4endl;
#endif
  m_particleGun->SetParticleDefinition(particle);
  m_particleGun->SetParticleMomentumDirection(p.v());
  m_particleGun->SetParticleEnergy(p.e() - m0);
  m_particleGun->SetParticlePosition(v.v());
  m_particleGun->GeneratePrimaryVertex(anEvent);
  anaMan.SetPrimaryParticle(0, pdg, p, v);
}

//_____________________________________________________________________________
void
S2SPrimaryGeneratorAction::GenerateMonochromeBeam(G4Event* anEvent)
{
  const G4int n_particle = 1;
  m_particleGun = new G4ParticleGun(n_particle);
  static const G4String name = "kaon+";
  // static const G4String name = "proton";
  static const auto particle = particleTable->FindParticle(name);
  static const auto pdg = particle->GetPDGEncoding();
  const G4double m0 = particle->GetPDGMass();
  // const G4double p0 = 1.37*GeV;
  const G4double p0 = 0.7*GeV;
  const auto& target_pos = geomMan.GetGlobalPosition("Target");
  G4LorentzVector p(0, 0, p0, TMath::Sqrt(p0*p0 + m0*m0));
  G4LorentzVector v(target_pos, 0);
  m_particleGun->SetParticleDefinition(particle);
  m_particleGun->SetParticleMomentumDirection(p.v());
  m_particleGun->SetParticleEnergy(p.e() - m0);
  m_particleGun->SetParticlePosition(v.v());
  m_particleGun->GeneratePrimaryVertex(anEvent);
  anaMan.SetPrimaryParticle(0, pdg, p, v);
}

//_____________________________________________________________________________
void
S2SPrimaryGeneratorAction::GenerateUniformSpherical(G4Event* anEvent)
{
  static const G4int n_particle = 1;
  m_particleGun = new G4ParticleGun(n_particle);
  static const G4String name = "kaon+";
  // static const G4String name = "proton";
  static const auto particle = particleTable->FindParticle(name);
  static const auto pdg = particle->GetPDGEncoding();
  static const auto& target_pos = geomMan.GetGlobalPosition("Target")*mm;
  static const auto& target_size = sizeMan.GetSize("Target")*mm/2;
  static const G4double m0 = particle->GetPDGMass();
  static const G4int experiment = confMan.Get<G4int>("Experiment");
  G4double p0 = (experiment == 10)
    ? G4RandFlat::shoot(0.6, 1.2)*GeV
    : G4RandFlat::shoot(1.0, 1.8)*GeV;
  G4double theta =
    std::acos(G4RandFlat::shoot(std::cos(0*deg), std::cos(20*deg)))*radian;
  G4double phi = G4RandFlat::shoot(0., 360.)*deg;
  G4LorentzVector p(0, 0, 0, TMath::Sqrt(p0*p0 + m0*m0));
  p.setRThetaPhi(p0, theta, phi);
  G4double x0 =  G4RandFlat::shoot(-target_size.x(), target_size.x());
  G4double y0 =  G4RandFlat::shoot(-target_size.y(), target_size.y());
  G4double z0 =  G4RandFlat::shoot(-target_size.z(), target_size.z());
  G4LorentzVector v(target_pos, 0);
  // G4LorentzVector v(target_pos + G4ThreeVector(x0, y0, z0), 0);
#if 0
  G4cout << FUNC_NAME << G4endl
         << " " << p0 << " " << theta/deg << " " << phi/deg << G4endl
         << " " << p << " " << p.theta()/deg << " " << v << G4endl;
#endif
  m_particleGun->SetParticleDefinition(particle);
  m_particleGun->SetParticleMomentumDirection(p.v());
  m_particleGun->SetParticleEnergy(p.e() - m0);
  m_particleGun->SetParticlePosition(v.v());
  m_particleGun->GeneratePrimaryVertex(anEvent);
  anaMan.SetPrimaryParticle(0, pdg, p, v);
}
