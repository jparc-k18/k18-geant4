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
#include <TFile.h>
#include <TTree.h>
#include <TGraph.h>

#include "BeamMan.hh"
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
const auto& beamMan = BeamMan::GetInstance();
const auto& confMan = ConfMan::GetInstance();
const auto& geomMan = DCGeomMan::GetInstance();
const auto& sizeMan = DetSizeMan::GetInstance();
auto& anaMan = S2SAnaManager::GetInstance();
const auto particleTable = G4ParticleTable::GetParticleTable();
const auto& zK18Target = geomMan.LocalZ("K18Target");
BeamInfo beam;
}

//_____________________________________________________________________________
S2SPrimaryGeneratorAction::S2SPrimaryGeneratorAction()
  : G4VUserPrimaryGeneratorAction(),
    m_particleGun(nullptr),
    m_generator(confMan.Get<G4int>("Generator"))
{
  auto igene = confMan.Get<G4int>("Generator");  
  if(igene==4 ||igene==7501 ){
    auto ifsK18name = confMan.Get<G4String>("K18ROOT");  
    profileK18 = new TFile(ifsK18name);
    if(!profileK18){G4cout<< ifsK18name << " is not found." <<G4endl;}
    k18track = (TTree*)profileK18 ->Get("k18track");
    k18track ->SetBranchAddress("p_3rd",&p_3rd);
    k18track ->SetBranchAddress("xtgtK18",&xtgt);
    k18track ->SetBranchAddress("ytgtK18",&ytgt);
    k18track ->SetBranchAddress("utgtK18",&utgt);
    k18track ->SetBranchAddress("vtgtK18",&vtgt);
    k18track ->SetBranchAddress("trigflag",&trigflag);
    k18track ->SetBranchAddress("CBtof0",&CBtof0);
    k18track ->SetBranchAddress("ntK18",&ntK18);
    k18track ->SetBranchAddress("chisqrK18",&chisqrK18);
  }

  if(igene==7501 ){
    auto ifsTheoname = confMan.Get<G4String>("TheoCalc");  
    gr = new TGraph(ifsTheoname, "%lg %*lg %*lg %*lg %lg");
    for(Int_t i=0;i<gr->GetN();++i){
      G4double x,y;
      gr-> GetPoint(i,x,y);
      Sum += y;
    }
  }
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
  beam = beamMan.Get();

  if(m_particleGun) delete m_particleGun;
  switch(m_generator){
  case 0: GenerateDemo(anEvent); break;
  case 1: GenerateMonochromeBeam(anEvent); break;
  case 2: GenerateUniformSpherical(anEvent); break;
  case 3: GenerateBeam(anEvent); break;
  case 4: GenerateBeamThrough(anEvent); break;
  case 7001: Generate12XiBe(anEvent); break;
  case 7002: GenerateElementaryXiMinus(anEvent); break;
  case 7003: GenerateElementarySigmaMinus(anEvent); break;
  case 7004: GenerateElementarySigmaPlus(anEvent); break;
  //case 7005: Generate12XiBePeakStructure(anEvent); break;
  case 7501: GenerateKH7XiHSpectrum(anEvent); break;
  default:
    G4cerr << " * Generator number error : " << m_generator << G4endl;
    break;
  }
}

//_____________________________________________________________________________
void // 0
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
  G4double p0 = (experiment == 10) ? 0.9*GeV : 1.4*GeV;
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
  // G4double x0 =  G4RandFlat::shoot(-target_size.x(), target_size.x());
  // G4double y0 =  G4RandFlat::shoot(-target_size.y(), target_size.y());
  // G4double z0 =  G4RandFlat::shoot(-target_size.z(), target_size.z());
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
void // 1
S2SPrimaryGeneratorAction::GenerateMonochromeBeam(G4Event* anEvent)
{
  const G4int n_particle = 1;
  m_particleGun = new G4ParticleGun(n_particle);
  static const G4String name = "kaon+";
  // static const G4String name = "proton";
  static const auto particle = particleTable->FindParticle(name);
  static const auto pdg = particle->GetPDGEncoding();
  static const G4int experiment = confMan.Get<G4int>("Experiment");
  const G4double m0 = particle->GetPDGMass();
  G4double p0 = (experiment == 10) ? 0.9*GeV : 1.4*GeV;
  const auto& target_pos = geomMan.GetGlobalPosition("Target");
  beam.pos.setX(target_pos.x());
  beam.pos.setY(target_pos.y());
  beam.pos.setZ(target_pos.z()-956.*mm);
  G4LorentzVector p(0, 0, p0, TMath::Sqrt(p0*p0 + m0*m0));
  G4LorentzVector v(beam.pos, 0);
  m_particleGun->SetParticleDefinition(particle);
  m_particleGun->SetParticleMomentumDirection(p.v());
  m_particleGun->SetParticleEnergy(p.e() - m0);
  m_particleGun->SetParticlePosition(v.v());
  m_particleGun->GeneratePrimaryVertex(anEvent);
  anaMan.SetPrimaryParticle(0, pdg, p, v);
  anaMan.SetPrimaryData(beam.pos.x(),beam.pos.y(),beam.pos.z(),
			0.,0.,0.,0.,p0,0.,9999);
}

//_____________________________________________________________________________
void // 2
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
    : G4RandFlat::shoot(1.37, 1.38)*GeV;
  G4double theta =
    std::acos(G4RandFlat::shoot(std::cos(0*deg), std::cos(20*deg)))*radian;
  G4double phi = G4RandFlat::shoot(0., 360.)*deg;
  G4LorentzVector p(0, 0, 0, TMath::Sqrt(p0*p0 + m0*m0));
  p.setRThetaPhi(p0, theta, phi);
  G4double u0, v0;
  u0 = TMath::Tan(theta)*TMath::Cos(phi);
  v0 = TMath::Tan(theta)*TMath::Sin(phi);
  beam.VO(zK18Target);
  beam.pos.setZ(0);
  G4LorentzVector v(target_pos + beam.pos, 0);
#if 0
  beam.Print();
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
  anaMan.SetPrimaryData(0,0,0,u0,v0,phi,theta,p0,1.8*GeV,9999);
}

//_____________________________________________________________________________
void // 3
S2SPrimaryGeneratorAction::GenerateBeam(G4Event* anEvent)
{
  const G4int n_particle = 1;
  m_particleGun = new G4ParticleGun(n_particle);
  static const G4int experiment = confMan.Get<G4int>("Experiment");
  G4String name;
  const G4double p0 = confMan.Get<G4double>("PK18")*CLHEP::GeV;
  if(experiment == 10){
    name = "pi-";
  }else{
    name = "kaon-";
  }
  static const auto particle = particleTable->FindParticle(name);
  static const auto pdg = particle->GetPDGEncoding();
  const G4double m0 = particle->GetPDGMass();
  const auto& target_pos = geomMan.GetGlobalPosition("Target");
  beam.VO(zK18Target);
  beam.pos.setZ(0);
  beam.mom.setMag(p0);
  G4LorentzVector p(beam.mom, TMath::Sqrt(p0*p0 + m0*m0));
  G4LorentzVector v(target_pos + beam.pos, 0);
#if 0
  beam.Print();
  G4cout << FUNC_NAME << G4endl
         << " " << p << " " << v << G4endl;
#endif
  m_particleGun->SetParticleDefinition(particle);
  m_particleGun->SetParticleMomentumDirection(p.v());
  m_particleGun->SetParticleEnergy(p.e() - m0);
  m_particleGun->SetParticlePosition(v.v());
  m_particleGun->GeneratePrimaryVertex(anEvent);
  anaMan.SetPrimaryParticle(0, pdg, p, v);
}

void // 4
S2SPrimaryGeneratorAction::GenerateBeamThrough(G4Event* anEvent)
{
  static const G4int n_particle = 1;
  m_particleGun = new G4ParticleGun(n_particle);
  static const G4String name = "kaon+";
  static const auto particle = particleTable->FindParticle(name);
  static const auto pdg = particle->GetPDGEncoding();
  static const auto& target_pos = geomMan.GetGlobalPosition("Target")*mm;
  static const auto& bac1_pos = geomMan.GetGlobalPosition("BAC1")*mm;
  static const G4double m0 = particle->GetPDGMass();
  static const G4int experiment = confMan.Get<G4int>("Experiment");
  G4double theta;
  //  std::acos(G4RandFlat::shoot(std::cos(0*deg), std::cos(20*deg)))*radian;
  //
  G4double p0, x0, y0, z0;
  G4double u0, v0;
  G4int count =0;
  while(count<1){
    G4int i = int(G4UniformRand()*(k18track->GetEntries()));
    k18track->GetEntry(i);
    if(/*trigflag[21]>0 &&*/ ntK18==1 && abs(CBtof0[0])<0.2 && chisqrK18[0]<10.){
    z0 = bac1_pos.z() - target_pos.z() -556.*mm;
    // z0 = 556 mm upstrm from BAC1 (approximately dwnstr surface of BH2)
    x0 = xtgt[0] + utgt[0]*z0; // Horizontal direction
    y0 = ytgt[0] + vtgt[0]*z0; // vertical direction
    u0 = utgt[0];
    v0 = vtgt[0];
    p0 = p_3rd[0]*GeV;
    count++;
   }
  }
  beam.VO(zK18Target);
  beam.pos.setX(x0);
  beam.pos.setY(y0);
  beam.pos.setZ(z0);
  G4LorentzVector p(0, 0, 0, TMath::Sqrt(p0*p0 + m0*m0));
  p.setX(p0*u0/TMath::Sqrt(1+u0*u0+v0*v0));
  p.setY(p0*v0/TMath::Sqrt(1+u0*u0+v0*v0));
  p.setZ(p0/TMath::Sqrt(1+u0*u0+v0*v0));
  G4LorentzVector v(beam.pos+target_pos, 0);
#if 0
  beam.Print();
  G4cout << FUNC_NAME << G4endl
         << " " << p << " " << v << G4endl;
#endif
  m_particleGun->SetParticleDefinition(particle);
  m_particleGun->SetParticleMomentumDirection(p.v());
  m_particleGun->SetParticleEnergy(p.e() - m0);
  m_particleGun->SetParticlePosition(v.v());
  m_particleGun->GeneratePrimaryVertex(anEvent);
  anaMan.SetPrimaryParticle(0, pdg, p, v);
  anaMan.SetPrimaryData(x0,y0,z0,u0,v0,0.,0.,p0,p0,9999);
}

void // 7001 E70 12C(KK)12XiBe kinematics 
S2SPrimaryGeneratorAction::Generate12XiBe(G4Event* anEvent)
{
  static const G4int n_particle = 1;
  m_particleGun = new G4ParticleGun(n_particle);
  static const G4String name = "kaon+";
  static const auto particle = particleTable->FindParticle(name);
  static const auto pdg = particle->GetPDGEncoding();
  static const auto& target_pos = geomMan.GetGlobalPosition("Target")*mm;
  static const auto& target_size = sizeMan.GetSize("Target")*mm/2;
  static const G4double m0 = particle->GetPDGMass();
  static const G4int experiment = confMan.Get<G4int>("Experiment");
  G4double theta;
  //  std::acos(G4RandFlat::shoot(std::cos(0*deg), std::cos(20*deg)))*radian;
  G4double phi = G4RandFlat::shoot(0., 360.)*deg;
  G4double pB = 1.8*GeV;

  G4double mass_12C = 11.177929*GeV; // 7Li
  G4double mass_12XiBe = (10.2551+1.32171)*GeV; // 6He + Xi
  G4double m_tgt = mass_12C; // mass of target nucleus (7Li)
  G4double m_hyp = mass_12XiBe; // mass of target nucleus (7Li)
  G4double p0;
  G4double cost = 0.92;

  while(cost>1. || cost<0.93){
    p0 = G4RandFlat::shoot(1.33, 1.39)*GeV; // scatter momentum
    G4double Energy_B = TMath::Sqrt(pB*pB + m0*m0);
    G4double Energy_S = TMath::Sqrt(p0*p0 + m0*m0);
    G4double Energy_hyp = Energy_B + m_tgt - Energy_S;
  
    cost = (m_hyp*m_hyp - Energy_hyp*Energy_hyp + pB*pB + p0*p0)/(2.*pB*p0);
  }

  theta = std::acos(cost)*radian;

  G4LorentzVector p(0, 0, 0, TMath::Sqrt(p0*p0 + m0*m0));
  p.setRThetaPhi(p0, theta, phi);
  beam.VO(zK18Target);
  beam.pos.setZ(0);
  G4LorentzVector v(target_pos + beam.pos, 0);
#if 0
  beam.Print();
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

void // 7002 E70 p(KK)Xi kinematics 
S2SPrimaryGeneratorAction::GenerateElementaryXiMinus(G4Event* anEvent)
{
  static const G4int n_particle = 1;
  m_particleGun = new G4ParticleGun(n_particle);
  static const G4String name = "kaon+";
  static const auto particle = particleTable->FindParticle(name);
  static const auto pdg = particle->GetPDGEncoding();
  static const auto& target_pos = geomMan.GetGlobalPosition("Target")*mm;
  static const auto& target_size = sizeMan.GetSize("Target")*mm/2;
  static const G4double m0 = particle->GetPDGMass();
  static const G4int experiment = confMan.Get<G4int>("Experiment");
  G4double theta;
  //  std::acos(G4RandFlat::shoot(std::cos(0*deg), std::cos(20*deg)))*radian;
  G4double phi = G4RandFlat::shoot(0., 360.)*deg;
  G4double pB = 1.8*GeV;

  G4double m_tgt = particleTable->FindParticle("proton") ->GetPDGMass(); // mass of target nucleus (p)
  G4double m_hyp = particleTable->FindParticle("xi-") ->GetPDGMass(); // mass of target nucleus (Xi)
  G4double p0;
  G4double cost = 0.92;

  while(cost>1. || cost<0.93){
    p0 = G4RandFlat::shoot(1.1, 1.32)*GeV; // scatter momentum
    G4double Energy_B = TMath::Sqrt(pB*pB + m0*m0);
    G4double Energy_S = TMath::Sqrt(p0*p0 + m0*m0);
    G4double Energy_hyp = Energy_B + m_tgt - Energy_S;
  
    cost = (m_hyp*m_hyp - Energy_hyp*Energy_hyp + pB*pB + p0*p0)/(2.*pB*p0);
  }

  theta = std::acos(cost)*radian;

  G4LorentzVector p(0, 0, 0, TMath::Sqrt(p0*p0 + m0*m0));
  p.setRThetaPhi(p0, theta, phi);
  beam.VO(zK18Target);
  beam.pos.setZ(0);
  G4LorentzVector v(target_pos + beam.pos, 0);
#if 0
  beam.Print();
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

void // 7003 E70 p(Kpi)SigmaMinus kinematics 
S2SPrimaryGeneratorAction::GenerateElementarySigmaMinus(G4Event* anEvent)
{
  static const G4int n_particle = 1;
  m_particleGun = new G4ParticleGun(n_particle);
  static const G4String name = "pi+";
  static const auto particle = particleTable->FindParticle(name);
  static const auto pdg = particle->GetPDGEncoding();
  static const auto& target_pos = geomMan.GetGlobalPosition("Target")*mm;
  static const auto& target_size = sizeMan.GetSize("Target")*mm/2;
  static const G4double m0 = particle->GetPDGMass();
  static const G4int experiment = confMan.Get<G4int>("Experiment");
  G4double theta;
  //  std::acos(G4RandFlat::shoot(std::cos(0*deg), std::cos(20*deg)))*radian;
  G4double phi = G4RandFlat::shoot(0., 360.)*deg;
  static const G4double pB = 1.8*GeV;
  static const G4double mB = particleTable->FindParticle("kaon-") ->GetPDGMass(); // mass of target nucleus (p)

  static const G4double m_tgt = particleTable->FindParticle("proton") ->GetPDGMass(); // mass of target nucleus (p)
  static const G4double m_hyp = particleTable->FindParticle("sigma-") ->GetPDGMass(); // mass of target nucleus (Sigma)
  G4double p0;
  G4double cost = 0.92;

  while(cost>1. || cost<0.93){
    p0 = G4RandFlat::shoot(1.4, 1.6)*GeV; // scatter momentum
    G4double Energy_B = TMath::Sqrt(pB*pB + mB*mB);
    G4double Energy_S = TMath::Sqrt(p0*p0 + m0*m0);
    G4double Energy_hyp = Energy_B + m_tgt - Energy_S;
  
    cost = (m_hyp*m_hyp - Energy_hyp*Energy_hyp + pB*pB + p0*p0)/(2.*pB*p0);
  }

  theta = std::acos(cost)*radian;

  G4LorentzVector p(0, 0, 0, TMath::Sqrt(p0*p0 + m0*m0));
  p.setRThetaPhi(p0, theta, phi);
  beam.VO(zK18Target);
  beam.pos.setZ(0);
  G4LorentzVector v(target_pos + beam.pos, 0);
#if 0
  beam.Print();
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

void // 7004 E70 p(pi,K)SigmaPlus kinematics 
S2SPrimaryGeneratorAction::GenerateElementarySigmaPlus(G4Event* anEvent)
{
  static const G4int n_particle = 1;
  m_particleGun = new G4ParticleGun(n_particle);
  static const G4String name = "kaon+";
  static const auto particle = particleTable->FindParticle(name);
  static const auto pdg = particle->GetPDGEncoding();
  static const auto& target_pos = geomMan.GetGlobalPosition("Target")*mm;
  static const auto& target_size = sizeMan.GetSize("Target")*mm/2;
  static const G4double m0 = particle->GetPDGMass();
  static const G4int experiment = confMan.Get<G4int>("Experiment");
  G4double theta;
  //  std::acos(G4RandFlat::shoot(std::cos(0*deg), std::cos(20*deg)))*radian;
  G4double phi = G4RandFlat::shoot(0., 360.)*deg;
  G4double pB = 1.8*GeV;
  G4double mB = particleTable->FindParticle("pi+") ->GetPDGMass(); // mass of target nucleus (p)

  G4double m_tgt = particleTable->FindParticle("proton") ->GetPDGMass(); // mass of target nucleus (p)
  G4double m_hyp = particleTable->FindParticle("sigma+") ->GetPDGMass(); // mass of target nucleus (Sigma)
  G4double p0;
  G4double cost = 0.92;

  while(cost>1. || cost<0.93){
    p0 = G4RandFlat::shoot(1.2, 1.42)*GeV; // scatter momentum
    G4double Energy_B = TMath::Sqrt(pB*pB + mB*mB);
    G4double Energy_S = TMath::Sqrt(p0*p0 + m0*m0);
    G4double Energy_hyp = Energy_B + m_tgt - Energy_S;
  
    cost = (m_hyp*m_hyp - Energy_hyp*Energy_hyp + pB*pB + p0*p0)/(2.*pB*p0);
  }

  theta = std::acos(cost)*radian;

  G4LorentzVector p(0, 0, 0, TMath::Sqrt(p0*p0 + m0*m0));
  p.setRThetaPhi(p0, theta, phi);
  beam.VO(zK18Target);
  beam.pos.setZ(0);
  G4LorentzVector v(target_pos + beam.pos, 0);
#if 0
  beam.Print();
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


void // 7501 E75 phase-1 7XiH+6XiH spectrum kinematics 
S2SPrimaryGeneratorAction::GenerateKH7XiHSpectrum(G4Event* anEvent)
{
  static const G4int n_particle = 1;
  m_particleGun = new G4ParticleGun(n_particle);
  static const G4String name = "kaon+";
  static const auto particle = particleTable->FindParticle(name);
  static const auto pdg = particle->GetPDGEncoding();
  static const auto& target_pos = geomMan.GetGlobalPosition("Target")*mm;
  static const auto& target_size = sizeMan.GetSize("Target")*mm/2;
  static const G4double m0 = particle->GetPDGMass();
  static const G4int experiment = confMan.Get<G4int>("Experiment");
  G4double theta;
  //  std::acos(G4RandFlat::shoot(std::cos(0*deg), std::cos(20*deg)))*radian;
  //
  G4double pB, x0, y0, z0;
  G4double u0, v0;
  G4int count =0;
  while(count<1){
    G4int i = int(G4UniformRand()*(k18track->GetEntries()));
    k18track->GetEntry(i);
    if(trigflag[21]>0 && ntK18>0 && abs(CBtof0[0])<0.2 && chisqrK18[0]<10.){
    //z0 = G4RandFlat::shoot(-target_size.z()/2, target_size.z()/2);
    z0 = G4RandFlat::shoot(-187.*mm/2, 187*mm/2);
    x0 = xtgt[0] + utgt[0]*z0; // Horizontal direction
    y0 = ytgt[0] + vtgt[0]*z0; // vertical direction
    pB = p_3rd[0]*GeV;
    count++;
   }
  }
  G4double phi = G4RandFlat::shoot(0., 360.)*deg;
  G4double p0=0.;
  G4double BE=0.;

  G4double r = Sum*G4UniformRand(); 
  G4double yint =0;
  for(Int_t i=0;i<gr->GetN();++i){
    Double_t x, y;
    gr->GetPoint(i,x,y);
    yint += y;
    if(yint>r){
      BE=x;
      break;
    }
  }

  G4double mass_7Li = 6.53468*GeV; // 7Li
  G4double mass_7XiH = (5.61153+1.32171+BE*1e-3)*GeV; // 6He + Xi
  G4double mT = mass_7Li; // mass of target nucleus (7Li)
  G4double mY = mass_7XiH; // mass of target nucleus (7Li)

  G4double cost = 0.92;
  G4int count_pt=0;

  while(cost>1. || cost<0.93){
    p0 = G4RandFlat::shoot(1.27, 1.47)*GeV; // scatter momentum
    G4double EB = TMath::Sqrt(pB*pB + m0*m0);
    G4double ES = TMath::Sqrt(p0*p0 + m0*m0);
    G4double EY = EB + mT - ES;
    cost = (mY*mY - EY*EY + pB*pB + p0*p0)/(2.*pB*p0);
    count_pt ++;
    if(count_pt >10){
      cost =0.5;
      p0 = 0.0*GeV;
      break;
    }
  }
  std::cout <<"BE: "<< BE <<", pB: "<< pB<<", pS: "<< p0 <<", Theta: "<<cost<< std::endl;  
  theta = std::acos(cost)*radian;

  G4LorentzVector p(0, 0, 0, TMath::Sqrt(p0*p0 + m0*m0));
  p.setRThetaPhi(p0, theta, phi);
  beam.VO(zK18Target);
  
  beam.pos.setX(x0);
  beam.pos.setY(y0);
  beam.pos.setZ(z0);
  G4LorentzVector v(target_pos + beam.pos, 0);
  u0 = TMath::Tan(theta)*TMath::Cos(phi);
  v0 = TMath::Tan(theta)*TMath::Sin(phi);
#if 0
  beam.Print();
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
  anaMan.SetPrimaryData(x0,y0,z0,u0,v0,phi,theta,p0,pB,9999);
}
