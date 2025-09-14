// -*- C++ -*-

#include "S2SPrimaryGeneratorAction.hh"

#include <G4Event.hh>
#include <G4ParticleGun.hh>
#include <G4ParticleTable.hh>
#include <G4ParticleDefinition.hh>
#include <G4LorentzVector.hh>
#include <G4ThreeVector.hh>
#include <Randomize.hh>
#include <G4RandomDirection.hh> // E63

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
#include "FdComplex.hh"
#include "FermiMotion.hh"
#include "CMSMomentum.hh"


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
  case 5: GenerateBeamGausProfile(anEvent); break;
  case 6: GenerateBeamFixSeed(anEvent); break;
  case 7: GenerateScatParticles(anEvent); break;
  case 7001: Generate12XiBeryllium(anEvent); break;
  case 7002: GenerateElementaryXiMinus(anEvent); break;
  case 7003: GenerateElementarySigmaMinus(anEvent); break;
  case 7004: GenerateElementarySigmaPlus(anEvent); break;
    //case 7005: Generate12XiBePeakStructure(anEvent); break;
  case 7501: GenerateKH7XiHSpectrum(anEvent); break;
  case 6301: GenerateE63_7LambdaLi(anEvent, 7); break;
  case 6302: GenerateE63_7LambdaLi(anEvent, 10); break;
  case 6303: GenerateE63_7LambdaLi(anEvent, 12); break;
  case 9001: GenerateSigmaNCusp(anEvent); break;
  case 9002: GenerateQFLambda(anEvent); break;
  case 9003: GenerateQFSigmaZ(anEvent); break;
  case 9004: GenerateQFSigmaP(anEvent); break;
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
  // p0 = (experiment == 90) ?
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
  if(experiment==90)
    p0 = G4RandFlat::shoot(0.9, 1.5)*GeV;

  G4double theta =
    std::acos(G4RandFlat::shoot(std::cos(0*deg), std::cos(20*deg)))*radian;
  G4double phi = G4RandFlat::shoot(0., 360.)*deg;
  G4LorentzVector p(0, 0, 0, TMath::Sqrt(p0*p0 + m0*m0));
  p.setRThetaPhi(p0, theta, phi);
  G4double u0, v0;
  u0 = TMath::Tan(theta)*TMath::Cos(phi);
  v0 = TMath::Tan(theta)*TMath::Sin(phi);
  beam.VO(zK18Target);

  if(experiment!=90){
    beam.pos.setX(target_pos.x());
    beam.pos.setY(target_pos.y());
    beam.pos.setZ(target_pos.z());
  }
  if(experiment==90){
    double beam_x = G4RandGauss::shoot(target_pos.x(),23.);
    while(1){
      if(fabs(beam_x)<(54./2.))
	break;
      else
	beam_x = G4RandGauss::shoot(target_pos.x(),23.);
    }
    beam.pos.setX(beam_x);
    beam.pos.setY(target_pos.y());
    beam.pos.setZ(target_pos.z());
  }

  G4LorentzVector v(beam.pos, 0);
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
  anaMan.SetPrimaryData(beam.pos.x(),beam.pos.y(),beam.pos.z(),u0,v0,phi,theta,p0,1.8*GeV,9999);
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
    if(p0 < 0) name = "pi-";
    if(p0 > 0) name = "pi+";
  }else{
    if(p0 < 0) name = "kaon-";
    if(p0 > 0) name = "kaon+";
  }
  static const auto particle = particleTable->FindParticle(name);
  static const auto pdg = particle->GetPDGEncoding();
  const G4double m0 = particle->GetPDGMass();
  const auto& target_pos = geomMan.GetGlobalPosition("Target");
  beam.VO(zK18Target);
  beam.pos.setX(target_pos.x());
  beam.pos.setY(target_pos.y());
  beam.pos.setZ(target_pos.z());
  beam.mom.setMag(p0);
  G4LorentzVector p(beam.mom, TMath::Sqrt(p0*p0 + m0*m0));
  G4LorentzVector v(beam.pos, 0);
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
  beam.pos.setX(x0+target_pos.x());
  beam.pos.setY(y0+target_pos.y());
  beam.pos.setZ(z0+target_pos.z());
  G4LorentzVector p(0, 0, 0, TMath::Sqrt(p0*p0 + m0*m0));
  p.setX(p0*u0/TMath::Sqrt(1+u0*u0+v0*v0));
  p.setY(p0*v0/TMath::Sqrt(1+u0*u0+v0*v0));
  p.setZ(p0/TMath::Sqrt(1+u0*u0+v0*v0));
  G4LorentzVector v(beam.pos, 0);
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
  anaMan.SetPrimaryData(beam.pos.x(),beam.pos.y(),beam.pos.z(),u0,v0,0.,0.,p0,p0,9999);
}

void // 5
S2SPrimaryGeneratorAction::GenerateBeamGausProfile(G4Event* anEvent)
{
  static const G4int n_particle = 1;
  m_particleGun = new G4ParticleGun(n_particle);
  static const G4String name = "kaon+";
  static const auto particle = particleTable->FindParticle(name);
  static const auto pdg = particle->GetPDGEncoding();
  static const auto& target_pos = geomMan.GetGlobalPosition("Target")*mm;
  static const G4double m0 = particle->GetPDGMass();
  static const G4int experiment = confMan.Get<G4int>("Experiment");

  G4double BeamSizeX = 10.*mm;
  G4double BeamSizeY = 5.*mm;
  G4double BeamSizeU = 5.e-3;
  G4double BeamSizeV = 3.e-3;
  G4double MomCenter = 1.4*GeV ;
  G4double MomSize = 0.1*GeV ;

  G4double x0 = 0.0;
  G4double y0 = 0.0;
  G4double z0 = 0.0;
  G4double u0 = 0.0;
  G4double v0 = 0.0;
  G4double p0 = 0.0;

  x0 = G4RandGauss::shoot(target_pos.x(),BeamSizeX);
  y0 = G4RandGauss::shoot(target_pos.y(),BeamSizeY);
  z0 = target_pos.z();

  u0 = G4RandGauss::shoot(0.0,BeamSizeU);
  v0 = G4RandGauss::shoot(0.0,BeamSizeV);

  p0 = G4RandFlat::shoot(MomCenter-MomSize, MomCenter+MomSize);

  beam.VO(zK18Target);
  beam.pos.setX(x0);
  beam.pos.setY(y0);
  beam.pos.setZ(z0);
  G4LorentzVector p(0, 0, 0, TMath::Sqrt(p0*p0 + m0*m0));
  p.setX(p0*u0/TMath::Sqrt(1+u0*u0+v0*v0));
  p.setY(p0*v0/TMath::Sqrt(1+u0*u0+v0*v0));
  p.setZ(p0/TMath::Sqrt(1+u0*u0+v0*v0));
  G4LorentzVector v(beam.pos, 0);
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
  anaMan.SetPrimaryData(beam.pos.x(),beam.pos.y(),beam.pos.z(),u0,v0,0.,0.,p0,p0,9999);
}

void // 6
S2SPrimaryGeneratorAction::GenerateBeamFixSeed(G4Event* anEvent)
{
  static const G4int n_particle = 1;
  m_particleGun = new G4ParticleGun(n_particle);
  static const G4String name = "kaon+";
  static const auto particle = particleTable->FindParticle(name);
  static const auto pdg = particle->GetPDGEncoding();
  static const auto& target_pos = geomMan.GetGlobalPosition("Target")*mm;
  static const G4double m0 = particle->GetPDGMass();
  static const G4int experiment = confMan.Get<G4int>("Experiment");

  // G4double BeamSizeX = 10.*mm;
  // G4double BeamSizeY = 5.*mm;
  G4double BeamSizeX = 30.*mm;
  G4double BeamSizeY = 30.*mm;

  // G4double BeamSizeU = 5.e-3;
  // G4double BeamSizeV = 3.e-3;
  G4double BeamSizeU = 0.1;
  G4double BeamSizeV = 0.1;
  G4double ThetaSize =3.;
  G4double MomCenter = 1.35*GeV ;
  G4double MomSize = 0.*GeV ;

  G4double x0 = 0.0;
  G4double y0 = 0.0;
  G4double z0 = 0.0;
  G4double u0 = 0.0;
  G4double v0 = 0.0;
  G4double p0 = 0.0;

  ++inum;
  G4Random::setTheSeed(inum);
  //  G4cout << "   Seed = " << G4Random::getTheSeed() << G4endl;
  x0 = G4RandGauss::shoot(target_pos.x(),BeamSizeX);
  y0 = G4RandGauss::shoot(target_pos.y(),BeamSizeY);
  z0 = target_pos.z();

  // G4double theta =
  //   std::acos(G4RandFlat::shoot(std::cos(0*deg), std::cos(ThetaSize*deg)))*radian;
  // G4double phi = G4RandFlat::shoot(0., 360.)*deg;

  G4double theta =G4RandFlat::shoot(0., ThetaSize)*deg;
  G4double phi = G4RandFlat::shoot(0., 360.)*deg;

  //  std::cout<<"x:"<<x0<<", y:"<<y0<<", theta:"<<theta/deg<<", phi:"<<phi/deg<<std::endl;
  // u0 = G4RandGauss::shoot(0.0,BeamSizeU);
  // v0 = G4RandGauss::shoot(0.0,BeamSizeV);

  u0 = TMath::Tan(theta)*TMath::Cos(phi);
  v0 = TMath::Tan(theta)*TMath::Sin(phi);

  p0 = G4RandFlat::shoot(MomCenter-MomSize, MomCenter+MomSize);


  beam.VO(zK18Target);
  beam.pos.setX(x0);
  beam.pos.setY(y0);
  beam.pos.setZ(z0);
  G4LorentzVector p(0, 0, 0, TMath::Sqrt(p0*p0 + m0*m0));
  p.setX(p0*u0/TMath::Sqrt(1+u0*u0+v0*v0));
  p.setY(p0*v0/TMath::Sqrt(1+u0*u0+v0*v0));
  p.setZ(p0/TMath::Sqrt(1+u0*u0+v0*v0));
  //  G4LorentzVector v(beam.pos+target_pos, 0);
  G4LorentzVector v(beam.pos, 0);
  //  std::cout<<"target:"<<target_pos<<", beam:"<<beam.pos<<", v:"<<v.v()<<std::endl;
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
  //anaMan.SetPrimaryData(x0,y0,z0,u0,v0,0.,0.,p0,p0,9999);
  anaMan.SetPrimaryData(x0,y0,z0,u0,v0,phi,theta,p0,p0,9999);
}

void // 7
S2SPrimaryGeneratorAction::GenerateScatParticles(G4Event* anEvent)
{
  static const G4int n_particle = 1;
  m_particleGun = new G4ParticleGun(n_particle);
  static const G4String name = "kaon+";
  static const auto particle = particleTable->FindParticle(name);
  static const auto pdg = particle->GetPDGEncoding();
  static const auto& target_pos = geomMan.GetGlobalPosition("Target")*mm;
  static const G4double m0 = particle->GetPDGMass();
  static const G4int experiment = confMan.Get<G4int>("Experiment");

  G4double BeamSizeX = 10.*mm;
  G4double BeamSizeY = 5.*mm;
  // G4double BeamSizeU = 5.e-3;
  // G4double BeamSizeV = 3.e-3;
  G4double BeamSizeU = 0.1;
  G4double BeamSizeV = 0.1;
  G4double ThetaSize =15.;
  G4double MomCenter = 1.2*GeV ;
  G4double MomSize = 0.3*GeV ;

  G4double x0 = 0.0;
  G4double y0 = 0.0;
  G4double z0 = 0.0;
  G4double u0 = 0.0;
  G4double v0 = 0.0;
  G4double p0 = 0.0;

  x0 = G4RandGauss::shoot(target_pos.x(),BeamSizeX);
  y0 = G4RandGauss::shoot(target_pos.y(),BeamSizeY);
  z0 = target_pos.z();

  // G4double theta =
  //   std::acos(G4RandFlat::shoot(std::cos(0*deg), std::cos(ThetaSize*deg)))*radian;
  // G4double phi = G4RandFlat::shoot(0., 360.)*deg;

  G4double theta =G4RandFlat::shoot(0., ThetaSize)*deg;
  G4double phi = G4RandFlat::shoot(0., 360.)*deg;

  // u0 = G4RandGauss::shoot(0.0,BeamSizeU);
  // v0 = G4RandGauss::shoot(0.0,BeamSizeV);

  u0 = TMath::Tan(theta)*TMath::Cos(phi);
  v0 = TMath::Tan(theta)*TMath::Sin(phi);

  p0 = G4RandFlat::shoot(MomCenter-MomSize, MomCenter+MomSize);

  beam.VO(zK18Target);
  beam.pos.setX(x0);
  beam.pos.setY(y0);
  beam.pos.setZ(z0);
  G4LorentzVector p(0, 0, 0, TMath::Sqrt(p0*p0 + m0*m0));
  p.setX(p0*u0/TMath::Sqrt(1+u0*u0+v0*v0));
  p.setY(p0*v0/TMath::Sqrt(1+u0*u0+v0*v0));
  p.setZ(p0/TMath::Sqrt(1+u0*u0+v0*v0));
  //  G4LorentzVector v(beam.pos+target_pos, 0);
  G4LorentzVector v(beam.pos, 0);
  //  std::cout<<"target:"<<target_pos<<", beam:"<<beam.pos<<", v:"<<v.v()<<std::endl;
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
  //anaMan.SetPrimaryData(x0,y0,z0,u0,v0,0.,0.,p0,p0,9999);
  anaMan.SetPrimaryData(x0,y0,z0,u0,v0,phi,theta,p0,p0,9999);
}

void // 7001 E70 12C(KK)12XiBe kinematics
S2SPrimaryGeneratorAction::Generate12XiBeryllium(G4Event* anEvent)
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
  beam.pos.setX(target_pos.x());
  beam.pos.setY(target_pos.y());
  beam.pos.setZ(target_pos.z());
  G4LorentzVector v(beam.pos, 0);
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
  beam.pos.setX(target_pos.x());
  beam.pos.setY(target_pos.y());
  beam.pos.setZ(target_pos.z());
  G4LorentzVector v(beam.pos, 0);
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
  beam.pos.setX(target_pos.x());
  beam.pos.setY(target_pos.y());
  beam.pos.setZ(target_pos.z());
  G4LorentzVector v(beam.pos, 0);
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
  beam.pos.setX(target_pos.x());
  beam.pos.setY(target_pos.y());
  beam.pos.setZ(target_pos.z());
  G4LorentzVector v(beam.pos, 0);
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

  beam.pos.setX(x0+target_pos.x());
  beam.pos.setY(y0+target_pos.y());
  beam.pos.setZ(z0+target_pos.z());
  G4LorentzVector v(beam.pos, 0);
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
  anaMan.SetPrimaryData(beam.pos.x(),beam.pos.y(),beam.pos.z(),u0,v0,phi,theta,p0,pB,9999);
}

//_____________________________________________________________________________
// E63
void // 6301~  [ E63 A(K-,pi-)lambda_hyper kinematics ] 
S2SPrimaryGeneratorAction::GenerateE63_7LambdaLi(G4Event* anEvent, G4int MassNum)
{
  static const G4int n_particle = 1;  // should 1 even if you generate weak pion
  m_particleGun = new G4ParticleGun(n_particle);
  static const auto& target_pos = geomMan.GetGlobalPosition("Target")*mm;
  static const auto& target_size = sizeMan.GetSize("Target")*mm/2;
  static const G4int experiment = confMan.Get<G4int>("Experiment");
  auto WeakParticle = confMan.Get<G4String>("WeakDecayParticle");
  int dummy_event_flag = 0;
  int dummy_event_flag_out_of_tgt = 0; // out of target
  int dummy_event_flag_out_of_cs = 0; // out of cross section
  
  // ***************
  // *** beam K- ***
  static const auto beam_particle = particleTable->FindParticle("kaon-");
  static const auto beam_pdg =      beam_particle->GetPDGEncoding();
  static const G4double m_beam =    beam_particle->GetPDGMass();
  // momentum
  G4double p_beam = confMan.Get<G4double>("PK18")*CLHEP::GeV;
  if(1){ // include BeamMomentuBite
    p_beam += G4RandGauss::shoot( 0.0,  p_beam*0.0134 );  // momentum bite [202501 data]
  }
  // momentum direction
  G4double beam_u_rms = 0.0164; // [202501 data]
  G4double beam_v_rms = 0.0045;
  G4double beam_u = G4RandGauss::shoot( 0.0, beam_u_rms );
  G4double beam_v = G4RandGauss::shoot( 0.0, beam_v_rms );
  G4double beam_z = 1./sqrt(1. + beam_u*beam_u + beam_v*beam_v);
  G4ThreeVector BeamMomDir = G4ThreeVector( beam_u*beam_z, beam_v*beam_z, beam_z );
  // profile at target center -> vertex at terget coodinate
  G4double beam_x_rms = 16.94; // [202501 data]
  G4double beam_y_rms =  7.01;
  G4double vertex_x = G4RandGauss::shoot( 0.0, beam_x_rms );
  G4double vertex_y = G4RandGauss::shoot( 0.0, beam_y_rms );
  G4double vertex_z = G4RandFlat::shoot( -target_size.z(), target_size.z() );
  G4ThreeVector VertexPos = G4ThreeVector( vertex_x, vertex_y, vertex_z );
  G4LorentzVector VertexLv(VertexPos, 0);

  if(1){  // InsideTarget cut
    if( fabs(vertex_x)>fabs(target_size.x()) || fabs(vertex_y)>fabs(target_size.y()) ){
      //G4cerr << "out of target" << G4endl;
      dummy_event_flag_out_of_tgt = 1;  // set dummy vertex if vertex is outof target_size
    }
  }

  // ****************
  // *** scat pi- ***
  static const auto scat_particle = particleTable->FindParticle("pi-");
  static const auto scat_pdg =      scat_particle->GetPDGEncoding();
  static const G4double m_scat =    scat_particle->GetPDGMass();
  G4double costLab;
  // limit of scat angle
  {
    //costLab = 1. - G4RandFlat::shoot(0., 0.003805); // (0- 5 deg.)
    costLab = 1. - G4RandFlat::shoot(0., 0.015192); // (0-10 deg.)
    //costLab = 1. - G4RandFlat::shoot(0., 0.034074); // (0-15 deg.)
    //costLab = 1. - G4RandFlat::shoot(0., 0.060307); // (0-20 deg.)
    //costLab = 1. - G4RandFlat::shoot(0., 0.093692); // (0-25 deg.)
    //costLab = 1. - G4RandFlat::shoot(0., 0.133974); // (0-30 deg.)
  }
  
  if(1){ // cross section shape cut
    G4int DeltaL = 0;  // need 0 or 1 or 2
    G4double p[6]; // f_cross = pol(6)
    G4double MaxCrossSection = 1100; // [a.u.]
    if(DeltaL==0){  // table for 0.9 GeV/c Li (k.pi)
      p[0]=1033; p[1]=70.55; p[2]=-59.28;
      p[3]=7.2325; p[4]=-0.37978; p[5]=0.00927192;
      p[6]=-8.56922e-5;
    }
    else if(DeltaL==1){
      p[0]=259.9; p[1]=-38.1581; p[2]=27.6301;
      p[3]=-3.38422; p[4]=0.149127; p[5]=-0.00233459;
      p[6]=3.37185e-6;
    }
    else if(DeltaL==2){
      p[0]=140.35; p[1]=31.2678; p[2]=-17.7024;
      p[3]=4.27016; p[4]=-0.3732; p[5]=0.0135775;
      p[6]=-0.000177677;
    }
    else{
      G4cerr << "DeltaL setting is wrong" << G4endl;
      exit(-1);
    }

    G4double thetaLab = acos(costLab)*(180./3.141592); // degree
    G4double RandValue = G4RandFlat::shoot(0., MaxCrossSection);
    G4double cross_section = p[0];
    for(int n=1; n<7; n++){
      cross_section += p[n]*pow(thetaLab,n);
    }
    if( RandValue>cross_section ){
      dummy_event_flag_out_of_cs = 1; // set dummy vertex if RandValue > cross section table
    }
  } // cross section cut
  
  // ****************
  // *** weak pi- ***
  static const auto weak_particle = particleTable->FindParticle("pi-");
  static const auto weak_pdg =      weak_particle->GetPDGEncoding();
  static const G4double m_weak =    weak_particle->GetPDGMass();
  G4double WeakT;
  {
    if(WeakParticle == "3LH") WeakT = 40.88 *CLHEP::MeV ; // for 3LH
    if(WeakParticle == "4LH") WeakT = 53.25 *CLHEP::MeV ; // for 4LH
    if(WeakParticle == "6LH") WeakT = 37.20 *CLHEP::MeV ; // for 6LHe
    //else WeakT = 40. *CLHEP::MeV ; //
  }
  G4double WeakMom = sqrt( pow(WeakT+m_weak,2) -m_scat*m_scat );
  
  
  // *************************************
  // *** target nuclei and hypernuclei ***
  G4double AtomicMassUnit = 0.93149432;
  G4double LambdaMass = particleTable->FindParticle("lambda")->GetPDGMass();
  G4double m_tgt, m_hyp;
  {
    // 7Li -----------
    G4double mass_7Li = (7.0*AtomicMassUnit+0.014908)*GeV;
    //G4double mass_7LambdaLi = (6.0*AtomicMassUnit+0.014086-0.00522+0.000)*GeV + LambdaMass; // Ex=0 MeV
    G4double mass_7LambdaLi = (6.0*AtomicMassUnit+0.014086-0.00522+0.020)*GeV + LambdaMass; // Ex=20 MeV
    
    // 10B -----------
    G4double mass_10B = (10.0*AtomicMassUnit+0.0120508)*GeV;
    G4double mass_10LambdaB = (9.0*AtomicMassUnit+0.0113477-0.0081+0.000)*GeV + LambdaMass; // Ex=0 MeV

    // 12C -----------
    G4double mass_12C = (12.0*AtomicMassUnit+0.0)*GeV;
    G4double mass_12LambdaC = (11.0*AtomicMassUnit+0.010650-0.0108+0.000)*GeV + LambdaMass; // Ex=0 MeV

    if(MassNum==7){
      m_tgt = mass_7Li; 
      m_hyp = mass_7LambdaLi;
    }
    else if(MassNum==10){
      m_tgt = mass_10B; 
      m_hyp = mass_10LambdaB;
    }
    else if(MassNum==12){
      m_tgt = mass_12C; 
      m_hyp = mass_12LambdaC;
    }
    else{
      G4cerr << "not regstered MassNumber  MassNum = " << MassNum << G4endl;
      exit(-1);
    }
  }
  
  
  // ************************
  // **** calculate scat ****
  // ************************
  
  //Kaon 1.5GeV/c
  G4LorentzVector BeamLv( p_beam*BeamMomDir, 
			  sqrt( m_beam*m_beam+p_beam*p_beam ) );
  //Neutron 0.0GeV/c
  G4double NuclMom = 0.0;
  G4ThreeVector NuclMomDir( 0., 0., NuclMom );
  G4LorentzVector NuclLv( NuclMom*NuclMomDir,
			  sqrt( m_tgt*m_tgt+NuclMom*NuclMom )  );

  //Primary frame
  G4LorentzVector PrimaryLv =  BeamLv+NuclLv;
  G4double TotalEnergyCM = PrimaryLv.mag();
  G4ThreeVector beta( PrimaryLv.vect()/PrimaryLv.e() ); 

  //scat CM
  G4double ScatMomCM 
    = 0.5*sqrt(( TotalEnergyCM*TotalEnergyCM
  		 -( m_scat+m_hyp )*( m_scat+m_hyp ))
  	       *( TotalEnergyCM*TotalEnergyCM
  		  -( m_scat-m_hyp )*( m_scat-m_hyp )))/TotalEnergyCM;

  G4double cottLab=costLab/sqrt(1.-costLab*costLab);
  G4double bt=beta.mag(), gamma=1./sqrt(1.-bt*bt);
  G4double gbep=gamma*bt*sqrt(ScatMomCM*ScatMomCM+m_scat*m_scat)/ScatMomCM;
  G4double a  = gamma*gamma+cottLab*cottLab;
  G4double bp = gamma*gbep;
  G4double c  = gbep*gbep-cottLab*cottLab;

  G4double dd=bp*bp-a*c;
  if( dd<0. ){
    G4cerr << "dd<0." << G4endl;
    exit(-1);
  }

  G4double costCM=(sqrt(dd)-bp)/a;
  if( costCM>1. || costCM<-1. ){
    G4cerr << "costCM>1. || costCM<-1." << G4endl;
    exit(-1);
  }

  G4double sintCM=sqrt(1.-costCM*costCM);
  G4double phiCM=G4RandFlat::shoot(0., 360.)*deg;
  G4ThreeVector ScatMomCM_vector( ScatMomCM*sintCM*cos(phiCM),
				  ScatMomCM*sintCM*sin(phiCM), 
				  ScatMomCM*costCM );
  //ScatMomCM_vector.rotateY(KaonMomDir.theta()); 
  //ScatMomCM_vector.rotateZ(KaonMomDir.phi());
  ScatMomCM_vector.rotateUz(BeamMomDir);
  
  G4LorentzVector ScatLv( ScatMomCM_vector, 
			  sqrt( ScatMomCM*ScatMomCM + m_scat*m_scat ));
  ScatLv.boost(beta);
  
  G4ThreeVector ScatMom_vector = ScatLv.vect();
  G4double ScatMom = ScatMom_vector.mag();
  G4ThreeVector ScatMomDir = ScatMom_vector/ScatMom;
  G4double ScatT = sqrt( ScatMom*ScatMom + m_scat*m_scat ) - m_scat;

  
  // ************************
  // **** calculate weak ****
  // ************************
  G4ThreeVector WeakMomDir = G4RandomDirection();
  G4ThreeVector WeakMom_vector = WeakMom * WeakMomDir;
  G4LorentzVector WeakLv( WeakMom_vector, WeakT);

  
  // ***************************
  // **** particle generate ****
  // ***************************

  // ****************************
  // ** scat particle generate **
  if(dummy_event_flag_out_of_tgt || dummy_event_flag_out_of_cs) dummy_event_flag = 1;
  if(dummy_event_flag){ // set primary vertex at out of world
    ScatMomDir = G4ThreeVector(0, 0, -1);
    ScatMom = 0; ScatT = 0;
    VertexPos = G4ThreeVector(-200*CLHEP::m, -200*CLHEP::m, -200*CLHEP::m);
  }
  m_particleGun->SetParticleDefinition(scat_particle);
  m_particleGun->SetParticleMomentumDirection(ScatMomDir);
  m_particleGun->SetParticleEnergy(ScatT);
  m_particleGun->SetParticlePosition(VertexPos + target_pos);
  m_particleGun->GeneratePrimaryVertex(anEvent);
  // analyzer fill : scat particle
  {
    G4double x0 = VertexPos.x(); G4double y0 = VertexPos.y(); G4double z0 = VertexPos.z();
    G4double u0 = ScatMomDir.x()/ScatMomDir.z(); G4double v0 = ScatMomDir.y()/ScatMomDir.z();
    G4double phi = ScatMomDir.phi(); G4double theta = ScatMomDir.theta();
    G4double p0 = ScatMom; G4double pB = ScatT;
    anaMan.SetPrimaryData(x0,y0,z0,u0,v0,phi,theta,p0,pB,9999);
    // unit: x0,y0,z0=mm, phi,theta=rad, p0[momentum]=MeV/c, pB[KineticEnergy]=MeV
    anaMan.SetPrimaryParticle(0, scat_pdg, ScatLv, VertexLv);
  }

  // **********************************
  // ** weak decay particle generate **
  if(1){
    if(dummy_event_flag){ // set primary vertex at out of world
      WeakMomDir = G4ThreeVector(0, 0, -1);
      WeakMom = 0; WeakT = 0;
      if(dummy_event_flag_out_of_tgt){
	VertexPos = G4ThreeVector(-200*CLHEP::m, -200*CLHEP::m, -200*CLHEP::m);
	if(dummy_event_flag_out_of_cs)
	VertexPos = G4ThreeVector(-190*CLHEP::m, -190*CLHEP::m, -190*CLHEP::m);
      }
      if(dummy_event_flag_out_of_cs)
	VertexPos = G4ThreeVector(-180*CLHEP::m, -180*CLHEP::m, -180*CLHEP::m);
    }

    m_particleGun->SetParticleDefinition(weak_particle);
    m_particleGun->SetParticleMomentumDirection(WeakMomDir);
    m_particleGun->SetParticleEnergy(WeakT);
    m_particleGun->SetParticlePosition(VertexPos + target_pos);
    m_particleGun->GeneratePrimaryVertex(anEvent);
    // analyzer fill : weak decay particle
    {
      G4double x1 = VertexPos.x(); G4double y1 = VertexPos.y(); G4double z1 = VertexPos.z();
      G4double px1 = WeakMom_vector.x(); G4double py1 = WeakMom_vector.y(); G4double pz1 = WeakMom_vector.z();
      G4double p1 = WeakMom; G4double t1 = WeakT;
      anaMan.SetSecondaryData(x1,y1,z1,px1,py1,pz1,p1,t1);
      // unit: x1,y1,z1=mm, px1,py1,pz1,p1[momentum]=MeV/c, t1[KineticEnergy]=MeV
    }
  } // if weak decay particle
  
}

//_____________________________________________________________________________
void // 9001 E90 SigmaN Cusp Scattering
S2SPrimaryGeneratorAction::GenerateSigmaNCusp(G4Event* anEvent)
{
  static const G4int n_particle = 1;
  m_particleGun = new G4ParticleGun(n_particle);

  //========== Define particles ==========
  const auto pTable = G4ParticleTable::GetParticleTable();
  const auto kaon_minus = pTable->FindParticle("kaon-");
  const auto pi_minus   = pTable->FindParticle("pi-");
  const auto proton     = pTable->FindParticle("proton");
  const auto neutron    = pTable->FindParticle("neutron");
  const auto deuteron   = pTable->FindParticle("deuteron");
  const auto lambda     = pTable->FindParticle("lambda");
  const auto sigma_plus = pTable->FindParticle("sigma+");
  
  const G4double M_Kaon     = kaon_minus->GetPDGMass();
  const G4double M_PiM       = pi_minus->GetPDGMass();
  const G4double M_Proton   = proton->GetPDGMass();
  const G4double M_Deuteron = deuteron->GetPDGMass();
  const G4double M_Lambda   = lambda->GetPDGMass();
  const G4double M_SigmaP   = sigma_plus->GetPDGMass();
  const G4double M_Neutron  = neutron->GetPDGMass();

  //========== Calculate maximum ΣN cusp mass ==========
  const G4double scat_a = 2.06;
  const G4double scat_b = 4.64;
  const G4double threshold_mass = M_SigmaP + M_Neutron;
  const G4double mass_min = threshold_mass - 100. * CLHEP::MeV;
  const G4double mass_max = threshold_mass + 100. * CLHEP::MeV;

  G4double fcusp_max = 0.;
  for(int i=0; i<200; ++i){
    double mass = mass_min + (mass_max-mass_min)*i/199.;
    double val = FdComplex::f_single(mass/CLHEP::MeV, scat_a, scat_b, 186.);
    if(val > fcusp_max) fcusp_max = val;
  }

  // loop until proper event is generated
  while(true){

    //--- Generate ΣN (X) mass ---
    G4double CuspM = 0.;
    while(true){
      G4double MM = G4RandFlat::shoot(mass_min, mass_max);
      G4double ds_MM = FdComplex::f_single(MM/CLHEP::MeV, scat_a, scat_b, 186.);
      G4double Rand = G4RandFlat::shoot(0., fcusp_max * 1.3);
      if(Rand <= ds_MM){
        CuspM = MM;
        break;
      }
    }

    //========== Kinematics ==========
    //--- 1. K- + d -> pi- + X ---
    G4double beam_mom_mean = 1.4 * CLHEP::GeV;
    G4double beam_mom_sigma = beam_mom_mean * 3.3e-4 / 2.355;
    G4double beam_mom = G4RandGauss::shoot(beam_mom_mean, beam_mom_sigma);
    G4LorentzVector beam_lv(0, 0, beam_mom, sqrt(beam_mom*beam_mom + M_Kaon*M_Kaon));
    TVector3 p_fermi_vec = FermiMotion::GetMomentum();
    G4ThreeVector p_deuteron_lab(p_fermi_vec.X(), p_fermi_vec.Y(), p_fermi_vec.Z());
    G4double E_deuteron = sqrt(M_Deuteron*M_Deuteron + p_deuteron_lab.mag2());
    G4LorentzVector target_lv(p_deuteron_lab, E_deuteron);

    G4LorentzVector W = beam_lv + target_lv;
    G4ThreeVector beta_cm = W.boostVector();

    G4double p_cm = CMSMomentum::GetCMSMomentum(W.m(), M_PiM, CuspM);
    if (p_cm == 0.0) continue;

    G4ThreeVector p_pi_cm_dir = G4RandomDirection();
    G4LorentzVector pi_lv_cm(p_pi_cm_dir * p_cm, sqrt(p_cm*p_cm + M_PiM*M_PiM));

    G4LorentzVector pi_lv_lab = pi_lv_cm;
    pi_lv_lab.boost(beta_cm);

    if (pi_lv_lab.theta() > 15. * CLHEP::deg) continue; // S2s acceptance cut
    
    G4LorentzVector X_lv_lab = W - pi_lv_lab;

    //--- 2. X -> Λ + p ---
    G4ThreeVector beta_X = X_lv_lab.boostVector();
    G4double p_X_rest = CMSMomentum::GetCMSMomentum(X_lv_lab.m(), M_Lambda, M_Proton);
    if (p_X_rest == 0.0) continue;

    G4ThreeVector p_lambda_rest_dir = G4RandomDirection();
    G4LorentzVector lambda_lv_rest(p_lambda_rest_dir * p_X_rest, sqrt(p_X_rest*p_X_rest + M_Lambda*M_Lambda));
    G4LorentzVector p_secondary_lv_rest(-p_lambda_rest_dir * p_X_rest, sqrt(p_X_rest*p_X_rest + M_Proton*M_Proton));

    G4LorentzVector lambda_lv_lab = lambda_lv_rest;
    lambda_lv_lab.boost(beta_X);
    G4LorentzVector p_secondary_lv_lab = p_secondary_lv_rest;
    p_secondary_lv_lab.boost(beta_X);

    //--- 3. Λ -> p + pi- ---
    G4ThreeVector beta_lambda = lambda_lv_lab.boostVector();
    G4double p_lambda_rest = CMSMomentum::GetCMSMomentum(M_Lambda, M_Proton, M_PiM);

    G4ThreeVector p_decay_rest_dir = G4RandomDirection();
    G4LorentzVector p_decay_lv_rest(p_decay_rest_dir * p_lambda_rest, sqrt(p_lambda_rest*p_lambda_rest + M_Proton*M_Proton));
    G4LorentzVector pi_decay_lv_rest(-p_decay_rest_dir * p_lambda_rest, sqrt(p_lambda_rest*p_lambda_rest + M_PiM*M_PiM));

    G4LorentzVector p_decay_lv_lab = p_decay_lv_rest;
    p_decay_lv_lab.boost(beta_lambda);
    G4LorentzVector pi_decay_lv_lab = pi_decay_lv_rest;
    pi_decay_lv_lab.boost(beta_lambda);

    //========== Calculate Vertex position ==========
    //--- 1. Primary reaction ---
    static const auto& target_pos = geomMan.GetGlobalPosition("Target")*mm;
    static const auto& target_size = sizeMan.GetSize("Target");

    G4double target_radius = target_size.y() / 2.0 * mm;
    G4double target_height_half = target_size.z() / 2.0 * mm;

    G4double r_vtx_x, r_vtx_y, r_vtx_z;
    G4double beam_x_sigma = 0.2 * CLHEP::mm;
    G4double beam_y_sigma = 0.2 * CLHEP::mm;
    while (true) {
        r_vtx_x = G4RandGauss::shoot(0, beam_x_sigma);
        r_vtx_y = G4RandGauss::shoot(0, beam_y_sigma);
        r_vtx_z = G4RandFlat::shoot(-target_radius, target_radius);
        if (r_vtx_x*r_vtx_x + r_vtx_z*r_vtx_z < target_radius*target_radius && r_vtx_y < target_height_half && r_vtx_y > -target_height_half) break;
    }

    G4ThreeVector primary_vertex_pos = target_pos + G4ThreeVector(r_vtx_x, r_vtx_y, r_vtx_z);

    //--- 2. Secondary reaction ---
    G4double lambda_lifetime = lambda->GetPDGLifeTime(); 
    G4double decay_time_rest = -lambda_lifetime * log(G4RandFlat::shoot());
    G4double decay_time_lab = lambda_lv_lab.gamma() * decay_time_rest;

    G4ThreeVector velocity_vec = lambda_lv_lab.vect() * (CLHEP::c_light * CLHEP::c_light / lambda_lv_lab.e());
    G4ThreeVector flight_distance = velocity_vec * decay_time_lab;

    G4ThreeVector secondary_vertex_pos = primary_vertex_pos + flight_distance;

    //========== Particle Generation ==========
    //--- scattering pion  ---
    m_particleGun->SetParticleDefinition(pi_minus);
    m_particleGun->SetParticleMomentumDirection(pi_lv_lab.vect().unit());
    m_particleGun->SetParticleEnergy(pi_lv_lab.e() - M_PiM);
    m_particleGun->SetParticlePosition(primary_vertex_pos);
    m_particleGun->GeneratePrimaryVertex(anEvent);
    anaMan.SetPrimaryParticle(0, pi_minus->GetPDGEncoding(), pi_lv_lab, G4LorentzVector(primary_vertex_pos, 0));

    //--- secondary proton ---
    m_particleGun->SetParticleDefinition(proton);
    m_particleGun->SetParticleMomentumDirection(p_secondary_lv_lab.vect().unit());
    m_particleGun->SetParticleEnergy(p_secondary_lv_lab.e() - M_Proton);
    m_particleGun->SetParticlePosition(primary_vertex_pos);
    m_particleGun->GeneratePrimaryVertex(anEvent);

    // //--- lambda decay proton ---
    // m_particleGun->SetParticleDefinition(proton);
    // m_particleGun->SetParticleMomentumDirection(p_decay_lv_lab.vect().unit());
    // m_particleGun->SetParticleEnergy(p_decay_lv_lab.e() - M_Proton);
    // m_particleGun->SetParticlePosition(secondary_vertex_pos);
    // m_particleGun->GeneratePrimaryVertex(anEvent);

    // //--- decay pion ---
    // m_particleGun->SetParticleDefinition(pi_minus);
    // m_particleGun->SetParticleMomentumDirection(pi_decay_lv_lab.vect().unit());
    // m_particleGun->SetParticleEnergy(pi_decay_lv_lab.e() - M_PiM);
    // m_particleGun->SetParticlePosition(secondary_vertex_pos);
    // m_particleGun->GeneratePrimaryVertex(anEvent);
    
    break; 
  }
}

//_____________________________________________________________________________
void // 9002 E90 QF Lambda reaction
S2SPrimaryGeneratorAction::GenerateQFLambda(G4Event* anEvent)
{
  static const G4int n_particle = 1;
  m_particleGun = new G4ParticleGun(n_particle);

  //========== Define particles ==========
  const auto pTable = G4ParticleTable::GetParticleTable();
  const auto kaon_minus = pTable->FindParticle("kaon-");
  const auto pi_minus   = pTable->FindParticle("pi-");
  const auto proton     = pTable->FindParticle("proton");
  const auto neutron    = pTable->FindParticle("neutron");
  const auto lambda     = pTable->FindParticle("lambda");

  const G4double M_Kaon     = kaon_minus->GetPDGMass();
  const G4double M_PiM      = pi_minus->GetPDGMass();
  const G4double M_Proton   = proton->GetPDGMass();
  const G4double M_Neutron  = neutron->GetPDGMass();
  const G4double M_Lambda   = lambda->GetPDGMass();

  // loop until proper event is generated
  while(true){

    //--- Generate Fermi momentum for the target "n" ---
    TVector3 p_fermi_vec = FermiMotion::GetMomentum();
    G4ThreeVector target_mom_lab(p_fermi_vec.X(), p_fermi_vec.Y(), p_fermi_vec.Z());
    G4double E_target = sqrt(M_Neutron*M_Neutron + target_mom_lab.mag2());
    
    //========== Kinematics ==========
    //--- 1. K- + "n" -> Λ + pi- ---
    G4double beam_mom_mean = 1.4 * CLHEP::GeV;
    G4double beam_mom_sigma = beam_mom_mean * 3.3e-4 / 2.355;
    G4double beam_mom = G4RandGauss::shoot(beam_mom_mean, beam_mom_sigma);
    G4LorentzVector beam_lv(0, 0, beam_mom, sqrt(beam_mom*beam_mom + M_Kaon*M_Kaon));
    
    G4LorentzVector target_lv(target_mom_lab, E_target);
    G4LorentzVector W = beam_lv + target_lv;
    G4ThreeVector beta_cm = W.boostVector();

    G4double p_cm = CMSMomentum::GetCMSMomentum(W.m(), M_Lambda, M_PiM);
    if (p_cm == 0.0) continue;

    G4ThreeVector p_pi_cm_dir = G4RandomDirection();
    G4LorentzVector pi_lv_cm(p_pi_cm_dir * p_cm, sqrt(p_cm*p_cm + M_PiM*M_PiM));
    G4LorentzVector lambda_lv_cm(-p_pi_cm_dir * p_cm, sqrt(p_cm*p_cm + M_Lambda*M_Lambda));

    G4LorentzVector pi_lv_lab = pi_lv_cm;
    pi_lv_lab.boost(beta_cm);
    G4LorentzVector lambda_lv_lab = lambda_lv_cm;
    lambda_lv_lab.boost(beta_cm);
    
    if (pi_lv_lab.theta() > 15. * CLHEP::deg) continue; // S2s acceptance cut
    
    //--- 2. Λ -> p + pi- ---
    G4ThreeVector beta_lambda = lambda_lv_lab.boostVector();
    G4double p_lambda_rest = CMSMomentum::GetCMSMomentum(M_Lambda, M_Proton, M_PiM);

    G4ThreeVector p_decay_rest_dir = G4RandomDirection();
    G4LorentzVector p_decay_lv_rest(p_decay_rest_dir * p_lambda_rest, sqrt(p_lambda_rest*p_lambda_rest + M_Proton*M_Proton));
    G4LorentzVector pi_decay_lv_rest(-p_decay_rest_dir * p_lambda_rest, sqrt(p_lambda_rest*p_lambda_rest + M_PiM*M_PiM));

    G4LorentzVector p_decay_lv_lab = p_decay_lv_rest;
    p_decay_lv_lab.boost(beta_lambda);
    G4LorentzVector pi_decay_lv_lab = pi_decay_lv_rest;
    pi_decay_lv_lab.boost(beta_lambda);

    //========== Calculate Vertex position ==========
    static const auto& target_pos = geomMan.GetGlobalPosition("Target")*mm;
    static const auto& target_size = sizeMan.GetSize("Target");
    G4double target_radius = target_size.y() / 2.0 * mm;
    G4double target_height_half = target_size.z() / 2.0 * mm;

    G4double r_vtx_x, r_vtx_y, r_vtx_z;
    G4double beam_x_sigma = 0.2 * CLHEP::mm;
    G4double beam_y_sigma = 0.2 * CLHEP::mm;
    while (true) {
        r_vtx_x = G4RandGauss::shoot(0, beam_x_sigma);
        r_vtx_y = G4RandGauss::shoot(0, beam_y_sigma);
        r_vtx_z = G4RandFlat::shoot(-target_radius, target_radius);
        if (r_vtx_x*r_vtx_x + r_vtx_z*r_vtx_z < target_radius*target_radius && r_vtx_y < target_height_half && r_vtx_y > -target_height_half) break;
    }
    G4ThreeVector primary_vertex_pos = target_pos + G4ThreeVector(r_vtx_x, r_vtx_y, r_vtx_z);

    //--- 2. Secondary reaction (Λ decay) ---
    G4double lambda_lifetime = lambda->GetPDGLifeTime(); 
    G4double decay_time_rest = -lambda_lifetime * log(G4RandFlat::shoot());
    G4double decay_time_lab = lambda_lv_lab.gamma() * decay_time_rest;

    G4ThreeVector velocity_vec = lambda_lv_lab.vect() * (CLHEP::c_light * CLHEP::c_light / lambda_lv_lab.e());
    G4ThreeVector flight_distance = velocity_vec * decay_time_lab;

    G4ThreeVector secondary_vertex_pos = primary_vertex_pos + flight_distance;

    //========== Particle Generation ==========
    //--- scattering pion  ---
    m_particleGun->SetParticleDefinition(pi_minus);
    m_particleGun->SetParticleMomentumDirection(pi_lv_lab.vect().unit());
    m_particleGun->SetParticleEnergy(pi_lv_lab.e() - M_PiM);
    m_particleGun->SetParticlePosition(primary_vertex_pos);
    m_particleGun->GeneratePrimaryVertex(anEvent);
    
    //--- spectator proton   ---
    G4ThreeVector total_mom_before = beam_lv.vect() + target_lv.vect();
    G4ThreeVector total_mom_after_reaction = pi_lv_lab.vect() + lambda_lv_lab.vect();
    G4ThreeVector spectator_mom = total_mom_before - total_mom_after_reaction;

    G4double E_spectator = sqrt(M_Proton*M_Proton + spectator_mom.mag2());
    m_particleGun->SetParticleDefinition(proton);
    m_particleGun->SetParticleMomentumDirection(spectator_mom.unit());
    m_particleGun->SetParticleEnergy(E_spectator - M_Proton);
    m_particleGun->SetParticlePosition(primary_vertex_pos);
    m_particleGun->GeneratePrimaryVertex(anEvent);
    
    //--- Λ decay proton ---
    m_particleGun->SetParticleDefinition(proton);
    m_particleGun->SetParticleMomentumDirection(p_decay_lv_lab.vect().unit());
    m_particleGun->SetParticleEnergy(p_decay_lv_lab.e() - M_Proton);
    m_particleGun->SetParticlePosition(secondary_vertex_pos);
    m_particleGun->GeneratePrimaryVertex(anEvent);

    //--- Λ decay pion ---
    m_particleGun->SetParticleDefinition(pi_minus);
    m_particleGun->SetParticleMomentumDirection(pi_decay_lv_lab.vect().unit());
    m_particleGun->SetParticleEnergy(pi_decay_lv_lab.e() - M_PiM);
    m_particleGun->SetParticlePosition(secondary_vertex_pos);
    m_particleGun->GeneratePrimaryVertex(anEvent);
    
    break; 
  }
}

//_____________________________________________________________________________
void // 9003 E90 QF Sigma0 reaction
S2SPrimaryGeneratorAction::GenerateQFSigmaZ(G4Event* anEvent)
{
  static const G4int n_particle = 1;
  m_particleGun = new G4ParticleGun(n_particle);

  //========== Define particles ==========
  const auto pTable = G4ParticleTable::GetParticleTable();
  const auto kaon_minus = pTable->FindParticle("kaon-");
  const auto pi_minus   = pTable->FindParticle("pi-");
  const auto proton     = pTable->FindParticle("proton");
  const auto neutron    = pTable->FindParticle("neutron");
  const auto lambda     = pTable->FindParticle("lambda");
  const auto sigma0     = pTable->FindParticle("sigma0");
  const G4double M_Kaon     = kaon_minus->GetPDGMass();
  const G4double M_PiM      = pi_minus->GetPDGMass();
  const G4double M_Proton   = proton->GetPDGMass();
  const G4double M_Neutron  = neutron->GetPDGMass();
  const G4double M_Lambda   = lambda->GetPDGMass();
  const G4double M_Sigma0   = sigma0->GetPDGMass();
  
  // loop until proper event is generated
  while(true){
    //--- Generate Fermi momentum for the target "n" ---
    TVector3 p_fermi_vec = FermiMotion::GetMomentum();
    G4ThreeVector target_mom_lab(p_fermi_vec.X(), p_fermi_vec.Y(), p_fermi_vec.Z());
    G4double E_target = sqrt(M_Neutron*M_Neutron + target_mom_lab.mag2());
    
    //========== Kinematics ==========
    //--- 1. K- + "n" -> Σ0 + pi- ---
    G4double beam_mom_mean = 1.4 * CLHEP::GeV;
    G4double beam_mom_sigma = beam_mom_mean * 3.3e-4 / 2.355;
    G4double beam_mom = G4RandGauss::shoot(beam_mom_mean, beam_mom_sigma);
    G4LorentzVector beam_lv(0, 0, beam_mom, sqrt(beam_mom*beam_mom + M_Kaon*M_Kaon));
    G4LorentzVector target_lv(target_mom_lab, E_target);
    G4LorentzVector W = beam_lv + target_lv;
    G4ThreeVector beta_cm = W.boostVector();
    G4double p_cm = CMSMomentum::GetCMSMomentum(W.m(), M_Sigma0, M_PiM);
    if (p_cm == 0.0) continue; 
    G4ThreeVector p_pi_cm_dir = G4RandomDirection();
    G4LorentzVector pi_lv_cm(p_pi_cm_dir * p_cm, sqrt(p_cm*p_cm + M_PiM*M_PiM));
    G4LorentzVector sigma0_lv_cm(-p_pi_cm_dir * p_cm, sqrt(p_cm*p_cm + M_Sigma0*M_Sigma0));
    G4LorentzVector pi_lv_lab = pi_lv_cm;
    pi_lv_lab.boost(beta_cm);
    G4LorentzVector sigma0_lv_lab = sigma0_lv_cm;
    sigma0_lv_lab.boost(beta_cm);
    if (pi_lv_lab.theta() > 15. * CLHEP::deg) continue; // S2s acceptance cut
    
    //--- 2. Σ0 -> Λ + γ ---
    G4ThreeVector beta_sigma0 = sigma0_lv_lab.boostVector();
    G4double p_sigma0_rest = CMSMomentum::GetCMSMomentum(M_Sigma0, M_Lambda, 0);
    G4ThreeVector p_lambda_rest_dir = G4RandomDirection();
    G4LorentzVector lambda_lv_rest(p_lambda_rest_dir * p_sigma0_rest, sqrt(p_sigma0_rest*p_sigma0_rest + M_Lambda*M_Lambda));
    G4LorentzVector lambda_lv_lab = lambda_lv_rest;
    lambda_lv_lab.boost(beta_sigma0);
    
    //--- 3. Λ -> p + pi- ---
    G4ThreeVector beta_lambda = lambda_lv_lab.boostVector();
    G4double p_lambda_rest = CMSMomentum::GetCMSMomentum(M_Lambda, M_Proton, M_PiM);
    G4ThreeVector p_decay_rest_dir = G4RandomDirection();
    G4LorentzVector p_decay_lv_rest(p_decay_rest_dir * p_lambda_rest, sqrt(p_lambda_rest*p_lambda_rest + M_Proton*M_Proton));
    G4LorentzVector pi_decay_lv_rest(-p_decay_rest_dir * p_lambda_rest, sqrt(p_lambda_rest*p_lambda_rest + M_PiM*M_PiM));
    G4LorentzVector p_decay_lv_lab = p_decay_lv_rest;
    p_decay_lv_lab.boost(beta_lambda);
    G4LorentzVector pi_decay_lv_lab = pi_decay_lv_rest;
    pi_decay_lv_lab.boost(beta_lambda);
    
    //========== Calculate Vertex position ==========
    static const auto& target_pos = geomMan.GetGlobalPosition("Target")*mm;
    static const auto& target_size = sizeMan.GetSize("Target");
    G4double target_radius = target_size.y() / 2.0 * mm;
    G4double target_height_half = target_size.z() / 2.0 * mm;
    G4double r_vtx_x, r_vtx_y, r_vtx_z;
    G4double beam_x_sigma = 0.2 * CLHEP::mm;
    G4double beam_y_sigma = 0.2 * CLHEP::mm;
    while (true) {
        r_vtx_x = G4RandGauss::shoot(0, beam_x_sigma);
        r_vtx_y = G4RandGauss::shoot(0, beam_y_sigma);
        r_vtx_z = G4RandFlat::shoot(-target_radius, target_radius);
        if (r_vtx_x*r_vtx_x + r_vtx_z*r_vtx_z < target_radius*target_radius && r_vtx_y < target_height_half && r_vtx_y > -target_height_half) break;
    }
    G4ThreeVector primary_vertex_pos = target_pos + G4ThreeVector(r_vtx_x, r_vtx_y, r_vtx_z);

    //--- Λ decay vertex (Σ0 lifetime is extremely short) ---
    G4double lambda_lifetime = lambda->GetPDGLifeTime(); 
    G4double decay_time_rest = -lambda_lifetime * log(G4RandFlat::shoot());
    G4double decay_time_lab = lambda_lv_lab.gamma() * decay_time_rest;
    G4ThreeVector velocity_vec = lambda_lv_lab.vect() * (CLHEP::c_light * CLHEP::c_light / lambda_lv_lab.e());
    G4ThreeVector flight_distance = velocity_vec * decay_time_lab;
    G4ThreeVector secondary_vertex_pos = primary_vertex_pos + flight_distance;

    //========== Particle Generation ==========
    //--- scattering pion  ---
    m_particleGun->SetParticleDefinition(pi_minus);
    m_particleGun->SetParticleMomentumDirection(pi_lv_lab.vect().unit());
    m_particleGun->SetParticleEnergy(pi_lv_lab.e() - M_PiM);
    m_particleGun->SetParticlePosition(primary_vertex_pos);
    m_particleGun->GeneratePrimaryVertex(anEvent);
    
    //--- spectator proton   ---
    G4ThreeVector total_mom_before = beam_lv.vect() + target_lv.vect();
    G4ThreeVector total_mom_after_reaction = pi_lv_lab.vect() + sigma0_lv_lab.vect();
    G4ThreeVector spectator_mom = total_mom_before - total_mom_after_reaction;

    G4double E_spectator = sqrt(M_Proton*M_Proton + spectator_mom.mag2());
    m_particleGun->SetParticleDefinition(proton);
    m_particleGun->SetParticleMomentumDirection(spectator_mom.unit());
    m_particleGun->SetParticleEnergy(E_spectator - M_Proton);
    m_particleGun->SetParticlePosition(primary_vertex_pos);
    m_particleGun->GeneratePrimaryVertex(anEvent);
    
    //--- Λ decay proton ---
    m_particleGun->SetParticleDefinition(proton);
    m_particleGun->SetParticleMomentumDirection(p_decay_lv_lab.vect().unit());
    m_particleGun->SetParticleEnergy(p_decay_lv_lab.e() - M_Proton);
    m_particleGun->SetParticlePosition(secondary_vertex_pos);
    m_particleGun->GeneratePrimaryVertex(anEvent);

    //--- Λ decay pion ---
    m_particleGun->SetParticleDefinition(pi_minus);
    m_particleGun->SetParticleMomentumDirection(pi_decay_lv_lab.vect().unit());
    m_particleGun->SetParticleEnergy(pi_decay_lv_lab.e() - M_PiM);
    m_particleGun->SetParticlePosition(secondary_vertex_pos);
    m_particleGun->GeneratePrimaryVertex(anEvent);
    
    break; 
  }
}

//_____________________________________________________________________________
void // 9004 E90 QF Sigma+ reaction
S2SPrimaryGeneratorAction::GenerateQFSigmaP(G4Event* anEvent)
{
  static const G4int n_particle = 1;
  m_particleGun = new G4ParticleGun(n_particle);

  //========== Define particles ==========
  const auto pTable = G4ParticleTable::GetParticleTable();
  const auto kaon_minus = pTable->FindParticle("kaon-");
  const auto pi_minus   = pTable->FindParticle("pi-");
  const auto sigma_plus = pTable->FindParticle("sigma+");
  const auto proton     = pTable->FindParticle("proton");
  const auto neutron    = pTable->FindParticle("neutron");

  const G4double M_Kaon     = kaon_minus->GetPDGMass();
  const G4double M_PiM      = pi_minus->GetPDGMass();
  const G4double M_SigmaP   = sigma_plus->GetPDGMass();
  const G4double M_Proton   = proton->GetPDGMass();
  const G4double M_Neutron  = neutron->GetPDGMass();

  // loop until proper event is generated
  while(true){

    //--- Generate Fermi momentum for the target "p" ---
    TVector3 p_fermi_vec = FermiMotion::GetMomentum();
    G4ThreeVector target_mom_lab(p_fermi_vec.X(), p_fermi_vec.Y(), p_fermi_vec.Z());
    G4double E_target = sqrt(M_Proton*M_Proton + target_mom_lab.mag2());
    
    //========== Kinematics ==========
    //--- 1. K- + "p" -> Σ+ + pi- ---
    G4double beam_mom_mean = 1.4 * CLHEP::GeV;
    G4double beam_mom_sigma = beam_mom_mean * 3.3e-4 / 2.355;
    G4double beam_mom = G4RandGauss::shoot(beam_mom_mean, beam_mom_sigma);
    G4LorentzVector beam_lv(0, 0, beam_mom, sqrt(beam_mom*beam_mom + M_Kaon*M_Kaon));
    G4LorentzVector target_lv(target_mom_lab, E_target); // Target is a proton
    G4LorentzVector W = beam_lv + target_lv;
    G4ThreeVector beta_cm = W.boostVector();
    G4double p_cm = CMSMomentum::GetCMSMomentum(W.m(), M_SigmaP, M_PiM);
    if (p_cm == 0.0) continue; 
    G4ThreeVector p_pi_cm_dir = G4RandomDirection();
    G4LorentzVector pi_lv_cm(p_pi_cm_dir * p_cm, sqrt(p_cm*p_cm + M_PiM*M_PiM));
    G4LorentzVector sigma_plus_lv_cm(-p_pi_cm_dir * p_cm, sqrt(p_cm*p_cm + M_SigmaP*M_SigmaP));
    G4LorentzVector pi_lv_lab = pi_lv_cm;
    pi_lv_lab.boost(beta_cm);
    G4LorentzVector sigma_plus_lv_lab = sigma_plus_lv_cm;
    sigma_plus_lv_lab.boost(beta_cm);
    if (pi_lv_lab.theta() > 15. * CLHEP::deg) continue; // S2s acceptance cut
    
    //--- 2. Σ+ -> n + π+ OR p + π0 ---
    G4ThreeVector beta_sigma_plus = sigma_plus_lv_lab.boostVector();
    G4double p_sigma_plus_rest_npi_plus = CMSMomentum::GetCMSMomentum(M_SigmaP, M_Neutron, pTable->FindParticle("pi+")->GetPDGMass());
    G4double p_sigma_plus_rest_ppi0 = CMSMomentum::GetCMSMomentum(M_SigmaP, M_Proton, pTable->FindParticle("pi0")->GetPDGMass());

    G4int decay_mode = (G4RandFlat::shoot() < 0.5) ? 1 : 2;
    
    G4ThreeVector primary_vertex_pos; // This needs to be defined
    
    if (decay_mode == 1) {
      G4ThreeVector n_decay_rest_dir = G4RandomDirection();
      G4LorentzVector n_decay_lv_rest(n_decay_rest_dir * p_sigma_plus_rest_npi_plus, sqrt(p_sigma_plus_rest_npi_plus*p_sigma_plus_rest_npi_plus + M_Neutron*M_Neutron));
      G4LorentzVector pi_plus_decay_lv_rest(-n_decay_rest_dir * p_sigma_plus_rest_npi_plus, sqrt(p_sigma_plus_rest_npi_plus*p_sigma_plus_rest_npi_plus + pTable->FindParticle("pi+")->GetPDGMass()*pTable->FindParticle("pi+")->GetPDGMass()));

      G4LorentzVector pi_plus_decay_lv_lab = pi_plus_decay_lv_rest;
      pi_plus_decay_lv_lab.boost(beta_sigma_plus);
      
      m_particleGun->SetParticleDefinition(pTable->FindParticle("pi+"));
      m_particleGun->SetParticleMomentumDirection(pi_plus_decay_lv_lab.vect().unit());
      m_particleGun->SetParticleEnergy(pi_plus_decay_lv_lab.e() - pTable->FindParticle("pi+")->GetPDGMass());
      m_particleGun->SetParticlePosition(primary_vertex_pos);
      m_particleGun->GeneratePrimaryVertex(anEvent);
    } else {
      G4ThreeVector p_decay_rest_dir = G4RandomDirection();
      G4LorentzVector p_decay_lv_rest(p_decay_rest_dir * p_sigma_plus_rest_ppi0, sqrt(p_sigma_plus_rest_ppi0*p_sigma_plus_rest_ppi0 + M_Proton*M_Proton));
      G4LorentzVector p_decay_lv_lab = p_decay_lv_rest;
      p_decay_lv_lab.boost(beta_sigma_plus);
      
      m_particleGun->SetParticleDefinition(pTable->FindParticle("proton"));
      m_particleGun->SetParticleMomentumDirection(p_decay_lv_lab.vect().unit());
      m_particleGun->SetParticleEnergy(p_decay_lv_lab.e() - M_Proton);
      m_particleGun->SetParticlePosition(primary_vertex_pos);
      m_particleGun->GeneratePrimaryVertex(anEvent);
    }

    //========== Calculate Vertex position ==========
    static const auto& target_pos = geomMan.GetGlobalPosition("Target")*mm;
    static const auto& target_size = sizeMan.GetSize("Target");
    G4double target_radius = target_size.y() / 2.0 * mm;
    G4double target_height_half = target_size.z() / 2.0 * mm;
    G4double r_vtx_x, r_vtx_y, r_vtx_z;
    G4double beam_x_sigma = 0.2 * CLHEP::mm;
    G4double beam_y_sigma = 0.2 * CLHEP::mm;
    while (true) {
        r_vtx_x = G4RandGauss::shoot(0, beam_x_sigma);
        r_vtx_y = G4RandGauss::shoot(0, beam_y_sigma);
        r_vtx_z = G4RandFlat::shoot(-target_radius, target_radius);
        if (r_vtx_x*r_vtx_x + r_vtx_z*r_vtx_z < target_radius*target_radius && r_vtx_y < target_height_half && r_vtx_y > -target_height_half) break;
    }
    primary_vertex_pos = target_pos + G4ThreeVector(r_vtx_x, r_vtx_y, r_vtx_z);

    //========== Particle Generation ==========
    //--- scattering pion  ---
    m_particleGun->SetParticleDefinition(pi_minus);
    m_particleGun->SetParticleMomentumDirection(pi_lv_lab.vect().unit());
    m_particleGun->SetParticleEnergy(pi_lv_lab.e() - M_PiM);
    m_particleGun->SetParticlePosition(primary_vertex_pos);
    m_particleGun->GeneratePrimaryVertex(anEvent);
    
    //--- spectator neutron   ---
    G4ThreeVector total_mom_before = beam_lv.vect() + target_lv.vect();
    G4ThreeVector total_mom_after_reaction = pi_lv_lab.vect() + sigma_plus_lv_lab.vect();
    G4ThreeVector spectator_mom = total_mom_before - total_mom_after_reaction;
    
    G4double E_spectator = sqrt(M_Neutron*M_Neutron + spectator_mom.mag2());
    m_particleGun->SetParticleDefinition(neutron);
    m_particleGun->SetParticleMomentumDirection(spectator_mom.unit());
    m_particleGun->SetParticleEnergy(E_spectator - M_Neutron);
    m_particleGun->SetParticlePosition(primary_vertex_pos);
    m_particleGun->GeneratePrimaryVertex(anEvent);
    
    break; 
  }
}