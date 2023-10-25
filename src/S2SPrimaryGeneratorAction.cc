// -*- C++ -*-

#include "S2SPrimaryGeneratorAction.hh"

#include <G4Event.hh>
#include <G4ParticleGun.hh>
#include <G4ParticleTable.hh>
#include <G4ParticleDefinition.hh>
#include <G4SystemOfUnits.hh>
#include <G4LorentzVector.hh>
#include <G4ThreeVector.hh>
#include <Randomize.hh>

#include <TMath.h>

#include "ConfMan.hh"
#include "DCGeomMan.hh"
#include "FuncName.hh"
#include "MagnetConstant.hh"
#include "S2SAnaManager.hh"
#include "DetSizeMan.hh"

namespace
{
const auto& confMan = ConfMan::GetInstance();
const auto& geomMan = DCGeomMan::GetInstance();
const auto& sizeMan = DetSizeMan::GetInstance();
auto& anaMan = S2SAnaManager::GetInstance();
const auto particleTable = G4ParticleTable::GetParticleTable();
}

//_____________________________________________________________________________
S2SPrimaryGeneratorAction::S2SPrimaryGeneratorAction()
  : G4VUserPrimaryGeneratorAction(),
    T2Distance(),
    momcent(),
    mombite(),
    GenPID(1),
    generator(confMan.Get<G4int>("Generator"))
{
}

//_____________________________________________________________________________
S2SPrimaryGeneratorAction::~S2SPrimaryGeneratorAction()
{
  if(particleGun) delete particleGun;
}

//_____________________________________________________________________________
void
S2SPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  switch(generator){
  case 0: break;//no generation
  case 1: GenerateUniform0(anEvent); break; //Uniform generation
  case 2: GenerateFocusCheck(anEvent); break;//Mom=0, +-97.5MeV/c (0, 2, 4) KanatsukiM Fig.4.1
  case 3: GenerateMonoBeam(anEvent); break;
  case 4: GenerateAcceptance(anEvent); break;
  default:
    G4cerr << " * Generator number error : " << generator << G4endl;
    break;
  }
}

//_____________________________________________________________________________
void
S2SPrimaryGeneratorAction::GenerateUniform0(G4Event* anEvent)
{
  G4int n_particle = 1;
  particleGun = new G4ParticleGun(n_particle);
  if(GenPID==1){
    m_particle = particleTable->FindParticle("kaon+");
  }
  else if(GenPID==2){
    m_particle = particleTable->FindParticle("kaon-");
  }
  else if(GenPID==3){
    m_particle = particleTable->FindParticle("pi+");
  }
  else if(GenPID==4){
    m_particle = particleTable->FindParticle("pi-");
  }
  else if(GenPID==5){
    m_particle = particleTable->FindParticle("proton");
  }
  else if(GenPID==6){
    m_particle = particleTable->FindParticle("e-");
  }
  else if(GenPID==7){
    m_particle = particleTable->FindParticle("mu-");
  }
  else if(GenPID==8){
    m_particle = particleTable->FindParticle("xi-");
  }
  else {
    G4cout << " Sorry, I do not know GenPID=" << GenPID << G4endl;
    G4cout << " --> Set particle type: kaon+" << G4endl;
    m_particle = particleTable->FindParticle("kaon+");
  }
  particleGun -> SetParticleDefinition(m_particle);

  // ~~~~~~~~~ Set Particle Momentum ~~~~~~~~~~~~~~~~
  G4double parMom = 0.0;
  parMom = (momcent+ mombite*(G4UniformRand()-0.5)*2.0 ) * GeV;
  particleGun->SetParticleMomentum(parMom);

  //~~~~~~~~~~ Set Particle Position ~~~~~~~~~~~~~~~
  //G4double xlim = 10*cm;
  //G4double ylim = 2.5*cm;
  //G4double zlim = 1.5*cm;
  G4double xorg = -rhoD*tan(bendAngleD*TMath::DegToRad()/2.)-driftL2-Q2z-driftL1-Q1z;
  G4double x0 = xorg-T2Distance; // Beam direction
  G4double y0 = 0.0; // Horizontal direction
  G4double z0 = 0.0; // Vertical direction
  x0 = x0 + beamz*( G4UniformRand()-0.5 ); // Beam direction
  y0 = CLHEP::RandGauss::shoot(y0,beamx); // Horizontal direction
  z0 = CLHEP::RandGauss::shoot(z0,beamy); // Vertical direction
  //   x0 += xlim*int(G4UniformRand()*3);//mm From -T2Distance to -T2Distance+xlim
  //   y0 = 2*ylim*(G4UniformRand()-0.5);//mm
  //   z0 = 2*zlim*(G4UniformRand()-0.5);//mm
  //   while(1){
  //     y0 = 2*ylim*(G4UniformRand()-0.5);//mm
  //     z0 = 2*zlim*(G4UniformRand()-0.5);//mm
  //     //     y0 = CLHEP::RandGauss::shoot(0,2.0*cm);
  //     //     z0 = CLHEP::RandGauss::shoot(0,0.5*cm);
  //     if(fabs(y0)<ylim && fabs(z0)<zlim){
  //       break;
  //     }
  //   }
  G4ThreeVector vertex(x0,y0,z0);
  particleGun->SetParticlePosition(vertex);

  // ~~~~~~~~~~~ Particle Direction ~~~~~~~~~~~~~~~~~~
  G4double limitTheta = thetamax*TMath::DegToRad(); // [deg] --> [rad]
  //G4double limitTheta = 25*TMath::DegToRad(); //degree
  //G4double limituAng = 10*TMath::DegToRad();
  //G4double limitvAng = 20*TMath::DegToRad();
  G4double limitcos = cos(limitTheta);
  G4double phi   = 0.0;
  G4double theta = 0.0;
  G4double rand  = 0.0;
  G4double u0    = 0.0;
  G4double v0    = 0.0;
  phi = 2.0 * TMath::Pi() * G4UniformRand();
  rand = (1-limitcos)*G4UniformRand()+limitcos;
  theta = acos(rand);
  u0 = tan(theta)*cos(phi);
  v0 = tan(theta)*sin(phi);
  //  while(1){
  //    phi = 2*TMath::Pi()*G4UniformRand();
  //    rand = (1-limitcos)*G4UniformRand()+limitcos;
  //    theta = acos(rand);
  //    u0 = tan(theta)*cos(phi);
  //    v0 = tan(theta)*sin(phi);
  //    if(1){
  //    //if(fabs(atan(u0))<limituAng && fabs(atan(v0))<limitvAng){
  //      break;
  //    }
  //  }
  //   if(anEvent->GetEventID()/3==0)
  //     u0=-4*TMath::DegToRad();
  //   if(anEvent->GetEventID()/3==1)
  //     u0=-2*TMath::DegToRad();
  //   if(anEvent->GetEventID()/3==2)
  //     u0=-0*TMath::DegToRad();
  //   if(anEvent->GetEventID()/3==3)
  //     u0=1.8*TMath::DegToRad();
  //   if(anEvent->GetEventID()/3==4)
  //     u0=3.6*TMath::DegToRad();
  //   //(((anEvent->GetEventID())/3) -2) *TMath::DegToRad();
  //   v0=0;
  G4ThreeVector direction(1., u0, v0);
  particleGun->SetParticleMomentumDirection(direction);


  //particleGun->SetParticleEnergy(parKinE*MeV);
  //double parMom particleGun->GetParticleMomentum();
  double parKinE = particleGun->GetParticleEnergy();
  particleGun-> GeneratePrimaryVertex(anEvent);

  anaMan.SetPrimaryData(y0/mm,z0/mm,-x0/mm+xorg,
                        u0,v0,
                        phi,theta,
                        parMom/GeV,
                        parKinE/GeV,
                        9999);
  G4cout << FUNC_NAME << " gun." << G4endl;
}

//_____________________________________________________________________________
void
S2SPrimaryGeneratorAction::GenerateFocusCheck(G4Event* anEvent)
{
 G4int n_particle = 1;
  particleGun = new G4ParticleGun(n_particle);
  particleGun->SetParticleEnergy(0.0*GeV); // Reset Particle energy

  // ~~~~~~~~~~~ Set Particle Type ~~~~~~~~~~~~~~~~
  if(GenPID==1){
    m_particle = particleTable->FindParticle("kaon+");
  }
  else if(GenPID==2){
    m_particle = particleTable->FindParticle("kaon-");
  }
  else if(GenPID==3){
    m_particle = particleTable->FindParticle("pi+");
  }
  else if(GenPID==4){
    m_particle = particleTable->FindParticle("pi-");
  }
  else if(GenPID==5){
    m_particle = particleTable->FindParticle("proton");
  }
  else if(GenPID==6){
    m_particle = particleTable->FindParticle("e-");
  }
  else if(GenPID==7){
    m_particle = particleTable->FindParticle("mu-");
  }
  else if(GenPID==8){
    m_particle = particleTable->FindParticle("xi-");
  }
  else {
    G4cout << " Sorry, I do not know GenPID=" << GenPID << G4endl;
    G4cout << " --> Set particle type: kaon+" << G4endl;
    m_particle = particleTable->FindParticle("kaon+");
  }
  particleGun -> SetParticleDefinition( m_particle );


  // ~~~~~~~~~ Set Particle Momentum ~~~~~~~~~~~~~~~~
  G4double parMom = 0.0;
  //parMom = (momcent+ mombite*(G4UniformRand()-0.5)*2.0 ) * GeV;
  parMom = (momcent + (double)((int)(G4UniformRand()*3.)-1)*mombite)*GeV;
  std::cout<<"parMom="<<parMom<<std::endl;
  particleGun->SetParticleMomentum(parMom);


  //~~~~~~~~~~ Set Particle Position ~~~~~~~~~~~~~~~
  //G4double xlim = 10*cm;
  //G4double ylim = 2.5*cm;
  //G4double zlim = 1.5*cm;
  G4double xorg = -rhoD*tan(bendAngleD*TMath::DegToRad()/2.)-driftL2-Q2z-driftL1-Q1z;
  G4double x0 = xorg-T2Distance; // Beam direction
  G4double y0 = 0.0; // Horizontal direction
  G4double z0 = 0.0; // Vertical direction
  x0 = x0 + beamz*( G4UniformRand()-0.5 ); // Beam direction
  y0 = CLHEP::RandGauss::shoot(y0,beamx); // Horizontal direction
  z0 = CLHEP::RandGauss::shoot(z0,beamy); // Vertical direction
  //   x0 += xlim*int(G4UniformRand()*3);//mm From -T2Distance to -T2Distance+xlim
  //   y0 = 2*ylim*(G4UniformRand()-0.5);//mm
  //   z0 = 2*zlim*(G4UniformRand()-0.5);//mm
  //   while(1){
  //     y0 = 2*ylim*(G4UniformRand()-0.5);//mm
  //     z0 = 2*zlim*(G4UniformRand()-0.5);//mm
  //     //     y0 = CLHEP::RandGauss::shoot(0,2.0*cm);
  //     //     z0 = CLHEP::RandGauss::shoot(0,0.5*cm);
  //     if(fabs(y0)<ylim && fabs(z0)<zlim){
  //       break;
  //     }
  //   }
  G4ThreeVector vertex(x0,y0,z0);
  particleGun->SetParticlePosition(vertex);

  // ~~~~~~~~~~~ Particle Direction ~~~~~~~~~~~~~~~~~~
  //G4double limitTheta = thetamax*TMath::DegToRad(); // [deg] --> [rad]
  //G4double limitTheta = 25*TMath::DegToRad(); //degree
  //G4double limituAng = 10*TMath::DegToRad();
  //G4double limitvAng = 20*TMath::DegToRad();
  //G4double limitcos = cos(limitTheta);
  G4double phi   = 0.0;
  G4double theta = 0.0;
  // G4double rand  = 0.0;
  G4double u0    = 0.0;
  G4double v0    = 0.0;
  //  phi = 2.0 * TMath::Pi() * G4UniformRand();
  phi = (double)((int)(G4UniformRand()*4.))*0.5*TMath::Pi();
  //rand = (1-limitcos)*G4UniformRand()+limitcos;
  //theta = acos(rand);
  theta = (double)((int)(G4UniformRand()*3.)*2)*TMath::Pi()/180.;
  std::cout<<"theta="<<theta<<"rad, "<<theta*180./TMath::Pi()<<"degree"<<std::endl;
  std::cout<<"phi="<<phi<<"rad, "<<phi*180./TMath::Pi()<<"degree"<<std::endl;
  u0 = tan(theta)*cos(phi);
  v0 = tan(theta)*sin(phi);
  //  while(1){
  //    phi = 2*TMath::Pi()*G4UniformRand();
  //    rand = (1-limitcos)*G4UniformRand()+limitcos;
  //    theta = acos(rand);
  //    u0 = tan(theta)*cos(phi);
  //    v0 = tan(theta)*sin(phi);
  //    if(1){
  //    //if(fabs(atan(u0))<limituAng && fabs(atan(v0))<limitvAng){
  //      break;
  //    }
  //  }
  //   if(anEvent->GetEventID()/3==0)
  //     u0=-4*TMath::DegToRad();
  //   if(anEvent->GetEventID()/3==1)
  //     u0=-2*TMath::DegToRad();
  //   if(anEvent->GetEventID()/3==2)
  //     u0=-0*TMath::DegToRad();
  //   if(anEvent->GetEventID()/3==3)
  //     u0=1.8*TMath::DegToRad();
  //   if(anEvent->GetEventID()/3==4)
  //     u0=3.6*TMath::DegToRad();
  //   //(((anEvent->GetEventID())/3) -2) *TMath::DegToRad();
  //   v0=0;
  G4ThreeVector direction(1., u0, v0);
  particleGun->SetParticleMomentumDirection(direction);


  //particleGun->SetParticleEnergy(parKinE*MeV);
  //double parMom particleGun->GetParticleMomentum();
  double parKinE = particleGun->GetParticleEnergy();
  particleGun-> GeneratePrimaryVertex(anEvent);

  anaMan.SetPrimaryData(y0/mm,z0/mm,-x0/mm+xorg,
                        u0,v0,
                        phi,theta,
                        parMom/GeV,
                        parKinE/GeV,
                        9999);
}

//_____________________________________________________________________________
void
S2SPrimaryGeneratorAction::GenerateMonoBeam(G4Event* anEvent)
{
  const G4int n_particle = 1;
  particleGun = new G4ParticleGun(n_particle);
  static const G4String name = "kaon+";
  // static const G4String name = "proton";
  m_particle = particleTable->FindParticle(name);
  static const auto pdg = m_particle->GetPDGEncoding();
  const G4double m0 = m_particle->GetPDGMass();
  const G4double p0 = 1.37*CLHEP::GeV;
  const auto& target_pos = geomMan.GetGlobalPosition("Target");
  G4LorentzVector p(0, 0, p0, TMath::Sqrt(p0*p0 + m0*m0));
  G4LorentzVector v(target_pos, 0);
  particleGun->SetParticleDefinition(m_particle);
  particleGun->SetParticleMomentumDirection(p.v());
  particleGun->SetParticleEnergy(p.e() - m0);
  particleGun->SetParticlePosition(v.v());
  particleGun->GeneratePrimaryVertex(anEvent);
  anaMan.SetPrimaryParticle(0, pdg, p, v);
}

//_____________________________________________________________________________
void
S2SPrimaryGeneratorAction::GenerateAcceptance(G4Event* anEvent)
{
  const G4int n_particle = 1;
  particleGun = new G4ParticleGun(n_particle);
  static const G4String name = "kaon+";
  // static const G4String name = "proton";
  m_particle = particleTable->FindParticle(name);
  static const auto pdg = m_particle->GetPDGEncoding();
  static const auto& target_pos = geomMan.GetGlobalPosition("Target")*CLHEP::mm;
  static const auto& target_size = sizeMan.GetSize("Target")*CLHEP::mm/2;
  static const G4double m0 = m_particle->GetPDGMass();
  G4double p0 = G4RandFlat::shoot(1.0, 1.8)*CLHEP::GeV;
  G4double theta = std::acos(G4RandFlat::shoot(std::sqrt(3.)/2., 1.))*CLHEP::radian;
  G4double phi = G4RandFlat::shoot(0., 180.)*CLHEP::deg;
  G4LorentzVector p(0, 0, 0, TMath::Sqrt(p0*p0 + m0*m0));
  p.setRThetaPhi(p0, theta, phi);
  G4double x0 =  G4RandFlat::shoot(-target_size.x(), target_size.x());
  G4double y0 =  G4RandFlat::shoot(-target_size.y(), target_size.y());
  G4double z0 =  G4RandFlat::shoot(-target_size.z(), target_size.z());
  G4LorentzVector v(target_pos + G4ThreeVector(x0, y0, z0), 0);
#if 0
  G4cout << FUNC_NAME << G4endl
         << " " << p0 << " " << theta/CLHEP::deg << " " << phi/CLHEP::deg << G4endl
         << " " << p << " " << p.theta()/CLHEP::deg << " " << v << G4endl;
#endif
  particleGun->SetParticleDefinition(m_particle);
  particleGun->SetParticleMomentumDirection(p.v());
  particleGun->SetParticleEnergy(p.e() - m0);
  particleGun->SetParticlePosition(v.v());
  particleGun->GeneratePrimaryVertex(anEvent);
  anaMan.SetPrimaryParticle(0, pdg, p, v);
}
