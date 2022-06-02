/*
  "S2SPrimaryGeneratorAction.cc"

  Modified by Toshi Gogami , 10Nov2014
*/

#include "S2SPrimaryGeneratorAction.hh"
#include "RadDeg.hh"
#include "Area.hh"
#include "S2SAnalysis.hh"
#include "MagnetConstant.hh"
#include "ConfMan.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "Randomize.hh"
#include "G4SystemOfUnits.hh"
using namespace CLHEP;
 

S2SPrimaryGeneratorAction::S2SPrimaryGeneratorAction(S2SAnalysis *ana,
						     G4double parMomPre_,
						     G4double T2distance_)
  :  anaMan_(ana),parMomPre(parMomPre_),T2Distance(T2distance_),
     momcent(0.0),mombite(0), GenPID(1), generator(0)
{
  //G4int n_particle = 1;
  //default particle kinematics
  //particleGun= new G4ParticleGun(n_particle);
  //particleGun->SetParticleMomentumDirection(G4ThreeVector(1., 0., 0.));
  //particleGun->SetParticlePosition(G4ThreeVector(-4.92*m, 0., 0.));
  //particleGun->SetParticleMomentum(parMomPre);
}

//S2SPrimaryGeneratorAction::S2SPrimaryGeneratorAction(S2SAnalysis *ana,
//						     ConfMan* confman_,
//						     G4double T2distance_)
//  :  anaMan_(ana),parMomPre(1300.0),T2Distance(T2distance_),
//     confman(confman_),momcent(0.0),mombite(0),thetamax(30.0),
//     GenPID(1)
//{
//  momcent  = confman->GetGenMomCent();
//  mombite  = confman->GetGenMomBite();
//  thetamax = confman->GetGenTheta();
//  GenPID   = confman->GetGenPID();
//}

S2SPrimaryGeneratorAction::S2SPrimaryGeneratorAction(S2SAnalysis *ana, ConfMan* confman_)
  :  anaMan_(ana),parMomPre(1300.0),T2Distance(600.0*mm),
     confman(confman_),momcent(0.0),mombite(0),thetamax(30.0),
     GenPID(1), generator(0)
{
  momcent  = confman->GetGenMomCent();
  mombite  = confman->GetGenMomBite();
  thetamax = confman->GetGenTheta();
  GenPID   = confman->GetGenPID();
  generator   = confman->GetGenerator();
  T2Distance = confman->GetTargetPositionZ() * mm;
  beamx    = confman->GetBeamX() *mm;
  beamy    = confman->GetBeamY() *mm;
  beamz    = confman->GetBeamZ() *mm;
}

S2SPrimaryGeneratorAction::~S2SPrimaryGeneratorAction()
{
  delete particleGun;
  //  delete particle;
}

void S2SPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  switch(generator){
  case 0: break;//no generation
  case 1: GenerateUniform0(anEvent); break;//Uniform generation
  case 2: GenerateFocusCheck(anEvent); break;//Mom=0, +-97.5MeV/c (0, 2, 4) KanatsukiM Fig.4.1
  default:
    G4cerr << " * Generator number error : " << generator << G4endl;
    break;
  }
}

void S2SPrimaryGeneratorAction::GenerateUniform0(G4Event* anEvent)
{
 G4int n_particle = 1;
  particleGun = new G4ParticleGun(n_particle);
  particleGun->SetParticleEnergy(0.0*GeV); // Reset Particle energy
  
  // ~~~~~~~~~~~ Set Particle Type ~~~~~~~~~~~~~~~~
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  if(GenPID==1){
    particle = particleTable->FindParticle("kaon+");
  }
  else if(GenPID==2){ 
    particle = particleTable->FindParticle("kaon-");
  }
  else if(GenPID==3){
    particle = particleTable->FindParticle("pi+");
  }
  else if(GenPID==4){
    particle = particleTable->FindParticle("pi-");
  }
  else if(GenPID==5){
    particle = particleTable->FindParticle("proton");
  }
  else if(GenPID==6){
    particle = particleTable->FindParticle("e-");
  }
  else if(GenPID==7){
    particle = particleTable->FindParticle("mu-");
  }
  else if(GenPID==8){
    particle = particleTable->FindParticle("xi-");
  }
  else {
    G4cout << " Sorry, I do not know GenPID=" << GenPID << G4endl;
    G4cout << " --> Set particle type: kaon+" << G4endl;
    particle = particleTable->FindParticle("kaon+");
  }
  particleGun -> SetParticleDefinition( particle );
  
  
  // ~~~~~~~~~ Set Particle Momentum ~~~~~~~~~~~~~~~~
  G4double parMom = 0.0;
  parMom = (momcent+ mombite*(G4UniformRand()-0.5)*2.0 ) * GeV;
  particleGun->SetParticleMomentum(parMom);

  
  //~~~~~~~~~~ Set Particle Position ~~~~~~~~~~~~~~~
  //G4double xlim = 10*cm;
  //G4double ylim = 2.5*cm;
  //G4double zlim = 1.5*cm;
  G4double xorg = -rhoD*tan(bendAngleD*Deg2Rad/2.)-driftL2-Q2z-driftL1-Q1z;
  G4double x0 = xorg-T2Distance; // Beam direction
  G4double y0 = 0.0; // Horizontal direction
  G4double z0 = 0.0; // Vertical direction
  x0 = x0 + beamz*( G4UniformRand()-0.5 ); // Beam direction
  y0 = RandGauss::shoot(y0,beamx); // Horizontal direction
  z0 = RandGauss::shoot(z0,beamy); // Vertical direction
  //   x0 += xlim*int(G4UniformRand()*3);//mm From -T2Distance to -T2Distance+xlim
  //   y0 = 2*ylim*(G4UniformRand()-0.5);//mm
  //   z0 = 2*zlim*(G4UniformRand()-0.5);//mm
  //   while(1){
  //     y0 = 2*ylim*(G4UniformRand()-0.5);//mm
  //     z0 = 2*zlim*(G4UniformRand()-0.5);//mm
  //     //     y0 = RandGauss::shoot(0,2.0*cm);
  //     //     z0 = RandGauss::shoot(0,0.5*cm);
  //     if(fabs(y0)<ylim && fabs(z0)<zlim){
  //       break;    
  //     }
  //   }
  G4ThreeVector vertex(x0,y0,z0);
  particleGun->SetParticlePosition(vertex);

  // ~~~~~~~~~~~ Particle Direction ~~~~~~~~~~~~~~~~~~
  G4double limitTheta = thetamax*Deg2Rad; // [deg] --> [rad]
  //G4double limitTheta = 25*Deg2Rad; //degree
  //G4double limituAng = 10*Deg2Rad;
  //G4double limitvAng = 20*Deg2Rad;
  G4double limitcos = cos(limitTheta);
  G4double phi   = 0.0;
  G4double theta = 0.0;
  G4double rand  = 0.0;
  G4double u0    = 0.0;
  G4double v0    = 0.0;
  phi = 2.0 * M_PI * G4UniformRand();
  rand = (1-limitcos)*G4UniformRand()+limitcos;
  theta = acos(rand);
  u0 = tan(theta)*cos(phi);
  v0 = tan(theta)*sin(phi);
  //  while(1){
  //    phi = 2*M_PI*G4UniformRand();
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
  //     u0=-4*Deg2Rad;
  //   if(anEvent->GetEventID()/3==1) 
  //     u0=-2*Deg2Rad;
  //   if(anEvent->GetEventID()/3==2) 
  //     u0=-0*Deg2Rad;
  //   if(anEvent->GetEventID()/3==3) 
  //     u0=1.8*Deg2Rad;
  //   if(anEvent->GetEventID()/3==4) 
  //     u0=3.6*Deg2Rad;
  //   //(((anEvent->GetEventID())/3) -2) *Deg2Rad;
  //   v0=0;
  G4ThreeVector direction(1., u0, v0);
  particleGun->SetParticleMomentumDirection(direction);
  
  
  //particleGun->SetParticleEnergy(parKinE*MeV);
  //double parMom particleGun->GetParticleMomentum();
  double parKinE = particleGun->GetParticleEnergy();
  particleGun-> GeneratePrimaryVertex(anEvent);
  
  anaMan_->SetPrimaryData(y0/mm,z0/mm,-x0/mm+xorg,
			  u0,v0,
			  phi,theta,
			  parMom/GeV,
			  parKinE/GeV,
			  9999);
}

void S2SPrimaryGeneratorAction::GenerateFocusCheck(G4Event* anEvent)
{
 G4int n_particle = 1;
  particleGun = new G4ParticleGun(n_particle);
  particleGun->SetParticleEnergy(0.0*GeV); // Reset Particle energy
  
  // ~~~~~~~~~~~ Set Particle Type ~~~~~~~~~~~~~~~~
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  if(GenPID==1){
    particle = particleTable->FindParticle("kaon+");
  }
  else if(GenPID==2){ 
    particle = particleTable->FindParticle("kaon-");
  }
  else if(GenPID==3){
    particle = particleTable->FindParticle("pi+");
  }
  else if(GenPID==4){
    particle = particleTable->FindParticle("pi-");
  }
  else if(GenPID==5){
    particle = particleTable->FindParticle("proton");
  }
  else if(GenPID==6){
    particle = particleTable->FindParticle("e-");
  }
  else if(GenPID==7){
    particle = particleTable->FindParticle("mu-");
  }
  else if(GenPID==8){
    particle = particleTable->FindParticle("xi-");
  }
  else {
    G4cout << " Sorry, I do not know GenPID=" << GenPID << G4endl;
    G4cout << " --> Set particle type: kaon+" << G4endl;
    particle = particleTable->FindParticle("kaon+");
  }
  particleGun -> SetParticleDefinition( particle );
  
  
  // ~~~~~~~~~ Set Particle Momentum ~~~~~~~~~~~~~~~~
  G4double parMom = 0.0;
  //parMom = (momcent+ mombite*(G4UniformRand()-0.5)*2.0 ) * GeV;
  parMom = (momcent + (double)((int)(G4UniformRand()*3.)-1)*0.0975)*GeV; 
  std::cout<<"parMom="<<parMom<<std::endl;
  particleGun->SetParticleMomentum(parMom);

  
  //~~~~~~~~~~ Set Particle Position ~~~~~~~~~~~~~~~
  //G4double xlim = 10*cm;
  //G4double ylim = 2.5*cm;
  //G4double zlim = 1.5*cm;
  G4double xorg = -rhoD*tan(bendAngleD*Deg2Rad/2.)-driftL2-Q2z-driftL1-Q1z;
  G4double x0 = xorg-T2Distance; // Beam direction
  G4double y0 = 0.0; // Horizontal direction
  G4double z0 = 0.0; // Vertical direction
  x0 = x0 + beamz*( G4UniformRand()-0.5 ); // Beam direction
  y0 = RandGauss::shoot(y0,beamx); // Horizontal direction
  z0 = RandGauss::shoot(z0,beamy); // Vertical direction
  //   x0 += xlim*int(G4UniformRand()*3);//mm From -T2Distance to -T2Distance+xlim
  //   y0 = 2*ylim*(G4UniformRand()-0.5);//mm
  //   z0 = 2*zlim*(G4UniformRand()-0.5);//mm
  //   while(1){
  //     y0 = 2*ylim*(G4UniformRand()-0.5);//mm
  //     z0 = 2*zlim*(G4UniformRand()-0.5);//mm
  //     //     y0 = RandGauss::shoot(0,2.0*cm);
  //     //     z0 = RandGauss::shoot(0,0.5*cm);
  //     if(fabs(y0)<ylim && fabs(z0)<zlim){
  //       break;    
  //     }
  //   }
  G4ThreeVector vertex(x0,y0,z0);
  particleGun->SetParticlePosition(vertex);

  // ~~~~~~~~~~~ Particle Direction ~~~~~~~~~~~~~~~~~~
  //G4double limitTheta = thetamax*Deg2Rad; // [deg] --> [rad]
  //G4double limitTheta = 25*Deg2Rad; //degree
  //G4double limituAng = 10*Deg2Rad;
  //G4double limitvAng = 20*Deg2Rad;
  //G4double limitcos = cos(limitTheta);
  G4double phi   = 0.0;
  G4double theta = 0.0;
  G4double rand  = 0.0;
  G4double u0    = 0.0;
  G4double v0    = 0.0;
  //  phi = 2.0 * M_PI * G4UniformRand();
  phi = (double)((int)(G4UniformRand()*4.))*0.5*M_PI;
  //rand = (1-limitcos)*G4UniformRand()+limitcos;
  //theta = acos(rand);
  theta = (double)((int)(G4UniformRand()*3.)*2)*M_PI/180.;
  std::cout<<"theta="<<theta<<"rad, "<<theta*180./M_PI<<"degree"<<std::endl;
  std::cout<<"phi="<<phi<<"rad, "<<phi*180./M_PI<<"degree"<<std::endl;
  u0 = tan(theta)*cos(phi);
  v0 = tan(theta)*sin(phi);
  //  while(1){
  //    phi = 2*M_PI*G4UniformRand();
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
  //     u0=-4*Deg2Rad;
  //   if(anEvent->GetEventID()/3==1) 
  //     u0=-2*Deg2Rad;
  //   if(anEvent->GetEventID()/3==2) 
  //     u0=-0*Deg2Rad;
  //   if(anEvent->GetEventID()/3==3) 
  //     u0=1.8*Deg2Rad;
  //   if(anEvent->GetEventID()/3==4) 
  //     u0=3.6*Deg2Rad;
  //   //(((anEvent->GetEventID())/3) -2) *Deg2Rad;
  //   v0=0;
  G4ThreeVector direction(1., u0, v0);
  particleGun->SetParticleMomentumDirection(direction);
  
  
  //particleGun->SetParticleEnergy(parKinE*MeV);
  //double parMom particleGun->GetParticleMomentum();
  double parKinE = particleGun->GetParticleEnergy();
  particleGun-> GeneratePrimaryVertex(anEvent);
  
  anaMan_->SetPrimaryData(y0/mm,z0/mm,-x0/mm+xorg,
			  u0,v0,
			  phi,theta,
			  parMom/GeV,
			  parKinE/GeV,
			  9999);
}
