// -*- C++ -*-

#include "S2SAnaManager.hh"

#include "G4Run.hh"
#include "G4Event.hh"
#include "G4Step.hh"
#include "G4ios.hh"
#include "G4SDManager.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
//#include "G4Poisson.h"

#include "Randomize.hh"

#include <TMath.h>
#include <TRandom3.h>
#include <TSystem.h>

#include "RootHelper.hh"

#include <iomanip>
#include <cmath>
#include <sstream>
#include <vector>

#include "ConfMan.hh"
#include "DCGeomMan.hh"
#include "FuncName.hh"
#include "DCHit.hh"
#include "TOFHit.hh"
#include "ACHit.hh"
#include "WCHit.hh"
#include "VPHit.hh"
#include "HistMan.hh"

namespace
{
using CLHEP::mm;
using CLHEP::MeV;
using CLHEP::ns;
const auto& confMan = ConfMan::GetInstance();
const auto& histMan = HistMan::GetInstance();
const auto qnan = TMath::QuietNaN();
Event event;
std::map<TString, TH1*> hmap;
}

S2SAnaManager::S2SAnaManager()
  : filename_("tmp.root"),
    fActive_(true),
    fTriggered(false),
    DataFile_(),
    m_tree(new TTree("g4s2s", "S-2S simulation"))
{

}

S2SAnaManager::~S2SAnaManager()
{
}

void
S2SAnaManager::BeginOfRun( const G4Run *aRun )
{
  fActive_=true;
  m_file = new TFile(filename_, "recreate");
  static auto obj = new TNamed("conf", confMan.ConfBuf());
  obj->Write();
  static auto git = new TNamed
    ("git", ("\n"+gSystem->GetFromPipe("git log -1")).Data());
  git->Write();
  m_tree->Reset();
  DefineTree();
  for(const auto& sd_name : std::vector<G4String>{
    "TOF" }// S2SDetectorConstruction::GetSDList()
      ){
    G4cout << "   make branch : " << sd_name << G4endl;
    MakeBranch(sd_name);
    MakeHistogram(sd_name);
  }
  for(auto& h: hmap){
    h.second->Reset();
  }
}

void S2SAnaManager::EndOfRun( const G4Run *aRun )
{
  m_file->cd();
  m_tree->Write();
  for(auto& h: hmap){
    h.second->Write();
  }
  m_file->Close();
}

void S2SAnaManager::BeginOfPrimaryAction()
{
  event.x0In = qnan;
  event.y0In = qnan;
  event.u0In = qnan;
  event.v0In = qnan;
  event.uDeg = qnan;
  event.vDeg = qnan;
  event.p0   = qnan;
  event.t0   = qnan;
  event.phi0 = qnan;
  event.theta0 = qnan;
  // event.Id = qnan;
  //  G4cout<<"BeginOfPrimaryAction"<<G4endl;
}

//_____________________________________________________________________________
void
S2SAnaManager::MakeBranch(const G4String& sd_name)
{
  static const Int_t bufsize = 32000;
  m_tree->Branch(sd_name.data(),
                 "std::vector<TParticle>",
                 &event.hits[sd_name], bufsize, -1);
}

//_____________________________________________________________________________
void
S2SAnaManager::MakeHistogram(const G4String& sd_name)
{
  for(const auto& suffix: std::vector<G4String>
        { "Nhits", "HitPat", "X", "Y", "Z", "U", "V",
          "Y%X", "V%U", "U%X", "V%Y" }){
    TString key = sd_name + suffix;
    TString title = sd_name + " " + suffix;
    const auto& params = histMan.Get(key);
    if(G4StrUtil::contains(suffix, "%")){
      hmap[key] = new TH2D(key, title,
                           params.at(0), params.at(1), params.at(2),
                           params.at(3), params.at(4), params.at(5));
    }else{
      hmap[key] = new TH1D(key, title,
                           params.at(0), params.at(1), params.at(2));
    }
  }
}


void S2SAnaManager::SetPrimaryData(double x0, double y0, double z0,
				 double u0, double v0, double phi, double theta,
				 double p0,double t0,int ParIdNb){
  event.x0In = x0; // generated position (x)
  event.y0In = y0; // generated position (y)
  event.z0In = z0; // generated position (z)
  //event.u0In = u0; // x' in rad
  event.u0In = -1.0 * u0; // x' in rad
  event.v0In = v0; // y' in rad
  //event.uDeg = atan(u0)*TMath::RadToDeg(); // x' in deg
  //event.vDeg = atan(v0)*TMath::RadToDeg(); // y' in deg
  event.uDeg = -1.0 * u0 * TMath::RadToDeg(); // x' in deg
  event.vDeg = v0*TMath::RadToDeg(); // y' in deg
  event.phi0 = phi; // Phi in rad
  event.theta0 = theta; // Theta in rad
  event.p0   = p0; // Momentum
  event.t0   = t0; // Kinetic energy
  //  event.Id = ParIdNb;
  //  G4cout<<"setPrimaryData"<<G4endl;
}

void S2SAnaManager::SetProcessData(G4int nP, G4int nN, G4int nL,
				 G4int nSm, G4int nSz, G4int nSp,
				 G4int nXm, G4int nXz, G4int nXsm,
				 G4int nXsz,G4int nPim,G4int nPiz,
				 G4int nPip,G4int nKm,G4int nKp){
  event.nP = nP;
  event.nN = nN;
  event.nL = nL;
  event.nSm = nSm;
  event.nSz = nSz;
  event.nSp = nSp;
  event.nXm = nXm;
  event.nXz = nXz;
  event.nXsm = nXsm;
  event.nXsz = nXsz;
  event.nPim = nPim;
  event.nPiz = nPiz;
  event.nPip = nPip;
  event.nKm = nKm;
  event.nKp = nKp;
}

void S2SAnaManager::BeginOfEvent( const G4Event *anEvent )
{
  //  G4cout<<"BeginOfEvent"<<G4endl;
}

void S2SAnaManager::EndOfEvent( const G4Event *anEvent )
{
  InitializeEvent();

  auto HCE = anEvent->GetHCofThisEvent();
  auto SDMan = G4SDManager::GetSDMpointer();

  //   G4int nhAc=0;
  G4int nhWC=0;
  G4int nhVP=0;

  static const G4int colIdDC = SDMan->GetCollectionID("SDC");
  // static const G4int colIdAC = SDMan->GetCollectionID("AC");
  static const G4int colIdWC = SDMan->GetCollectionID("WC");
  static const G4int colIdVP = SDMan->GetCollectionID("VP");

  auto DCHC = dynamic_cast<DCHitsCollection*>(HCE->GetHC(colIdDC));
  // auto ACHC = dynamic_cast<ACHitsCollection*>(HCE->GetHC(colIdAC));
  auto WCHC = dynamic_cast<WCHitsCollection*>(HCE->GetHC(colIdWC));
  auto VPHC = dynamic_cast<VPHitsCollection*>(HCE->GetHC(colIdVP));

  // if(ACHC) nhAc = ACHC ->entries();
  if(WCHC) nhWC = WCHC ->entries();
  if(VPHC) nhVP = VPHC ->entries();

  G4double pos_res = 0.0*mm; // 0 um
  // G4double pos_res = 0.2*mm; // 200 um
  if(DCHC){
    event.DCNhits = DCHC->entries();
    for(int i=0; i<event.DCNhits; ++i){
      DCHit *aHit = (*DCHC)[i];
      G4int dclayer = aHit->GetLayerID() - 101;
      G4int nh = event.DCNh[dclayer];
      double lx = aHit->GetXLocal()/mm;
      double ly = aHit->GetYLocal()/mm;
      double time = aHit->GetTime()/ns;
      double de = aHit->GetEdep();
      G4ThreeVector gPos = aHit->GetPos();
      G4ThreeVector pVec = aHit->GetMom();
      double mom = pVec.mag();
      event.DCXObs[dclayer] = CLHEP::RandGauss::shoot(lx, pos_res);
      event.DCYObs[dclayer] = CLHEP::RandGauss::shoot(ly, pos_res);
      event.DCX[dclayer] = lx;
      event.DCY[dclayer] = ly;
      event.DCt[dclayer] = time;
      event.DCp[dclayer] = mom;
      event.DCgPosx[dclayer][nh] = gPos.getY();
      event.DCgPosy[dclayer][nh] = gPos.getZ();
      event.DCgPosz[dclayer][nh] = gPos.getX();
      event.DCde[dclayer][nh] = de;
      event.DCNh[dclayer]++;
    }
  }

  for(int i=0;i<6;i++)
    if(event.DCt[i]>0) event.DC1Hit += 1;
  for(int i=6;i<12;i++)
    if(event.DCt[i]>0) event.DC2Hit += 1;
  for(int i=12;i<16;i++)
    if(event.DCt[i]>0) event.DC3Hit += 1;
  for(int i=16;i<22;i++)
    if(event.DCt[i]>0) event.DC4Hit += 1;
  for(int i=22;i<28;i++)
    if(event.DCt[i]>0) event.DC5Hit += 1;

  for(int i=0;i<4;i++){
    for(int j=0;j<6;j++){
      if(i==0 && event.DCt[i*6+j]>-999) event.DC1Hit = 1;
      if(i==1 && event.DCt[i*6+j]>-999) event.DC2Hit = 1;
      if(i==2 && event.DCt[i*6+j]>-999) event.DC3Hit = 1;
      if(i==3 && event.DCt[i*6+j]>-999) event.DC4Hit = 1;
    }
  }

  // ~~~~~~~~~ Q1Flag (T.Gogami, 23Mar2015) ~~~~~~~~~~~~~~~
  G4bool Q1Flag1 = false; // Q1 entrance
  G4bool Q1Flag2 = false; // Q1 exit
  G4bool Q1Flag  = false;
  G4double qx = event.SlitX[0][0]; // at Q1 entrance
  G4double qy = event.SlitY[0]; // at Q1 entrance
  G4double a,b,c;
  a = 8.5;
  b = 10429.0;
  c = 109261.0;
  if( ( qy<a+b/qx+c/qx/qx && qx>56.0 && qy>56.0 )     ||
      ( qy<a-b/qx+c/qx/qx && qx<-56.0 && qy>56.0 )    ||
      ( qy>-a-b/qx-c/qx/qx && qx>56.0 && qy<-56.0 )   ||
      ( qy>-a+b/qx-c/qx/qx && qx<-56.0 && qy<-56.0 )  ||
      ( -56.0<=qx && qx<=56.0 && -293.0<=qy && qy<=293.0 )||
      ( -56.0<=qy && qy<=56.0 && -293.0<=qx && qx<=293.0 )
      ){
    Q1Flag1=true;
  }
  else Q1Flag1=false;

  qx = event.SlitX[1][0]; // at Q1 exit
  qy = event.SlitY[1]; // at Q1 exit
  if( ( qy<a+b/qx+c/qx/qx && qx>56.0 && qy>56.0 )     ||
      ( qy<a-b/qx+c/qx/qx && qx<-56.0 && qy>56.0 )    ||
      ( qy>-a-b/qx-c/qx/qx && qx>56.0 && qy<-56.0 )   ||
      ( qy>-a+b/qx-c/qx/qx && qx<-56.0 && qy<-56.0 )  ||
      ( -56.0<=qx && qx<=56.0 && -293.0<=qy && qy<=293.0 )||
      ( -56.0<=qy && qy<=56.0 && -293.0<=qx && qx<=293.0 )
      ){
    Q1Flag2=true;
  }
  else Q1Flag2=false;

  if(Q1Flag1 && Q1Flag2) Q1Flag=true;
  else Q1Flag=false;

  // ~~~~~~~~~ Q2Flag (T.Gogami, 23Mar2015) ~~~~~~~~~~~~~~~
  G4bool Q2Flag1 = false;
  G4bool Q2Flag2 = false;
  G4bool Q2Flag  = false;
  qx = event.SlitX[2][0]; // at Q2 entrance
  qy = event.SlitY[2]; // at Q2 entrance
  a = 0.7;
  b = 16073.4;
  c = 5202.96;
  if( ( qy<a+b/qx+c/qx/qx && qx>60.0 && qy>32.4 )     ||
      ( qy<a-b/qx+c/qx/qx && qx<-60.0 && qy>32.4 )    ||
      ( qy>-a-b/qx-c/qx/qx && qx>60.0 && qy<-32.4 )   ||
      ( qy>-a+b/qx-c/qx/qx && qx<-60.0 && qy<-32.4 )  ||
      ( -60.0<=qx && qx<=60.0 && -270.0<=qy && qy<=270.0 )||
      ( -32.4<=qy && qy<=32.4 && -503.0<=qx && qx<=503.0 )
      ){
    Q2Flag1=true;
  }
  else Q2Flag1=false;

  qx = event.SlitX[3][0]; // at Q2 exit
  qy = event.SlitY[3]; // at Q2 exit
  if( ( qy<a+b/qx+c/qx/qx && qx>60.0 && qy>32.4 )     ||
      ( qy<a-b/qx+c/qx/qx && qx<-60.0 && qy>32.4 )    ||
      ( qy>-a-b/qx-c/qx/qx && qx>60.0 && qy<-32.4 )   ||
      ( qy>-a+b/qx-c/qx/qx && qx<-60.0 && qy<-32.4 )  ||
      ( -60.0<=qx && qx<=60.0 && -270.0<=qy && qy<=270.0 )||
      ( -32.4<=qy && qy<=32.4 && -503.0<=qx && qx<=503.0 )
      ){
    Q2Flag2=true;
  }
  else Q2Flag2=false;

  if(Q2Flag1 && Q2Flag2) Q2Flag=true;
  else Q2Flag=false;

  event.Q1Trig = Q1Flag;
  event.Q2Trig = Q2Flag;

  //G4cout<<"nhitSl1 = "<<nhitSl1<<G4endl;
  /*
    if(event.Slitp[0]<100) {
    G4cout<<"this event is wrong !!! Slitt[0]="<<event.Slitt[0]<<" Slitp[0]"<<event.Slitp[0]<<G4endl;
    }
  */

  {
    static const auto id = SDMan->GetCollectionID("TOF");
    auto HC = dynamic_cast<TOFHitsCollection*>(HCE->GetHC(id));
    for(G4int i=0, n=HC->entries(); i<n; ++i){
      SetHitData((*HC)[i]);
    }
    SetNhits("TOF", HC->entries());
  }

  // event.TOFNhits = (int)nhTof;
  // G4double t_res = 0.090*ns; // Sigma = 90 ps (From cosmic-ray test)
  // G4bool TOFTrig = false;
  // if( TOFHC ){
  //   for( int i=0; i<nhTof; ++i ){
  //     TOFHit *aHit = (*TOFHC)[i];
  //     G4int TOFHitSeg = aHit->GetLayerID(); // Hit Segment
  //     double TOFtime = aHit->GetTime()/ns;  // Hit timing
  //     double TOFdE   = aHit->GetEdep()/MeV; // Energy deposite
  //     event.toftime[TOFHitSeg] = event.toftime[TOFHitSeg] + TOFtime;
  //     event.toftime_reso[TOFHitSeg] = event.toftime_reso[TOFHitSeg] + CLHEP::RandGauss::shoot(TOFtime, t_res);
  //     event.tofdE[TOFHitSeg]   = event.tofdE[TOFHitSeg] + TOFdE;
  //     event.tofn[TOFHitSeg]++;
  //     //G4cout<<"TofSeg"<<TOFHitSeg<<" event.TOFC["<<TOFHitSeg<<"]"<<event.TOFt[TOFHitSeg]<<G4endl;
  //   }

  //   for(int i=0 ; i<18 ; i++){
  //     if(event.toftime[i]>0.0){
  //       // ----- Mean value of the timing -----
  //       event.toftime[i] = event.toftime[i] / event.tofn[i];  // w/o resolution
  //       event.toftime_reso[i] = event.toftime_reso[i] / event.tofn[i]; // w/ resolution

  //       //if(event.toftime[i]>10.0 && event.tofdE[i]>0.0){
  //       if( event.toftime[i]>10.0 ){
  //         TOFTrig = true;
  //       }
  //       else TOFTrig = false;

  //     }
  //     else{
  //       event.toftime[i] = qnan;
  //       event.tofdE[i]   = qnan;
  //     }
  //   }
  // }

  // event.TOFTrig = TOFTrig;


  // event.VDTrig  = SlitFlag;

  //for(int j=0;j<NumTOFSeg;j++){
  //  if(event.TOFt[j]>0) event.TOFHit = 1;
  //}

  //   //AC
  //   event.ACNhits = nhAc;
  //   if( ACHC ){
  //     for( int i=0;i<nhAc;++i){
  //       ACHit *aHit = (*ACHC)[i];
  //       double actime = aHit->GetTime()/ns;
  //       double lx = aHit->GetXLocal()/mm;
  //       double ly = aHit->GetYLocal()/mm;
  //       G4ThreeVector pVec = aHit->GetMom();
  //       G4int mom = pVec.mag();
  //       event.ACX = lx;
  //       event.ACY = ly;
  //       event.ACt = actime;
  //       event.ACp = mom;
  //     }
  //   }
  //   if(event.ACt>0) event.ACHit = 1;

  // // ~~~~~~~~~~~~~~ Water Cherenkov detecctor  ~~~~~~~~~~~~~~~~~~~~~~~~~
  // event.WCNhits = (int)nhWC;
  // //G4double t_res = 0.090*ns; // Sigma = 90 ps (From cosmic-ray test)
  // G4bool WCTrig = false;
  // //TRandom3* randtemp = new TRandom3();
  // //double aaatemp=0.0;
  // //aaatemp = randtemp->Uniform(1.0,65539.0);
  // //rrr = rrr%65539;
  // //TRandom3* nperand = new TRandom3(aaatemp);
  // double wcnpe1, wcnpe2;
  // if( WCHC ){
  //   for( int i=0; i<nhWC; i++ ){
  //     WCHit *aHit = (*WCHC)[i];
  //     G4int WCHitSeg = aHit->GetLayerID(); // Hit Segment
  //     double WCtime = aHit->GetTime()/ns;  // Hit timing
  //     double WCdE   = aHit->GetEdep()/MeV; // Energy deposite
  //     double WCNPE  = aHit->GetNPE();  // N.P.E.
  //     event.wctime[WCHitSeg] = event.wctime[WCHitSeg] + WCtime;
  //     //event.wctime_reso[WCHitSeg] = event.wctime_reso[WCHitSeg] + CLHEP::RandGauss::shoot(WCtime, t_res);
  //     event.wcdE[WCHitSeg]   = event.wcdE[WCHitSeg] + WCdE;
  //     event.wcnpe[WCHitSeg]  = event.wcnpe[WCHitSeg] + WCNPE;
  //     //event.wcnpe[WCHitSeg]  = event.wcnpe[WCHitSeg] + nperand->PoissonD(WCNPE);
  //     event.wcn[WCHitSeg]++;
  //     //G4cout<<"TofSeg"<<TOFHitSeg<<" event.TOFC["<<TOFHitSeg<<"]"<<event.TOFt[TOFHitSeg]<<G4endl;
  //   }
  //   for(int i=0 ; i<12 ; i++){
  //     if(event.wctime[i]>0.0){
  //       // ----- Mean value of the timing -----
  //       event.wctime[i] = event.wctime[i] / event.wcn[i];  // w/o resolution
  //       //event.wctime_reso[i] = event.wctime_reso[i] / event.wcn[i]; // w/ resolution
  //       ///event.wcnpe[i] = G4Poisson(event.wcnpe[i]) * 2.0; // up and down PMT --> * 2
  //       /*
  //         event.wcnpe[i] = nperand->PoissonD(event.wcnpe[i]) * 2.0;
  //       */
  //       //event.wcnpe[i] = nperand->PoissonD(event.wcnpe[i]);
  //       wcnpe1 = nperand->PoissonD(event.wcnpe[i]);
  //       wcnpe2 = nperand->PoissonD(event.wcnpe[i]);
  //       event.wcnpe[i] = wcnpe1+wcnpe2;
  //       //event.wcnpe[i] = nperand->Gaus(event.wcnpe[i],10.0); // 10:width from cosmic-ray data
  //       //event.wcnpe[i] = event.wcnpe[i] * 2.0;
  //       //if(event.wctime[i]>10.0 && event.wcdE[i]>0.0){
  //       if(event.wctime[i]>10.0 &&
  //          event.wcnpe[i]>0.0){// &&
  //         //event.wcnpe[i]>30.0){
  //         //event.wcnpe[i] = nperand->PoissonD(event.wcnpe[i]);
  //         WCTrig = true;
  //         if(i<6){
  //           event.wctime1[i] = event.wctime[i];
  //           event.wcnpe1[i]  = event.wcnpe[i];
  //           event.wcdE1[i]   = event.wcdE[i];
  //           event.wcnpe1[i] =  CLHEP::RandGauss::shoot(event.wcnpe1[i], event.wcnpe1[i]*0.06); // +/-6% (Gogami)
  //           //event.wcnpe1[i] = event.wcnpe1[i]
  //           //+ (event.wcnpe1[i] * (G4UniformRand()-0.5)*2.0 * 0.14);  // +/- 14% (Takenaka value)
  //           //+ (event.wcnpe1[i] * (G4UniformRand()-0.5)*2.0 * 0.2);  // +/- 20%
  //           //+ (event.wcnpe1[i] * (G4UniformRand()-0.5)*2.0 * 1.0);  // +/- 100%
  //           //if (i==5)G4cout << event.wcnpe1[i] << G4endl; // for check
  //         }
  //         else{
  //           event.wctime2[i-6] = event.wctime[i];
  //           event.wcnpe2[i-6]  = event.wcnpe[i];
  //           event.wcdE2[i-6]   = event.wcdE[i];
  //           event.wcnpe2[i-6] =  CLHEP::RandGauss::shoot(event.wcnpe2[i-6], event.wcnpe2[i-6]*0.06); // +/-6% (Gogami)
  //           //event.wcnpe2[i-6] = event.wcnpe2[i-6]
  //           //+ (event.wcnpe2[i-6] * (G4UniformRand()-0.5)*2.0 * 0.14); // +/- 14% (Takenaka value)
  //           //+ (event.wcnpe2[i-6] * (G4UniformRand()-0.5)*2.0 * 0.2); // +/- 20%
  //           //+ (event.wcnpe2[i-6] * (G4UniformRand()-0.5)*2.0 * 1.0); // +/- 100%
  //         }
  //       }
  //       else WCTrig = false;
  //     }
  //     else{
  //       event.wctime[i] = qnan;
  //       event.wcdE[i]   = qnan;
  //       event.wcnpe[i]  = qnan;
  //     }
  //   }
  // }
  // event.WCTrig = WCTrig;

  //   event.WCNhits = nhWc;
  //   if( WCHC ){
  //     for( int i=0; i<nhWc; ++i ){
  //       WCHit *aHit = (*WCHC)[i];
  //       G4int WCHitSeg = aHit->GetLayerID();
  //       double wctime = aHit->GetTime()/ns;
  //       G4ThreeVector pVec = aHit->GetMom();
  //       G4int mom = pVec.mag();
  //       event.WCt[WCHitSeg] = wctime;
  //       event.WCp[WCHitSeg] = mom;
  //       //G4cout<<"TofSeg"<<TOFHitSeg<<" event.TOFC["<<TOFHitSeg<<"]"<<event.TOFt[TOFHitSeg]<<G4endl;
  //     }
  //   }
  //   for(int j=0;j<NumWCSeg;j++){
  //     if(event.WCt[j]>0) event.WCHit = 1;
  //   }

  //   if( DataFile_.is_open() ){
  //     PrintHitsInformation( anEvent, DataFile_ );
  //   }

  //  //if(event.DC3Hit>0){
  //  if(1){
  //    //    m_tree->Fill();
  //    //if( DataFile_.is_open() ){
  //    if(0){
  //      PrintHitsInformation( anEvent, DataFile_ );
  //    }
  //  }

  //  G4cout<<"EndOf S2SAnaManager EndOfEvent"<<G4endl;
  //if(m_tree->GetEntries()>1&&m_tree->GetEntries()%100==0)
  //if(1)
  // m_file->Write();
  //if(TOFTrig==true){

  // if(SlitFlag==true){
  //   //m_tree->Fill();
  // }
  m_tree->Fill();
}

//_____________________________________________________________________________
void
S2SAnaManager::SetNhits(const G4String& sd_name, G4int nhits)
{
  hmap[sd_name + "Nhits"]->Fill(nhits);
}

//_____________________________________________________________________________
void
S2SAnaManager::SetHitData(const VHitInfo* hit)
{
  if(hit && hit->GetParticle()){
    const auto& name = hit->GetDetectorName();
    const auto& p = hit->GetParticle();
    event.hits.at(name).push_back(*p);
    hmap[name + "HitPat"]->Fill(p->GetMother(1));
    hmap[name + "X"]->Fill(p->Vx());
    hmap[name + "Y"]->Fill(p->Vy());
    hmap[name + "Z"]->Fill(p->Vz());
    hmap[name + "U"]->Fill(p->Px()/p->Pz());
    hmap[name + "V"]->Fill(p->Py()/p->Pz());
    hmap[name + "Y%X"]->Fill(p->Vx(), p->Vy());
    hmap[name + "V%U"]->Fill(p->Px()/p->Pz(), p->Py()/p->Pz());
    hmap[name + "U%X"]->Fill(p->Vx(), p->Px()/p->Pz());
    hmap[name + "V%Y"]->Fill(p->Vy(), p->Py()/p->Pz());
  }
}


void S2SAnaManager::SaveFile() const
{
  if( fActive_ ){
    gFile->Write();
  }
}

void S2SAnaManager::Terminate() const
{
  if( fActive_ ){
    m_file->Write();
    m_file->Close();
    G4cout << "[S2SAnaManager] Terminate() done" << G4endl;
  }
}

void S2SAnaManager::ShowStatus() const
{
  G4cout << "Analyzer Status\n"
	 << "  File : " << filename_ << "\n"
	 << "  S2SAnaManager : " << fActive_ << "\n"
	 << G4endl;
}

void S2SAnaManager::InitializeEvent()
{
  for(int i=0;i<NumDC;i++){
    event.DCXObs[i] = qnan;
    event.DCYObs[i] = qnan;
    event.DCX[i] = qnan;
    event.DCY[i] = qnan;
    event.DCt[i] = qnan;
    event.DCp[i] = qnan;
    event.DCNh[i] = 0;
    for( int ihit = 0; ihit < MaxHits; ihit++ ){
      event.DCgPosx[i][ihit] = qnan;
      event.DCgPosy[i][ihit] = qnan;
      event.DCgPosz[i][ihit] = qnan;
      event.DCde[i][ihit] = qnan;
    }
  }
  //   for(int i=0;i<NumTOFSeg;i++){
  //     event.TOFt[i] = qnan;
  //     event.TOFtObs[i] = qnan;
  //   }
  // event.TOFAll = qnan;
  //  for(int i=0;i<8;i++){
  //    event.SlitX[i]    = qnan;
  //    event.SlitY[i]    = qnan;
  //    event.SlituDeg[i] = qnan;
  //    event.SlitvDeg[i] = qnan;
  //    event.Slitt[i]    = qnan;
  //    event.SlitMom[i]  = qnan;
  //    event.Slitp[i]    = qnan;
  //    event.SlitNh[i]   = 0;
  //    event.SlitNP[i]   = 0;
  //    event.SlitNK[i]   = 0;
  //    event.SlitNPi[i]  = 0;
  //    event.SlitF[i]    = 0;
  //  }
  for(int i=0 ; i<18 ; i++){
    // ------- Virtual detectors (slits) ------
    if(i<11){
      for( int ihit = 0; ihit < MaxHits; ihit++ ) event.SlitX[i][ihit]    = qnan;
      event.SlitY[i]    = qnan;
      event.SlituDeg[i] = qnan;
      event.SlitvDeg[i] = qnan;
      event.Slitt[i]    = qnan;
      event.SlitMom[i]  = qnan;
      event.Slitp[i]    = qnan;
      event.SlitNh[i]   = 0;
      event.SlitNP[i]   = 0;
      event.SlitNK[i]   = 0;
      event.SlitNPi[i]  = 0;
      event.SlitF[i]    = 0;
    }

    // ------- TOF detector ------
    //event.tofco[i]        = qnan;
    event.toftime[i]      = 0.0;
    event.toftime_reso[i] = 0.0;
    event.tofdE[i]        = 0.0;
    event.tofn[i]         = 0.0;

    // ------- Water Cherenkov detector ------
    if(i<12){
      event.wctime[i] = 0.0;
      event.wcdE[i]   = 0.0;
      event.wcn[i]    = 0.0;
      event.wcnpe[i]  = 0.0;
    }
    if(i<6){
      event.wctime1[i] = 0.0;
      event.wcdE1[i]   = 0.0;
      event.wcn1[i]    = 0.0;
      event.wcnpe1[i]  = 0.0;
      event.wctime2[i] = 0.0;
      event.wcdE2[i]   = 0.0;
      event.wcn2[i]    = 0.0;
      event.wcnpe2[i]  = 0.0;
    }
  }

  event.DCNhits = 0;
  event.DC1Hit = 0;
  event.DC2Hit = 0;
  event.DC3Hit = 0;
  event.DC4Hit = 0;
  event.DC5Hit = 0;
  // event.TOFNhits = -999;
  // event.ACNhits = -999;
  // event.WCNhits = -999;
  // event.TOFHit = -999;
  // event.ACHit = -999;
  // event.WCHit = -999;
  // event.ACX = qnan;
  // event.ACY = qnan;
  // event.ACt = qnan;
  // for(int i=0;i<NumWCSeg;i++){
  //   event.WCt[i] = qnan;
  //   event.WCp[i] = qnan;
  // }
  //  G4cout<<"S2SAnaManager InitializeEvent"<<G4endl;
}

//void S2SAnaManager::DefineHistograms()
void S2SAnaManager::DefineTree()
{
  m_tree->Branch("x0",&event.x0In, "x0/D");
  m_tree->Branch("y0",&event.y0In, "y0/D");
  m_tree->Branch("z0",&event.z0In, "z0/D");
  //   m_tree->Branch("u0In",&event.u0In, "u0In/D");
  //   m_tree->Branch("v0In",&event.v0In, "v0In/D");
  m_tree->Branch("xp0",&event.uDeg, "xp0/D");
  m_tree->Branch("yp0",&event.vDeg, "yp0/D");
  m_tree->Branch("phi0",&event.phi0, "phi0/D");
  m_tree->Branch("theta0",&event.theta0, "theta0/D");
  m_tree->Branch("p0",&event.p0,   "p0/D");

  //  m_tree->Branch("t0",&event.t0,   "t0/D");
  // m_tree->Branch("Id",&event.Id, "Id/I");
  //   m_tree->Branch("nP",&event.nP,"nP/I");
  //   m_tree->Branch("nN",&event.nN,"nN/I");
  //   m_tree->Branch("nL",&event.nL,"nL/I");
  //   m_tree->Branch("nSm",&event.nSm,"nSm/I");
  //   m_tree->Branch("nSz",&event.nSz,"nSz/I");
  //   m_tree->Branch("nSp",&event.nSp,"nSp/I");
  //   m_tree->Branch("nXm",&event.nXm,"nXm/I");
  //   m_tree->Branch("nXz",&event.nXz,"nXz/I");
  //   m_tree->Branch("nXsm",&event.nXsm,"nXsm/I");
  //   m_tree->Branch("nXsz",&event.nXsz,"nXsz/I");
  //   m_tree->Branch("nPim",&event.nPim,"nPim/I");
  //   m_tree->Branch("nPiz",&event.nPiz,"nPiz/I");
  //   m_tree->Branch("nPip",&event.nPip,"nPip/I");
  //   m_tree->Branch("nKm",&event.nKm,"nKm/I");
  //   m_tree->Branch("nKp",&event.nKp,"nKp/I");

  m_tree->Branch("DCXObs",event.DCXObs, "DCXObs[26]/D");
  m_tree->Branch("DCYObs",event.DCYObs, "DCYObs[26]/D");
  m_tree->Branch("DCX",event.DCX, "DCX[26]/D");
  m_tree->Branch("DCY",event.DCY, "DCY[26]/D");
  m_tree->Branch("DCt",event.DCt, "DCt[26]/D");
  m_tree->Branch("DCp",event.DCp, "DCp[26]/D");
  m_tree->Branch("DCNh",event.DCNh, "DCNh[26]/I");
  m_tree->Branch("DCgPosx",event.DCgPosx, Form("DCgPosx[26][%d]/D", MaxHits));
  m_tree->Branch("DCgPosy",event.DCgPosy, Form("DCgPosy[26][%d]/D", MaxHits));
  m_tree->Branch("DCgPosz",event.DCgPosz, Form("DCgPosz[26][%d]/D", MaxHits));
  m_tree->Branch("DCde",   event.DCde, Form("DCde[26][%d]/D", MaxHits));
  m_tree->Branch("DCNhits",&event.DCNhits, "DCNhits/I");
  m_tree->Branch("DC1Hit", &event.DC1Hit, "DC1Hit/I");
  m_tree->Branch("DC2Hit", &event.DC2Hit, "DC2Hit/I");
  m_tree->Branch("DC3Hit", &event.DC3Hit, "DC3Hit/I");
  m_tree->Branch("DC4Hit", &event.DC4Hit, "DC4Hit/I");
  m_tree->Branch("DC5Hit", &event.DC5Hit, "DC5Hit/I");
  m_tree->Branch("vdxp", event.SlituDeg, "vdxp[11]/D");
  m_tree->Branch("vdyp", event.SlitvDeg, "vdyp[11]/D");
  /*
    m_tree->Branch("SlitX", event.SlitX, "SlitX[8]/D");
    m_tree->Branch("SlitY", event.SlitY, "SlitY[8]/D");
    m_tree->Branch("Slitt", event.Slitt, "Slitt[8]/D");
    m_tree->Branch("SlitMom", event.SlitMom, "SlitMom[8]/D");
    m_tree->Branch("Slitp", event.Slitp, "Slitp[8]/D");
    m_tree->Branch("SlitNh", event.SlitNh, "SlitNh[8]/I");
    m_tree->Branch("SlitNP", event.SlitNP, "SlitNP[8]/I");
    m_tree->Branch("SlitNK", event.SlitNK, "SlitNK[8]/I");
    m_tree->Branch("SlitNPi", event.SlitNPi, "SlitNPi[8]/I");
    m_tree->Branch("SlitF", event.SlitF, "SlitF[8]/I");
  */

  // Order -->
  // vd8 | TOF | vd9 | AC | vd10 | WC | vd11 (vd[10]) | SDC1
  m_tree->Branch("vdx",   event.SlitX, Form("vdx[11][%d]/D", MaxHits));
  m_tree->Branch("vdy",   event.SlitY, "vdy[11]/D");
  m_tree->Branch("vdtime",event.Slitt, "vdtime[11]/D");
  m_tree->Branch("vdmom", event.SlitMom, "vdmom[11]/D");
  m_tree->Branch("vdpath",event.Slitp, "vdpath[11]/D");
  m_tree->Branch("vdNh",  event.SlitNh, "vdNh[11]/I");
  m_tree->Branch("vdNP",  event.SlitNP, "vdNP[11]/I");
  m_tree->Branch("vdNK",  event.SlitNK, "vdNK[11]/I");
  m_tree->Branch("vdNPi", event.SlitNPi, "vdNPi[11]/I");
  m_tree->Branch("vdPID", event.SlitF, "vdPID[11]/I");
  // m_tree->Branch("TOFAll", &event.TOFAll, "TOFAll/D");
  // m_tree->Branch("TOFt",event.TOFt,   "TOFt[17]/D");
  // m_tree->Branch("TOFtObs",event.TOFtObs,   "TOFtObs[17]/D");
  //m_tree->Branch("tofco",  &event.tofco,   "tofco[17]/I");

  m_tree->Branch("tofnhits",&event.TOFNhits, "tofnhits/I");
  m_tree->Branch("toftime", &event.toftime,  "toftime[18]/D");
  m_tree->Branch("tofdE",   &event.tofdE,    "tofdE[18]/D");

  m_tree->Branch("wcnhits",&event.WCNhits,  "wcnhits/I");
  //m_tree->Branch("wctime", &event.wctime,   "wctime[12]/D");
  //m_tree->Branch("wcdE",   &event.wcdE,     "wcdE[12]/D");
  //m_tree->Branch("wcnpe",  &event.wcnpe,    "wcnpe[12]/D");
  m_tree->Branch("wctime1", &event.wctime1,   "wctime1[6]/D");
  m_tree->Branch("wcdE1",   &event.wcdE1,     "wcdE1[6]/D");
  m_tree->Branch("wcnpe1",  &event.wcnpe1,    "wcnpe1[6]/D");
  m_tree->Branch("wctime2", &event.wctime2,   "wctime2[6]/D");
  m_tree->Branch("wcdE2",   &event.wcdE2,     "wcdE2[6]/D");
  m_tree->Branch("wcnpe2",  &event.wcnpe2,    "wcnpe2[6]/D");

  m_tree->Branch("Q1Trig", &event.Q1Trig,  "Q1Trig/B");
  m_tree->Branch("Q2Trig", &event.Q2Trig,  "Q2Trig/B");
  m_tree->Branch("VDTrig", &event.VDTrig,  "VDTrig/B");
  // m_tree->Branch("TOFTrig",&event.TOFTrig, "TOFTrig/B");

  m_tree->Branch("WCTrig", &event.WCTrig,  "WCTrig/B");
  // m_tree->Branch("TOFHit", &event.TOFHit, "TOFHit/I");
  // m_tree->Branch("ACX", &event.ACX, "ACX/D");
  // m_tree->Branch("ACY", &event.ACY, "ACY/D");
  // m_tree->Branch("ACt", &event.ACt, "ACt/D");
  // m_tree->Branch("ACp", &event.ACp, "ACp/D");
  // m_tree->Branch("ACNhits", &event.ACNhits, "ACNhits/I");
  // m_tree->Branch("ACHit", &event.ACHit, "ACHit/I");
  // m_tree->Branch("WCt", event.WCt,  "WCt[16]/D");
  // m_tree->Branch("WCp", event.WCp,  "WCp[16]/D");
  // m_tree->Branch("WCNhits", &event.WCNhits, "WCNhits/I");
  // m_tree->Branch("WCHit", &event.WCHit, "WCHit/I");
  //G4cout<<"tree and Branch is defined"<<G4endl;

}

void S2SAnaManager::
PrintHitsInformation( const G4Event *anEvent,
		      std::ostream &ost ) const
{
  //   G4cout<<"PrintHits is called"<<G4endl;

  //  int GeomFlag = confMan->GeomFlag();

  G4HCofThisEvent *HCE = anEvent->GetHCofThisEvent();
  G4SDManager *SDMan = G4SDManager::GetSDMpointer();

  G4int nhDc=0;
  DCHitsCollection    *DCHC;
  G4int colIdDC = SDMan->GetCollectionID( "BcSD"/*"DCCollection"*/ );
  DCHC =  dynamic_cast<DCHitsCollection *>(  HCE->GetHC( colIdDC ) );
  if( DCHC )     nhDc     = DCHC ->entries();

  std::ios::fmtflags oldFlags = ost.flags();
  std::size_t preSiz = ost.precision();
  ost.setf( std::ios::fixed );

  ost.precision(5);
  ost << event.x0In << std::setw(10)
      << event.y0In << std::setw(10)
      <<" "<< -4320.62-event.z0In << std::setw(10);
  ost.precision(5);

  ost << event.p0     << std::setw(10)
      << event.theta0 << std::setw(10)
      << event.phi0   << std::setw(10);
  ost.precision(7);

  ost << 1.8 <<std::setw(10);
  ost << 0   <<std::setw(10);
  ost << 0   <<std::setw(10);

  /*
    ost << event.DCXObs[2] << std::setw(15)
    << event.DCYObs[2] << std::setw(15)
    << event.DCXObs[8] << std::setw(15)
    << event.DCYObs[8] << std::setw(15)
    << event.DCXObs[14] << std::setw(15)
    << event.DCYObs[14] << std::setw(15)
    << event.DCXObs[20] << std::setw(15)
    << event.DCYObs[20] << std::setw(15);
  */
  ost.precision(5);

  //G4cout<<"Printed"<<G4endl;
  //G4cout<<"osf="<<ost<<G4endl;
  //G4cout<<"x0="<<event.x0In<<" y0="<<event.y0In<<G4endl;

  const DCGeomMan & geomMan=DCGeomMan::GetInstance();

  //DC information
  for( int i=0; i<nhDc; ++i ){
    DCHit *aHit = (*DCHC)[i];
    int layer = aHit->GetLayerID();
    double lx = aHit->GetXLocal()/mm;
    double ly = aHit->GetYLocal()/mm;

    double tilt  = geomMan.GetTiltAngle(layer-100);
    double reso  = geomMan.GetResolution(layer-100);

    double lp = lx*cos(tilt*TMath::DegToRad()) + ly*sin(tilt*TMath::DegToRad());
    lp += CLHEP::RandGauss::shoot(0., reso);

    ost << std::setw(12) << layer-100
  	<< std::setw(12) << lp;
  }

  // G4double nhTOF=0;
  // TOFHitsCollection    *TOFHC;
  // G4int colIdTOF = SDMan->GetCollectionID( "TOFSD"/*"DCCollection"*/ );
  // TOFHC =  dynamic_cast<TOFHitsCollection *>(  HCE->GetHC( colIdTOF ) );
  // if( TOFHC )     nhTOF     = TOFHC ->entries();

  // for( int i=0; i<nhTOF; ++i ){
  //   TOFHit *aHit = (*TOFHC)[i];
  //   double time = aHit->GetTime();
  //   double de = aHit->GetEdep();
  //   if(i==nhTOF-1){
  //     ost << std::setw(12) << 61
  //         << std::setw(12) << time;
  //     de=1.0;
  //     ost << std::setw(12) << 62 <<" "
  //         << std::setw(12) << de;
  //   }
  // }

  ost << std::setw(5) << -1 << std::endl;

  ost.flags( oldFlags );
  ost.precision( preSiz );

  ost << std::endl;

}

void S2SAnaManager::SetDataFile( const char *datafile )
{
  DataFile_.open( datafile );
}
