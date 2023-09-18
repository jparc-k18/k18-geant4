/*
  "S2SAnalysis.cc"
  
  2013/4  S. Kanatsuki
  Modified by Toshi Gogami , 10Nov2014
*/

#include "S2SAnalysis.hh"
#include "DCHit.hh"
#include "TOFHit.hh"
#include "ACHit.hh"
#include "WCHit.hh"
#include "SlHit.hh"

#include "RadDeg.hh"

#include "G4Run.hh"
#include "G4Event.hh"
#include "G4Step.hh"
#include "G4ios.hh"
#include "G4SDManager.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
//#include "G4Poisson.h"
#include "TRandom3.h"

#include "Randomize.hh"
#include "RootHelper.hh"

#include "ConfMan.hh"
#include "DCGeomMan.hh"

#include <iomanip>
#include <cmath>
#include <sstream>
#include <vector>

namespace
{
  using namespace root;
  using namespace CLHEP;
}

S2SAnalysis::S2SAnalysis( G4String file )
  : filename_(file), fActive_(true), fTriggered(false), DataFile_()
{
  //G4cout<<"S2SAnalysis constructor"<<G4endl;
  //DefineHistograms();
  //DefineTree();
  nperand = new TRandom3();
}

S2SAnalysis::S2SAnalysis()
  : filename_("../analysis/root/DRAGON.root"),
    fActive_(true), fTriggered(false), DataFile_()
{
  //G4cout<<"S2SAnalysis constructor"<<G4endl;
  //DefineHistograms();
  //nperand = new TRandom3();
}

S2SAnalysis::S2SAnalysis(ConfMan* conf_)
  : filename_("../analysis/root/DRAGON.root"),
    fActive_(true), fTriggered(false), DataFile_(),
    conf(conf_)
{
  //G4cout<<"S2SAnalysis constructor"<<G4endl;
  //DefineHistograms();
  filename_ = conf->GetOutROOTFile();
  //TRandom3* nperand = new TRandom3(aaatemp);
  nperand = new TRandom3();
}


S2SAnalysis::~S2SAnalysis()
{
  //  G4cout<<"!!! S2SAnalysis destructor !!!"<<G4endl;
  //   SaveFile();
}

void S2SAnalysis::BeginOfRun( const G4Run *aRun )
{
  //G4cout<<"S2SAnalysis BeginOfRun"<<G4endl;
  //for(int i=0 ; i<10000 ; i++){
  //  G4cout << " taaatattttttttAAAAtttttttttttt " << G4endl;
  //}
  DefineTree();
}

void S2SAnalysis::EndOfRun( const G4Run *aRun )
{
  //  G4cout<<"S2SAnalysis EndOfRun"<<G4endl;
 TTree *tree = dynamic_cast<TTree *>(anafile->Get("s2s_g"));
  tree->Write();
  anafile->Close();
}

void S2SAnalysis::BeginOfPrimaryAction()
{
  event.x0In = -999.;
  event.y0In = -999.;
  event.u0In = -999.;
  event.v0In = -999.;
  event.uDeg = -999.;
  event.vDeg = -999.;
  event.p0   = -999.;
  event.t0   = -999.;
  event.phi0 = -999.;
  event.theta0 = -999.;
  // event.Id = -999.;
  //  G4cout<<"BeginOfPrimaryAction"<<G4endl;
}

void S2SAnalysis::SetPrimaryData(double x0, double y0, double z0,
				 double u0, double v0, double phi, double theta,
				 double p0,double t0,int ParIdNb){
  event.x0In = x0; // generated position (x)
  event.y0In = y0; // generated position (y)
  event.z0In = z0; // generated position (z)
  //event.u0In = u0; // x' in rad
  event.u0In = -1.0 * u0; // x' in rad
  event.v0In = v0; // y' in rad
  //event.uDeg = atan(u0)*Rad2Deg; // x' in deg
  //event.vDeg = atan(v0)*Rad2Deg; // y' in deg
  event.uDeg = -1.0 * u0 * Rad2Deg; // x' in deg
  event.vDeg = v0*Rad2Deg; // y' in deg
  event.phi0 = phi; // Phi in rad
  event.theta0 = theta; // Theta in rad
  event.p0   = p0; // Momentum
  event.t0   = t0; // Kinetic energy 
  //  event.Id = ParIdNb;
  //  G4cout<<"setPrimaryData"<<G4endl;
}

void S2SAnalysis::SetProcessData(G4int nP, G4int nN, G4int nL,
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

void S2SAnalysis::BeginOfEvent( const G4Event *anEvent )
{
  //  G4cout<<"BeginOfEvent"<<G4endl;
}

void S2SAnalysis::EndOfEvent( const G4Event *anEvent )
{
  // =========== Tree =======================================
  TTree *tree = dynamic_cast<TTree *>(anafile->Get("s2s_g"));
  InitializeEvent();

  // ========== Hit Collection =========================
  G4HCofThisEvent *HCE = anEvent->GetHCofThisEvent();
  G4SDManager *SDMan = G4SDManager::GetSDMpointer();

  G4int nhDc=0;
  G4int nhTof=0;
  //   G4int nhAc=0;
  G4int nhWC=0;
  G4int nhSl=0;
  
  G4int colIdDC = SDMan->GetCollectionID( "BcSD"/*"DCCollection"*/ );
  DCHitsCollection *DCHC =  dynamic_cast<DCHitsCollection *>(  HCE->GetHC( colIdDC )  );
  
  G4int colIdTOF = SDMan->GetCollectionID("TOFSD");
  TOFHitsCollection *TOFHC = dynamic_cast<TOFHitsCollection *>( HCE->GetHC( colIdTOF ) );
  //   G4int colIdAC = SDMan->GetCollectionID("AC");
  //   ACHitsCollection *ACHC = dynamic_cast<ACHitsCollection *>( HCE->GetHC( colIdAC ) );
  G4int colIdWC = SDMan->GetCollectionID("WCSD");
  WCHitsCollection *WCHC = dynamic_cast<WCHitsCollection *>( HCE->GetHC( colIdWC ) );
  
  G4int colIdSl = SDMan->GetCollectionID("SlSD");
  SlHitsCollection *SlHC = dynamic_cast<SlHitsCollection *>( HCE->GetHC( colIdSl ) );
  
   if( DCHC )     nhDc     = DCHC ->entries();
  if( TOFHC )nhTof = TOFHC ->entries();
  //   if( ACHC )     nhAc     = ACHC ->entries();
  if( WCHC ) nhWC  = WCHC ->entries();
  if( SlHC ) nhSl  = SlHC ->entries();
  
  // G4cout<<"DCHC : true, colID="<<colIdDC<<" TOFHC : true, colID="<<colIdTOF<<G4endl;
  
  //DC    
   G4double pos_res = 0.0*mm; // 0 um
   // G4double pos_res = 0.2*mm; // 200 um
  event.DCNhits = nhDc;
  if( DCHC ){
    for( int i=0; i<nhDc; ++i ){
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
      event.DCXObs[dclayer] = RandGauss::shoot(lx, pos_res);
      event.DCYObs[dclayer] = RandGauss::shoot(ly, pos_res);
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

  //int nhitSl1=0;
  // ~~~~ Bending angle of the S-2S dipole magnet ~~~~
  const double BendingAngle = 70.0; // [deg]
  // ~~~~      ~~~~~   ~~~~~~ ~~~~  ~~~~   ~~~~~   ~~~
  if( SlHC ){
    for( int i=0; i<nhSl; ++i ){
      SlHit *aHit = (*SlHC)[i];
      G4int Sllayer = aHit->GetLayerID();
      G4int nh = event.SlitNh[Sllayer];
      G4double lx = aHit->GetXLocal()/mm;
      G4double ly = aHit->GetYLocal()/mm;
      G4double time = aHit->GetTime()/ns;
      G4ThreeVector pVec = aHit->GetMom();
      G4double path = aHit->GetPath()/mm;
      G4String name = aHit->GetDecayParticleName(0);
      G4double mom = pVec.mag();
      event.SlituDeg[Sllayer] = atan(pVec.y()/pVec.x())*Rad2Deg;
      event.SlitvDeg[Sllayer] = atan(pVec.z()/pVec.x())*Rad2Deg;
      event.SlitX[Sllayer][nh]  = lx;
      // ------ High momentum: Large x (Toshi, 19Mar2015) ------
      //G4cout<<"layer ="<<Sllayer<<" lx= "<<lx<<G4endl;
      event.SlitX[Sllayer][nh] = -1.0 * event.SlitX[Sllayer][nh];
      event.SlituDeg[Sllayer] = -1.0 * event.SlituDeg[Sllayer];
      if(Sllayer>4){ // ~~~~~ After the dipole magnet ~~~~~
	event.SlituDeg[Sllayer] = event.SlituDeg[Sllayer] + BendingAngle;
      }
      // -------------------------------------------------------
      event.SlitY[Sllayer]  = ly;
      event.Slitt[Sllayer]  = time;
      event.SlitMom[Sllayer]= mom;
      event.Slitp[Sllayer]  = path;
      event.SlitNh[Sllayer]++;
      if(name=="proton"){ 
	event.SlitNP[Sllayer]++; 
	event.SlitF[Sllayer]=1;
      }
      else if(name=="kaon+"){
	event.SlitNK[Sllayer]++; 
	event.SlitF[Sllayer]=2;
      }
      else if(name=="pi+"){ 
	event.SlitNPi[Sllayer]++;
	  event.SlitF[Sllayer]=3;
      }
      else if(name=="anti_proton"){
	event.SlitNPi[Sllayer]++;
	event.SlitF[Sllayer]=4;
      }
      else if(name=="kaon-"){
	event.SlitNPi[Sllayer]++; 
	event.SlitF[Sllayer]=5;
      }
      else if(name=="pi-"){
	event.SlitNPi[Sllayer]++;
	event.SlitF[Sllayer]=6;
      }
      else if(name=="e+"){
	event.SlitNPi[Sllayer]++;
	event.SlitF[Sllayer]=7;
      }
      else if(name=="e-"){
	event.SlitNPi[Sllayer]++;
	event.SlitF[Sllayer]=8;
      }
      else if(name=="mu+"){
	event.SlitNPi[Sllayer]++;
	event.SlitF[Sllayer]=9;
      }
      else if(name=="mu-"){
	event.SlitNPi[Sllayer]++;
	event.SlitF[Sllayer]=10;
      }
      else if(name=="neutron"){
	event.SlitNPi[Sllayer]++;
	event.SlitF[Sllayer]=11;
      }
      else if(name=="gamma"){
	event.SlitNPi[Sllayer]++;
	event.SlitF[Sllayer]=12;
      }
      else{
	//if(Sllayer==7) G4cerr<<"Enigma Particle Name="<<name<<G4endl;
      }
    }
  }
  G4bool SlitFlag = false;
  if ( event.SlitX[0][0]>-1000.0 && // Before Q1
       //event.SlitX[1]>-1000.0 && 
       //event.SlitX[2]>-1000.0 && 
       //event.SlitX[3]>-1000.0 && 
       //event.SlitX[4]>-1000.0 && 
       event.SlitX[5][0]>-1000.0 && //After Dipole
       event.SlitX[6][0]>-1000.0 && //After Dipole end guard
       event.SlitX[7][0]>-1000.0 && // After Chamber
       event.SlitX[8][0]>-1000.0 && // After TOF
       event.SlitX[9][0]>-1000.0 && // After AC
       event.SlitX[10][0]>-1000.0 // After WC
       ){
    SlitFlag = true;
  }
  else{ SlitFlag=false; }
  
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
  
  // ~~~~~~~~~~~~~~~~ TOF detector ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  event.TOFNhits = (int)nhTof;
  G4double t_res = 0.090*ns; // Sigma = 90 ps (From cosmic-ray test)
  G4bool TOFTrig = false;
  if( TOFHC ){
    for( int i=0; i<nhTof; ++i ){
      TOFHit *aHit = (*TOFHC)[i];
      G4int TOFHitSeg = aHit->GetLayerID(); // Hit Segment
      double TOFtime = aHit->GetTime()/ns;  // Hit timing
      double TOFdE   = aHit->GetEdep()/MeV; // Energy deposite
      event.toftime[TOFHitSeg] = event.toftime[TOFHitSeg] + TOFtime;
      event.toftime_reso[TOFHitSeg] = event.toftime_reso[TOFHitSeg] + RandGauss::shoot(TOFtime, t_res);
      event.tofdE[TOFHitSeg]   = event.tofdE[TOFHitSeg] + TOFdE;
      event.tofn[TOFHitSeg]++;
      //G4cout<<"TofSeg"<<TOFHitSeg<<" event.TOFC["<<TOFHitSeg<<"]"<<event.TOFt[TOFHitSeg]<<G4endl;
    }
    
    for(int i=0 ; i<18 ; i++){
      if(event.toftime[i]>0.0){
	// ----- Mean value of the timing -----
	event.toftime[i] = event.toftime[i] / event.tofn[i];  // w/o resolution
	event.toftime_reso[i] = event.toftime_reso[i] / event.tofn[i]; // w/ resolution
	
	//if(event.toftime[i]>10.0 && event.tofdE[i]>0.0){
	if( event.toftime[i]>10.0 ){
	  TOFTrig = true;
	}
	else TOFTrig = false;
	
      }
      else{
	event.toftime[i] = -2222.0;
	event.tofdE[i]   = -2222.0;
      }
    }
  }
  
  event.TOFTrig = TOFTrig;
  event.VDTrig  = SlitFlag;
  
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

   // ~~~~~~~~~~~~~~ Water Cherenkov detecctor  ~~~~~~~~~~~~~~~~~~~~~~~~~
  event.WCNhits = (int)nhWC;
  //G4double t_res = 0.090*ns; // Sigma = 90 ps (From cosmic-ray test)
  G4bool WCTrig = false;
  //TRandom3* randtemp = new TRandom3();
  //double aaatemp=0.0;
  //aaatemp = randtemp->Uniform(1.0,65539.0);
  //rrr = rrr%65539;
  //TRandom3* nperand = new TRandom3(aaatemp);
  double wcnpe1, wcnpe2;
  if( WCHC ){
    for( int i=0; i<nhWC; i++ ){
      WCHit *aHit = (*WCHC)[i];
      G4int WCHitSeg = aHit->GetLayerID(); // Hit Segment
      double WCtime = aHit->GetTime()/ns;  // Hit timing
      double WCdE   = aHit->GetEdep()/MeV; // Energy deposite
      double WCNPE  = aHit->GetNPE();  // N.P.E.
      event.wctime[WCHitSeg] = event.wctime[WCHitSeg] + WCtime;
      //event.wctime_reso[WCHitSeg] = event.wctime_reso[WCHitSeg] + RandGauss::shoot(WCtime, t_res);
      event.wcdE[WCHitSeg]   = event.wcdE[WCHitSeg] + WCdE;
      event.wcnpe[WCHitSeg]  = event.wcnpe[WCHitSeg] + WCNPE;
      //event.wcnpe[WCHitSeg]  = event.wcnpe[WCHitSeg] + nperand->PoissonD(WCNPE);
      event.wcn[WCHitSeg]++;
      //G4cout<<"TofSeg"<<TOFHitSeg<<" event.TOFC["<<TOFHitSeg<<"]"<<event.TOFt[TOFHitSeg]<<G4endl;
    }
    for(int i=0 ; i<12 ; i++){
      if(event.wctime[i]>0.0){
	// ----- Mean value of the timing -----
	event.wctime[i] = event.wctime[i] / event.wcn[i];  // w/o resolution
	//event.wctime_reso[i] = event.wctime_reso[i] / event.wcn[i]; // w/ resolution
	///event.wcnpe[i] = G4Poisson(event.wcnpe[i]) * 2.0; // up and down PMT --> * 2
	/*
	  event.wcnpe[i] = nperand->PoissonD(event.wcnpe[i]) * 2.0;
	*/
	//event.wcnpe[i] = nperand->PoissonD(event.wcnpe[i]);
	wcnpe1 = nperand->PoissonD(event.wcnpe[i]);
	wcnpe2 = nperand->PoissonD(event.wcnpe[i]);
	event.wcnpe[i] = wcnpe1+wcnpe2;
	//event.wcnpe[i] = nperand->Gaus(event.wcnpe[i],10.0); // 10:width from cosmic-ray data
	//event.wcnpe[i] = event.wcnpe[i] * 2.0;
	//if(event.wctime[i]>10.0 && event.wcdE[i]>0.0){
	if(event.wctime[i]>10.0 &&
	   event.wcnpe[i]>0.0){// && 
	   //event.wcnpe[i]>30.0){
	  //event.wcnpe[i] = nperand->PoissonD(event.wcnpe[i]);
	  WCTrig = true;
	  if(i<6){
	    event.wctime1[i] = event.wctime[i];
	    event.wcnpe1[i]  = event.wcnpe[i];
	    event.wcdE1[i]   = event.wcdE[i];
	    event.wcnpe1[i] =  RandGauss::shoot(event.wcnpe1[i], event.wcnpe1[i]*0.06); // +/-6% (Gogami)
	    //event.wcnpe1[i] = event.wcnpe1[i] 
	    //+ (event.wcnpe1[i] * (G4UniformRand()-0.5)*2.0 * 0.14);  // +/- 14% (Takenaka value)
	    //+ (event.wcnpe1[i] * (G4UniformRand()-0.5)*2.0 * 0.2);  // +/- 20%
	    //+ (event.wcnpe1[i] * (G4UniformRand()-0.5)*2.0 * 1.0);  // +/- 100%
	    //if (i==5)G4cout << event.wcnpe1[i] << G4endl; // for check
	  }
	  else{
	    event.wctime2[i-6] = event.wctime[i];
	    event.wcnpe2[i-6]  = event.wcnpe[i];
	    event.wcdE2[i-6]   = event.wcdE[i];
	    event.wcnpe2[i-6] =  RandGauss::shoot(event.wcnpe2[i-6], event.wcnpe2[i-6]*0.06); // +/-6% (Gogami)
	    //event.wcnpe2[i-6] = event.wcnpe2[i-6] 
	    //+ (event.wcnpe2[i-6] * (G4UniformRand()-0.5)*2.0 * 0.14); // +/- 14% (Takenaka value)
	    //+ (event.wcnpe2[i-6] * (G4UniformRand()-0.5)*2.0 * 0.2); // +/- 20%
	    //+ (event.wcnpe2[i-6] * (G4UniformRand()-0.5)*2.0 * 1.0); // +/- 100%
	  }
	}
	else WCTrig = false;
      }
      else{
	event.wctime[i] = -2222.0;
	event.wcdE[i]   = -2222.0;
	event.wcnpe[i]  = -2222.0;
      }
    }
  }
  event.WCTrig = WCTrig;
  
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
//    //    tree->Fill();
//    //if( DataFile_.is_open() ){
//    if(0){
//      PrintHitsInformation( anEvent, DataFile_ );
//    }
//  }

  //  G4cout<<"EndOf S2SAnalysis EndOfEvent"<<G4endl;
  //if(tree->GetEntries()>1&&tree->GetEntries()%100==0)
  //if(1)
  // anafile->Write();
  //if(TOFTrig==true){
  
  if(SlitFlag==true){
    //tree->Fill();
  }
  tree->Fill();

}

void S2SAnalysis::SaveFile( void ) const
{
  if( fActive_ ){
  gFile->Write();
  }
}

void S2SAnalysis::Terminate( void ) const
{
  if( fActive_ ){
    anafile->Write();
    anafile->Close();
    G4cout << "[S2SAnalysis] Terminate() done" << G4endl;
  }
}

void S2SAnalysis::ShowStatus( void ) const
{
  G4cout << "Analyzer Status\n"
	 << "  File : " << filename_ << "\n"
	 << "  S2SAnalysis : " << fActive_ << "\n"
	 << G4endl;
}

void S2SAnalysis::InitializeEvent( void )
{
  for(int i=0;i<NumDC;i++){
    event.DCXObs[i] = -999.;
    event.DCYObs[i] = -999.;
    event.DCX[i] = -999.;
    event.DCY[i] = -999.;
    event.DCt[i] = -999.;
    event.DCp[i] = -999.;
    event.DCNh[i] = 0;
    for( int ihit = 0; ihit < MaxHits; ihit++ ){
      event.DCgPosx[i][ihit] = -2222.0;
      event.DCgPosy[i][ihit] = -2222.0;
      event.DCgPosz[i][ihit] = -2222.0;
      event.DCde[i][ihit] = -2222.0;
    }
  }
//   for(int i=0;i<NumTOFSeg;i++){
//     event.TOFt[i] = -999.;
//     event.TOFtObs[i] = -999.;
//   }
  // event.TOFAll = -999.;
//  for(int i=0;i<8;i++){
//    event.SlitX[i]    = -2222.0;
//    event.SlitY[i]    = -2222.0;
//    event.SlituDeg[i] = -2222.0;
//    event.SlitvDeg[i] = -2222.0;
//    event.Slitt[i]    = -2222.0;
//    event.SlitMom[i]  = -2222.0;
//    event.Slitp[i]    = -2222.0;
//    event.SlitNh[i]   = 0;
//    event.SlitNP[i]   = 0;
//    event.SlitNK[i]   = 0;
//    event.SlitNPi[i]  = 0;
//    event.SlitF[i]    = 0;
//  }
  for(int i=0 ; i<18 ; i++){
    // ------- Virtual detectors (slits) ------ 
    if(i<11){
      for( int ihit = 0; ihit < MaxHits; ihit++ ) event.SlitX[i][ihit]    = -2222.0;
      event.SlitY[i]    = -2222.0;
      event.SlituDeg[i] = -2222.0;
      event.SlitvDeg[i] = -2222.0;
      event.Slitt[i]    = -2222.0;
      event.SlitMom[i]  = -2222.0;
      event.Slitp[i]    = -2222.0;
      event.SlitNh[i]   = 0;
      event.SlitNP[i]   = 0;
      event.SlitNK[i]   = 0;
      event.SlitNPi[i]  = 0;
      event.SlitF[i]    = 0;
    }
    
    // ------- TOF detector ------ 
    //event.tofco[i]        = -2222.0;
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

  event.DCNhits = -999;
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
  // event.ACX = -999.;
  // event.ACY = -999.;
  // event.ACt = -999.;
  // for(int i=0;i<NumWCSeg;i++){
  //   event.WCt[i] = -999.;
  //   event.WCp[i] = -999.;
  // }
  //  G4cout<<"S2SAnalysis InitializeEvent"<<G4endl;
}

//void S2SAnalysis::DefineHistograms()
void S2SAnalysis::DefineTree(){
  //G4cout<<"S2SAnalysis DefineHistograms"<<G4endl;
  anafile = new TFile(filename_,"recreate");
  TTree* tree = new TTree("s2s_g","S-2S simulation");
  fActive_=true;

  tree->Branch("x0",&event.x0In, "x0/D");
  tree->Branch("y0",&event.y0In, "y0/D");
  tree->Branch("z0",&event.z0In, "z0/D");
//   tree->Branch("u0In",&event.u0In, "u0In/D");
//   tree->Branch("v0In",&event.v0In, "v0In/D");
  tree->Branch("xp0",&event.uDeg, "xp0/D");
  tree->Branch("yp0",&event.vDeg, "yp0/D");
  tree->Branch("phi0",&event.phi0, "phi0/D");
  tree->Branch("theta0",&event.theta0, "theta0/D");
  tree->Branch("p0",&event.p0,   "p0/D");
  
  //  tree->Branch("t0",&event.t0,   "t0/D");
  // tree->Branch("Id",&event.Id, "Id/I");
//   tree->Branch("nP",&event.nP,"nP/I");
//   tree->Branch("nN",&event.nN,"nN/I");
//   tree->Branch("nL",&event.nL,"nL/I");
//   tree->Branch("nSm",&event.nSm,"nSm/I");
//   tree->Branch("nSz",&event.nSz,"nSz/I");
//   tree->Branch("nSp",&event.nSp,"nSp/I");
//   tree->Branch("nXm",&event.nXm,"nXm/I");
//   tree->Branch("nXz",&event.nXz,"nXz/I");
//   tree->Branch("nXsm",&event.nXsm,"nXsm/I");
//   tree->Branch("nXsz",&event.nXsz,"nXsz/I");
//   tree->Branch("nPim",&event.nPim,"nPim/I");
//   tree->Branch("nPiz",&event.nPiz,"nPiz/I");
//   tree->Branch("nPip",&event.nPip,"nPip/I");
//   tree->Branch("nKm",&event.nKm,"nKm/I");
//   tree->Branch("nKp",&event.nKp,"nKp/I");

  tree->Branch("DCXObs",event.DCXObs, "DCXObs[26]/D");
  tree->Branch("DCYObs",event.DCYObs, "DCYObs[26]/D");
  tree->Branch("DCX",event.DCX, "DCX[26]/D");
  tree->Branch("DCY",event.DCY, "DCY[26]/D");
  tree->Branch("DCt",event.DCt, "DCt[26]/D");
  tree->Branch("DCp",event.DCp, "DCp[26]/D");
  tree->Branch("DCNh",event.DCNh, "DCNh[26]/I");
  tree->Branch("DCgPosx",event.DCgPosx, Form("DCgPosx[26][%d]/D", MaxHits));
  tree->Branch("DCgPosy",event.DCgPosy, Form("DCgPosy[26][%d]/D", MaxHits));
  tree->Branch("DCgPosz",event.DCgPosz, Form("DCgPosz[26][%d]/D", MaxHits));
  tree->Branch("DCde",   event.DCde, Form("DCde[26][%d]/D", MaxHits));
  tree->Branch("DCNhits",&event.DCNhits, "DCNhits/I");
  tree->Branch("DC1Hit", &event.DC1Hit, "DC1Hit/I");
  tree->Branch("DC2Hit", &event.DC2Hit, "DC2Hit/I");
  tree->Branch("DC3Hit", &event.DC3Hit, "DC3Hit/I");
  tree->Branch("DC4Hit", &event.DC4Hit, "DC4Hit/I");
  tree->Branch("DC5Hit", &event.DC5Hit, "DC5Hit/I");
  tree->Branch("vdxp", event.SlituDeg, "vdxp[11]/D");
  tree->Branch("vdyp", event.SlitvDeg, "vdyp[11]/D");
  /*
    tree->Branch("SlitX", event.SlitX, "SlitX[8]/D");
    tree->Branch("SlitY", event.SlitY, "SlitY[8]/D");
    tree->Branch("Slitt", event.Slitt, "Slitt[8]/D");
    tree->Branch("SlitMom", event.SlitMom, "SlitMom[8]/D");
    tree->Branch("Slitp", event.Slitp, "Slitp[8]/D");
    tree->Branch("SlitNh", event.SlitNh, "SlitNh[8]/I");
    tree->Branch("SlitNP", event.SlitNP, "SlitNP[8]/I");
    tree->Branch("SlitNK", event.SlitNK, "SlitNK[8]/I");
    tree->Branch("SlitNPi", event.SlitNPi, "SlitNPi[8]/I");
    tree->Branch("SlitF", event.SlitF, "SlitF[8]/I");
  */
  
  // Order -->
  // vd8 | TOF | vd9 | AC | vd10 | WC | vd11 (vd[10]) | SDC1 
  tree->Branch("vdx",   event.SlitX, Form("vdx[11][%d]/D", MaxHits));
  tree->Branch("vdy",   event.SlitY, "vdy[11]/D");
  tree->Branch("vdtime",event.Slitt, "vdtime[11]/D");
  tree->Branch("vdmom", event.SlitMom, "vdmom[11]/D");
  tree->Branch("vdpath",event.Slitp, "vdpath[11]/D");
  tree->Branch("vdNh",  event.SlitNh, "vdNh[11]/I");
  tree->Branch("vdNP",  event.SlitNP, "vdNP[11]/I");
  tree->Branch("vdNK",  event.SlitNK, "vdNK[11]/I");
  tree->Branch("vdNPi", event.SlitNPi, "vdNPi[11]/I");
  tree->Branch("vdPID", event.SlitF, "vdPID[11]/I");
  // tree->Branch("TOFAll", &event.TOFAll, "TOFAll/D");
  // tree->Branch("TOFt",event.TOFt,   "TOFt[17]/D");
  // tree->Branch("TOFtObs",event.TOFtObs,   "TOFtObs[17]/D");
  //tree->Branch("tofco",  &event.tofco,   "tofco[17]/I");
  
  tree->Branch("tofnhits",&event.TOFNhits, "tofnhits/I");
  tree->Branch("toftime", &event.toftime,  "toftime[18]/D");
  tree->Branch("tofdE",   &event.tofdE,    "tofdE[18]/D");
  
  tree->Branch("wcnhits",&event.WCNhits,  "wcnhits/I");
  //tree->Branch("wctime", &event.wctime,   "wctime[12]/D");
  //tree->Branch("wcdE",   &event.wcdE,     "wcdE[12]/D");
  //tree->Branch("wcnpe",  &event.wcnpe,    "wcnpe[12]/D");
  tree->Branch("wctime1", &event.wctime1,   "wctime1[6]/D");
  tree->Branch("wcdE1",   &event.wcdE1,     "wcdE1[6]/D");
  tree->Branch("wcnpe1",  &event.wcnpe1,    "wcnpe1[6]/D");
  tree->Branch("wctime2", &event.wctime2,   "wctime2[6]/D");
  tree->Branch("wcdE2",   &event.wcdE2,     "wcdE2[6]/D");
  tree->Branch("wcnpe2",  &event.wcnpe2,    "wcnpe2[6]/D");
  
  tree->Branch("Q1Trig", &event.Q1Trig,  "Q1Trig/B");
  tree->Branch("Q2Trig", &event.Q2Trig,  "Q2Trig/B");
  tree->Branch("VDTrig", &event.VDTrig,  "VDTrig/B");
  tree->Branch("TOFTrig",&event.TOFTrig, "TOFTrig/B");
  tree->Branch("WCTrig", &event.WCTrig,  "WCTrig/B");
  // tree->Branch("TOFHit", &event.TOFHit, "TOFHit/I");
  // tree->Branch("ACX", &event.ACX, "ACX/D");
  // tree->Branch("ACY", &event.ACY, "ACY/D");
  // tree->Branch("ACt", &event.ACt, "ACt/D");
  // tree->Branch("ACp", &event.ACp, "ACp/D");
  // tree->Branch("ACNhits", &event.ACNhits, "ACNhits/I");
  // tree->Branch("ACHit", &event.ACHit, "ACHit/I");
  // tree->Branch("WCt", event.WCt,  "WCt[16]/D");
  // tree->Branch("WCp", event.WCp,  "WCp[16]/D");
  // tree->Branch("WCNhits", &event.WCNhits, "WCNhits/I");
  // tree->Branch("WCHit", &event.WCHit, "WCHit/I");
  //G4cout<<"tree and Branch is defined"<<G4endl;

}

void S2SAnalysis::
PrintHitsInformation( const G4Event *anEvent, 
		      std::ostream &ost ) const
{
  //   G4cout<<"PrintHits is called"<<G4endl;
  
//  ConfMan *confMan = ConfMan::GetConfManager();
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

    double lp = lx*cos(tilt*Deg2Rad) + ly*sin(tilt*Deg2Rad);
    lp += RandGauss::shoot(0., reso);
    
    ost << std::setw(12) << layer-100 
  	<< std::setw(12) << lp;
  }

  G4double nhTOF=0;
  TOFHitsCollection    *TOFHC;
  G4int colIdTOF = SDMan->GetCollectionID( "TOFSD"/*"DCCollection"*/ );
  TOFHC =  dynamic_cast<TOFHitsCollection *>(  HCE->GetHC( colIdTOF ) );
  if( TOFHC )     nhTOF     = TOFHC ->entries();

  for( int i=0; i<nhTOF; ++i ){
    TOFHit *aHit = (*TOFHC)[i];
    double time = aHit->GetTime(); 
    double de = aHit->GetEdep();
    if(i==nhTOF-1){
      ost << std::setw(12) << 61
	  << std::setw(12) << time;
      de=1.0;
      ost << std::setw(12) << 62 <<" "
	  << std::setw(12) << de;
    }
  }

  ost << std::setw(5) << -1 << std::endl;

  ost.flags( oldFlags );
  ost.precision( preSiz );

  ost << std::endl;

}

void S2SAnalysis::SetDataFile( const char *datafile )
{
  DataFile_.open( datafile );
}
