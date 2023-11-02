/*
  SksTracking.cc
*/

#include "FieldMan.hh"
#include "DCLocalTrack.hh"
#include "DCHit.hh"
#include "DCConfMan.hh"
#include "SksTrack.hh"
#include "TrackHit.hh"
#include "ThreeVector.hh"
#include "LorentzVector.hh"
#include "Kinematics.hh"
#include "YdScatDetectorConst.hh"
#include "RadDeg.hh"
#include "DCAnalyzer.hh" 
#include "CylindricalAna.hh" 

#include "RootHelper.hh"

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <string>
#include <vector>
#include <signal.h>


//#define Interactive 1

////////////////Physics Parameters///////////////////
const double AtomicMassUnit  = 0.93149432;
const double PionMass        = 0.1395701;
const double KaonMass        = 0.493677;
const double ProtonMass      = 0.93827200;
const double NeutronMass     = 0.93956563;
const double LambdaMass      = 1.115648;
const double SigmaMinusMass  = 1.197449;
//2H
const double DeuteronMass     = 2.*AtomicMassUnit+0.01313672;
//4He
const double Nucl4HeMass     = 4.0*AtomicMassUnit+0.0024249;
const double Delta3He        = 0.0149312;
const double CoreMass3He     = 3.*AtomicMassUnit+Delta3He;
//7Li
const double Nucl7LiMass     = 7.0*AtomicMassUnit+0.014908;
const double Delta6Li        = 0.014086;
const double CoreMass6Li     = 6.*AtomicMassUnit+Delta6Li;
//9Be
const double Nucl9BeMass     = 9.*AtomicMassUnit+0.011347;
const double Delta8Be        = 0.0049417;
const double CoreMass8Be     = 8.*AtomicMassUnit+Delta8Be;
//12C
const double Nucl12CMass     = 12.*AtomicMassUnit;
const double Delta11C        = 0.0106502;
const double CoreMassC11     = 11.*AtomicMassUnit+Delta11C;
//16O
const double Nucl16OMass     = 16.*AtomicMassUnit-0.004736998;
const double Delta15O        = 0.0028556;
const double CoreMass15O     = 15.*AtomicMassUnit+Delta15O;
//10B --> 10LLi
const double Nucl10BMass     = 10.*AtomicMassUnit+0.012050;
const double Delta9Li        = 0.024954;
const double CoreMass9Li     = 9.*AtomicMassUnit+Delta15O;
////////////////////////////////////////////////////


bool ProcessOneEvent( std::ifstream & );
void DefineHistograms( const char * );
void InitializeEvent( void );

const double MaxChiSqr = 100.;
const double PIni = 1.0;

#ifndef MaxHits 
#define MaxHits 30
#endif
#ifndef MaxHits2 
#define MaxHits2 60
#endif

struct Event{
  //Primary for Yd scattering
  double thetaMeson;
  double phiMeson;
  double thetaMesonCM;
  double phiMesonCM;
  double thetaScatHypCM;
  double phiScatHypCM;
  double momVectorScatMeson[3];
  double momVectorHypBeam[3];
  double momVectorHypScat[3];
  double momVectorProtonScat[3];
  double momVectorDecayPi[3];
  double momVectorDecayNucleon[3];
  double momScatMeson;
  double momHypBeam;
  double momHypScat;
  double momProtonScat;
  double momDecayPi;
  double momDecayNucleon;
  double primaryVertex[3];
  double scatPos0[3];
  double NNscatPos[3];
  double PiNscatPos[3];
  double decayPos[3];
  int    decayFlag;
  int    scatFlag;
  int    scatTarget;
  int    NNscatFlag;
  int    NNscatTarget;
  int    PiNscatFlag;
  int    PiNscatTarget;

  double vertex1[3];
  double theta1;
  double cdist1;

  double u0BcOut;
  double v0BcOut;
  double u0SdcIn;
  double v0SdcIn;

  double vertex2[3];
  double theta2;
  double cdist2;

  double vertex3[3];
  double theta3;
  double cdist3;

  double vertex4[3];
  double theta4;
  double cdist4;

  double priMomVec[3];
  double priPosVec[3];

  double priMomSksInVec[3];
  double priPosSksInVec[3];

  double vertex5[3];
  double theta5;
  double cdist5;
  double chisqrSksIn;

  double vertexScat[3];
  double scatVec[3];
  double scatPos[3];
  double hypBeamVec[3];
  double cdistScat;
  double thetaScat;
  double chisqrXY;
  double chisqrZ;
  double scatMomCal;
  double scatEkinCal;
  double thetaScatCM;

  double vertexDecayPi[3];
  double decayPiVec[3];
  double decayPiPos[3];
  double cdistDecayPi;
  double thetaDecayPi;
  double chisqrDecayPiXY;
  double chisqrDecayPiZ;
  double decayPiMomCal;	     
  double decayPiMomVecCal[3];
  double decayNMomCal;	     
  double decayNMomVecCal[3];

  double vertexNpScat[3];
  double cdistNpScat;
  double thetaNpScat;
  double scatNpMomCal;
  double scatNpEkinCal;
  double thetaNpScatCM;

  //SksTracking
  int Hits;
  double chisqr;

  double P;
  double Theta;
  double Phi;
  double X0;
  double Y0;
  double U0;
  double V0;
  double path;
  double Pb;
  double Ub;
  double Vb;

  double dP;
  double dTheta;
  double dPhi;
  double dX;
  double dY;
  double dU;
  double dV;

  double m2;
  double RecoilMom;
  double beta;
  double Ekin;
  double MissMass;
  double BE;
  double BE_SMF;

  int DC1Hits;
  int DC2Hits;
  int DC3Hits;
  int DC4Hits;
  double DC1Posl[MaxHits];
  double DC1Resd[MaxHits];
  double DC2Posl[MaxHits];
  double DC2Resd[MaxHits];
  double DC3Posl[MaxHits];
  double DC3Resd[MaxHits];
  double DC4Posl[MaxHits];
  double DC4Resd[MaxHits];

  int tofHits;
  int smfHits;
  double tofTime[MaxHits];
  double tofEdep[MaxHits];
  double smfTime[MaxHits];
  double smfEdep[MaxHits];

  //Fiber1-4
  int    FiberHits[FiberLayerNum+1];
  int    FiberSeg[FiberLayerNum+1][MaxHits2];
  double FiberTime[FiberLayerNum+1][MaxHits2];
  double FiberEdep[FiberLayerNum+1][MaxHits2];
  int    FiberPID[FiberLayerNum+1][MaxHits2];

  double fClEdepForScat[FiberLayerNum+1];
  double fClEdepNormForScat[FiberLayerNum+1];
  double fClPathLengthForScat[FiberLayerNum+1];

  double fClEdepForPi[FiberLayerNum+1];
  double fClEdepNormForPi[FiberLayerNum+1];
  double fClPathLengthForPi[FiberLayerNum+1];

  //Crystal
  int    CrystalHits[CrystalLayerNum+1];
  int    CrystalSeg[CrystalLayerNum+1][MaxHits2];
  double CrystalTime[CrystalLayerNum+1][MaxHits2];
  double CrystalEdep[CrystalLayerNum+1][MaxHits2];
  int    CrystalPID[CrystalLayerNum+1][MaxHits2];

  double cClEdepForScat;
  double cClEdepForPi;
  double cRcEdepForScat;
  double cRcEdepForPi;

};
static Event event;

const int MinHitForPi = 6;
const int MinHitForP = 4;

void closeFile(int sig)
{
  if(gFile)
    {
      gFile->Write();
      gFile->Close();
    }
}


int main( int argc, char **argv )
{
  char *inFile;
  if( argc<5 ) {
    std::cerr << "Usage: " << argv[0] 
	      << " confFile fieldMap RootFile dataFile"
	      << std::endl;
    exit(-1);
  }

  DCConfMan::GetInstance().Initialize(argv[1]);
  FieldMan::GetInstance().Initialize(argv[2]);
  DefineHistograms(argv[3]);

  signal(SIGINT,closeFile);
  
  if(argc>=5){
    int idt=4;
    while( idt<argc ){
      inFile = argv[idt++];
      std::ifstream InputData( inFile );
      std::cout<<"***************************************************"<<std::endl;
      std::cout<<"Data File "<< idt-4 << std::endl;
      while( ProcessOneEvent( InputData ) ){
#ifdef Interactive 
	int it; 
	std::cout << "# ";
	std::cin >> it;
	if( it<0 ) break;
#endif
      }
    }
  }

  std::cout<<"*************************************"<<std::endl;
  gFile->Write();
  gFile->Close();
  std::cout<<"*************************************"<<std::endl;

  return 0;
}

bool ProcessOneEvent( std::ifstream &In )
{
  const std::string funcname = "ProcessEvent";

  static int events=0;
  events++;
  if (events%1000==0)
    std::cout<<"Events:"<< events<<std::endl;

  //Tree
  TTree *tree = dynamic_cast<TTree *>(gFile->Get("tree"));
  InitializeEvent();

  if( !In ) return false;

  In >> event.thetaMeson >> event.phiMeson >> event.thetaMesonCM >> event.phiMesonCM >> event.thetaScatHypCM >> event.phiScatHypCM;
  In >>  event.momVectorScatMeson[0] >> event.momVectorScatMeson[1] >> event.momVectorScatMeson[2];
  In >>  event.momVectorHypBeam[0] >> event.momVectorHypBeam[1] >> event.momVectorHypBeam[2];
  In >>  event.momVectorHypScat[0] >> event.momVectorHypScat[1] >> event.momVectorHypScat[2];
  In >>  event.momVectorProtonScat[0] >> event.momVectorProtonScat[1] >> event.momVectorProtonScat[2];
  In >>  event.momVectorDecayPi[0] >> event.momVectorDecayPi[1] >> event.momVectorDecayPi[2];
  In >>  event.momVectorDecayNucleon[0] >> event.momVectorDecayNucleon[1] >> event.momVectorDecayNucleon[2];
  In >>  event.primaryVertex[0] >> event.primaryVertex[1] >> event.primaryVertex[2];
  In >>  event.scatPos0[0] >> event.scatPos0[1] >> event.scatPos0[2];
  In >>  event.NNscatPos[0] >> event.NNscatPos[1] >> event.NNscatPos[2];
  In >>  event.PiNscatPos[0] >> event.PiNscatPos[1] >> event.PiNscatPos[2];
  In >>  event.decayPos[0] >> event.decayPos[1] >> event.decayPos[2];
  In >>  event.decayFlag >> event.scatFlag >> event.scatTarget >> event.NNscatFlag >> event.NNscatTarget >> event.PiNscatFlag >> event.PiNscatTarget;

  event.momScatMeson = 
    sqrt(event.momVectorScatMeson[0]*event.momVectorScatMeson[0] + 
         event.momVectorScatMeson[1]*event.momVectorScatMeson[1] + 
         event.momVectorScatMeson[2]*event.momVectorScatMeson[2]);
  event.momHypBeam = 
    sqrt(event.momVectorHypBeam[0]*event.momVectorHypBeam[0] + 
         event.momVectorHypBeam[1]*event.momVectorHypBeam[1] + 
         event.momVectorHypBeam[2]*event.momVectorHypBeam[2]);
  event.momHypScat = 
    sqrt(event.momVectorHypScat[0]*event.momVectorHypScat[0] + 
         event.momVectorHypScat[1]*event.momVectorHypScat[1] + 
         event.momVectorHypScat[2]*event.momVectorHypScat[2]);
  event.momProtonScat = 
    sqrt(event.momVectorProtonScat[0]*event.momVectorProtonScat[0] + 
         event.momVectorProtonScat[1]*event.momVectorProtonScat[1] + 
         event.momVectorProtonScat[2]*event.momVectorProtonScat[2]);
  event.momDecayPi = 
    sqrt(event.momVectorDecayPi[0]*event.momVectorDecayPi[0] + 
         event.momVectorDecayPi[1]*event.momVectorDecayPi[1] + 
         event.momVectorDecayPi[2]*event.momVectorDecayPi[2]);
  event.momDecayNucleon = 
    sqrt(event.momVectorDecayNucleon[0]*event.momVectorDecayNucleon[0] + 
         event.momVectorDecayNucleon[1]*event.momVectorDecayNucleon[1] + 
         event.momVectorDecayNucleon[2]*event.momVectorDecayNucleon[2]);

#if 0
  std::cout << "thetaMeson = " << event.thetaMeson 
	    << " phiMeson = " << event.phiMeson 
	    << " thetaMesonCM = " << event.thetaMesonCM 
	    << " phiMesonCM = "   << event.phiMesonCM 
	    << " thetaScatHypCM = " << event.thetaScatHypCM 
	    << " phiScatHypCM = " << event.phiScatHypCM << std::endl;
  
  std::cout << "momVectorScatmeson = ( " 
	    << event.momVectorScatMeson[0] << ", " 
	    << event.momVectorScatMeson[1] << ", " 
	    << event.momVectorScatMeson[2] << ")" << std::endl ;
  std::cout << "momVectorHypBeam = ( " 
	    << event.momVectorHypBeam[0] << ",  " 
	    << event.momVectorHypBeam[1] << ", " 
	    << event.momVectorHypBeam[2] << ")"  << std::endl;
  std::cout << "momVectorHypScat = ( " 
	    << event.momVectorHypScat[0] << ", " 
	    << event.momVectorHypScat[1] << ", " 
	    << event.momVectorHypScat[2] << ")"  << std::endl;
  std::cout << "momVectorProtonScat = ( " 
	    << event.momVectorProtonScat[0] << ", " 
	    << event.momVectorProtonScat[1] << ", " 
	    << event.momVectorProtonScat[2] << ")" << std::endl; ;
  std::cout << "momVectorDecayPi = ( " 
	    << event.momVectorDecayPi[0] << ", " 
	    << event.momVectorDecayPi[1] << ", " 
	    << event.momVectorDecayPi[2] << ")" << std::endl ;
  std::cout << "momVectorDecayNucleon = ( " 
	    << event.momVectorDecayNucleon[0] << ", " 
	    << event.momVectorDecayNucleon[1] << ", " 
	    << event.momVectorDecayNucleon[2] << ")"  << std::endl;
  std::cout << "primaryVertex = ( " 
	    << event.primaryVertex[0] << ", " 
	    << event.primaryVertex[1] << ", " 
	    << event.primaryVertex[2] << ")"  << std::endl;
  std::cout << "scatPos = ( " 
	    << event.scatPos0[0] << ", " 
	    << event.scatPos0[1] << ", " 
	    << event.scatPos0[2] << ")"  << std::endl;
  std::cout << "NNscatPos = ( " 
	    << event.NNscatPos[0] << ", " 
	    << event.NNscatPos[1] << ", " 
	    << event.NNscatPos[2] << ")"  << std::endl;
  std::cout << "decayPos = ( " 
	    << event.decayPos[0] << ", " 
	    << event.decayPos[1] << ", " 
	    << event.decayPos[2] << ")"  << std::endl;
  
  std::cout << "decayFlag = " << event.decayFlag 
	    << ", scatFlag =  " << event.scatFlag 
	    << ", scatTarget =  " << event.scatTarget 
	    << ", NNscatFlag =  " << event.NNscatFlag 
	    << ", NNscatTarget = " << event.NNscatTarget << std::endl;
#endif
  
  double gx,gy,gz,gp,gt,gf,gpb,gub,gvb;
  In >> gx >> gy >> gz >> gp >> gt >> gf >> gpb >> gub >> gvb;

#if 0
  std::cout<< "*********************************************" <<std::endl;
  std::cout<<" "<< gx 
	   <<" "<< gy 
	   <<" "<< gz 
	   <<" "<< gp 
	   <<" "<< gt 
	   <<" "<< gf
	   <<" "<< gpb
	   <<" "<< gub
	   <<" "<< gvb
	   <<std::endl;
  std::cout<< "*********************************************" <<std::endl;
#endif

  DCAnalyzer DCAna;
  CylindricalAna CylAna;
  CylindricalAna CylAnaPi;
  int lnum;
  In >> lnum;
  volatile int lnumold=lnum;
  int tofseg, smfseg;
  double toftime = -1.0, tofedep = -1.0;
  double smftime = -1.0, smfedep = -1.0;
  int tofhits = 0;
  int smfhits = 0;

  while( lnum!=-1 ){
    double pos=0.0;
    int segment=-1, pid=-1;
    double edep, time;
    double x, y, z;

    if (lnum < 150)
      In >> pos >> lnum;
    else if ((lnum>=151 && lnum<=150+FiberLayerNum) || (lnum >= 171 && lnum<=170+CrystalLayerNum))
      In >> segment >> time >> edep >> pid >> lnum;
    else {
      std::cout << "Strange lnum : " << lnum << std::endl;
      return true;
    }

#if 0
      std::cout<< "*********************************************" <<std::endl;
  std::cout<<" "<< pos 
	   <<" "<< lnum 
	   <<" "<< lnumold 
	   <<std::endl;
  std::cout<< "*********************************************" <<std::endl;
#endif
    if( lnumold<60 ){
      DCAna.SetDCHits(lnumold, pos);
    }
    if( lnumold>112 && lnumold <150){
      DCAna.SetDCHits(lnumold, pos);
    }

    //TOF
    if( lnum == 62 && lnumold == 61 ){
      toftime = pos;
      if (tofhits<MaxHits && tofhits>=0) {
	event.tofTime[tofhits] = toftime; 
	tofhits++;
      } else {
	std::cerr << "Too many tofhits " << tofhits << std::endl;
      }
    }
    if( (lnum == 61 && lnumold == 62) ||
	(lnum == -1 && lnumold == 62) ||
	(lnum >= 151 && lnumold == 62)){
      tofedep = pos;
      if (tofhits-1>=0 && tofhits-1<MaxHits) {
	event.tofEdep[tofhits-1] = tofedep;
      } else {
	std::cerr << "Strange tofhits-1 " << tofhits-1 << std::endl;
      } 
    }
    //SMF
    if( lnum == 72 && lnumold == 71 ){
      smftime = pos;
      if (smfhits<MaxHits && smfhits>=0) {
	event.smfTime[smfhits] = smftime; 
	smfhits++;
      }else {
	std::cerr << "Too many smfhits " << smfhits << std::endl;
      } 
    }
    if( (lnum == 71 && lnumold == 72) ||
	(lnum == -1 && lnumold == 72) || 
	(lnum >= 151 && lnumold == 72) ){
      smfedep = pos;
      if (smfhits-1>=0 && smfhits-1<MaxHits) {
	event.smfEdep[smfhits-1] = smfedep;
      } else {
	std::cerr << "Strange smfhits-1 " << smfhits-1 << std::endl;
      }
    }
    // Fiber
    if( lnumold >= 151 && lnumold <= 150+FiberLayerNum){
      int flayer = lnumold-150;
      if (event.FiberHits[flayer] < MaxHits2) {
	event.FiberSeg [flayer][event.FiberHits[flayer]] = segment;
	event.FiberTime[flayer][event.FiberHits[flayer]] = time;
	event.FiberEdep[flayer][event.FiberHits[flayer]] = edep;
	event.FiberPID [flayer][event.FiberHits[flayer]] = pid;
	event.FiberHits[flayer]++;

	if (pid == 10 || pid == 11 || pid == 14)
	  CylAna.SetFiberDataForScatPart(flayer, segment, edep);
	else if (pid == 1 || pid == 2)
	  CylAnaPi.SetFiberDataForScatPart(flayer, segment, edep);
	/*
	  std::cout << "Flayer = " << flayer 
	  << ", Segment = " << segment
	  << ", Edep = " << edep
	  << ", PID = " << pid << std::endl;
	*/
      } else {
	std::cerr << "Too many fiber hits : " << event.FiberHits[flayer]
		  << std::endl;
      }
    }
    // Crystal
    if( lnumold >= 171 && lnumold <= 170+CrystalLayerNum){
      int clayer = lnumold-170;
      if (event.CrystalHits[clayer] < MaxHits2) {
	event.CrystalSeg [clayer][event.CrystalHits[clayer]] = segment;
	event.CrystalTime[clayer][event.CrystalHits[clayer]] = time;
	event.CrystalEdep[clayer][event.CrystalHits[clayer]] = edep;
	event.CrystalPID [clayer][event.CrystalHits[clayer]] = pid;
	event.CrystalHits[clayer]++;
	if ((clayer ==1 || clayer == 2) && (pid == 10 || pid == 11 || pid == 14))
	  CylAna.SetCrystalDataForScatPart(segment, edep);
	else if ((clayer ==3 || clayer == 4) && (pid == 10 || pid == 11 || pid == 14))
	  CylAna.SetRcDataForScatPart(segment, edep);
	else if ((clayer ==1 || clayer == 2) && (pid == 1 || pid == 2))
	  CylAnaPi.SetCrystalDataForScatPart(segment, edep);
	else if ((clayer ==3 || clayer == 4) && (pid == 1 || pid == 2))
	  CylAnaPi.SetRcDataForScatPart(segment, edep);
	/*
	  std::cout << "Crystal "
	  << ", Segment = " << segment
	  << ", Edep = " << edep
	  << ", PID = " << pid << std::endl;
	*/
      } else {
	std::cerr << "Too many crystal hits : " << event.CrystalHits
		  << std::endl;
      }
    }

    lnumold = lnum;
  }
  event.tofHits  = tofhits; 
  event.smfHits  = smfhits; 

#if 0
  std::cout<< "*****************************************************" <<std::endl;
  std::cout<< "Hits="<< tofhits
	   << " Hits="<< smfhits
	   <<std::endl;
  for( int i=0; i<tofhits; i++){
    std::cout<< " T="<<  event.tofTime[i]
	     << " E="<<  event.tofEdep[i] 
	     <<std::endl;
  }
  for( int i=0; i<smfhits; i++){
    std::cout<< " T="<< event.smfTime[i] 
	     << " E="<< event.smfEdep[i] 
	     <<std::endl;
  }
  std::cout<< "*****************************************************" <<std::endl;
#endif

  //HF1( 101, 1. );

  // Analysis of Cylindrical Detector
  int fHitNumPi = CylAnaPi.Evaluate();
  bool flagFitPiXY, flagFitPiZ;

  ThreeVector piDecayMomVec;
  ThreeVector nDecayMomVec;
  ThreeVector VertDecayPi;

  if (fHitNumPi >= MinHitForPi) {
    for (int i=1; i<=FiberLayerNum; i++)
      event.fClEdepForPi[i] = CylAnaPi.GetFiberClusterEdep(i);
    event.cClEdepForPi = CylAnaPi.GetCrystalClusterEdep();
    event.cRcEdepForPi = CylAnaPi.GetRcClusterEdep();

    flagFitPiXY = CylAnaPi.FitXYPlane();
    if (flagFitPiXY)
      flagFitPiZ = CylAnaPi.FitZPlane();

    if (flagFitPiXY && flagFitPiZ) {

      ThreeVector DecayPiPos = CylAnaPi.GetScatPosition();
      ThreeVector DecayPiVec = CylAnaPi.GetScatVector();

      for (int i=1; i<=FiberLayerNum; i++) {
	double edepNorm = CylAnaPi.GetFiberClusterEdepNorm(i);
	event.fClEdepNormForPi[i] = edepNorm;
	double pathLength = CylAnaPi.GetFiberPathLength(i);
	event.fClPathLengthForPi[i] = pathLength;
      }
    }
  }

  int fHitNum = CylAna.Evaluate();
  bool flagFitXY, flagFitZ;
  if (fHitNum >= MinHitForP) {
    for (int i=1; i<=FiberLayerNum; i++)
      event.fClEdepForScat[i] = CylAna.GetFiberClusterEdep(i);
    event.cClEdepForScat = CylAna.GetCrystalClusterEdep();
    event.cRcEdepForScat = CylAna.GetRcClusterEdep();

    flagFitXY = CylAna.FitXYPlane();
    if (flagFitXY)
      flagFitZ = CylAna.FitZPlane();

    if (flagFitXY && flagFitZ) {

      ThreeVector ScatPos = CylAna.GetScatPosition();
      ThreeVector ScatVec = CylAna.GetScatVector();
      double chisqrXY = CylAna.GetChisqrXY();
      double chisqrZ = CylAna.GetChisqrZ();

      event.chisqrXY = chisqrXY;
      event.chisqrZ = chisqrZ;

      for (int i=1; i<=FiberLayerNum; i++) {
	double res = CylAna.GetResidual(i);
	//HF1(610+i, res);
	double edepNorm = CylAna.GetFiberClusterEdepNorm(i);
	event.fClEdepNormForScat[i] = edepNorm;
	double pathLength = CylAna.GetFiberPathLength(i);
	event.fClPathLengthForScat[i] = pathLength;
      }
    }
  }

  if (events>500000)
    return true;


  tree->Fill();

  return true;
}

void InitializeEvent( void )
{
  //Primary for Yd scattering
  event.thetaMeson = -999.9;
  event.phiMeson = -999.9;
  event.thetaMesonCM = -999.9;
  event.phiMesonCM = -999.9;
  event.thetaScatHypCM = -999.9;
  event.phiScatHypCM = -999.9;
  event.momVectorScatMeson[0] = -999.9;
  event.momVectorScatMeson[1] = -999.9;
  event.momVectorScatMeson[2] = -999.9;
  event.momVectorHypBeam[0] = -999.9;
  event.momVectorHypBeam[1] = -999.9;
  event.momVectorHypBeam[2] = -999.9;
  event.momVectorHypScat[0] = -999.9;
  event.momVectorHypScat[1] = -999.9;
  event.momVectorHypScat[2] = -999.9;
  event.momVectorProtonScat[0] = -999.9;
  event.momVectorProtonScat[1] = -999.9;
  event.momVectorProtonScat[2] = -999.9;
  event.momVectorDecayPi[0] = -999.9;
  event.momVectorDecayPi[1] = -999.9;
  event.momVectorDecayPi[2] = -999.9;
  event.momVectorDecayNucleon[0] = -999.9;
  event.momVectorDecayNucleon[1] = -999.9;
  event.momVectorDecayNucleon[2] = -999.9;
  event.momScatMeson = -999.9;
  event.momHypBeam = -999.9;
  event.momHypScat = -999.9;
  event.momProtonScat = -999.9;
  event.momDecayPi = -999.9;
  event.momDecayNucleon = -999.9;
  event.primaryVertex[0] = -999.9;
  event.primaryVertex[1] = -999.9;
  event.primaryVertex[2] = -999.9;
  event.scatPos0[0] = -999.9;
  event.scatPos0[1] = -999.9;
  event.scatPos0[2] = -999.9;
  event.NNscatPos[0] = -999.9;
  event.NNscatPos[1] = -999.9;
  event.NNscatPos[2] = -999.9;
  event.PiNscatPos[0] = -999.9;
  event.PiNscatPos[1] = -999.9;
  event.PiNscatPos[2] = -999.9;
  event.decayPos[0] = -999.9;
  event.decayPos[1] = -999.9;
  event.decayPos[2] = -999.9;
  event.decayFlag = -999;
  event.scatFlag = -999;
  event.scatTarget = -999;
  event.NNscatFlag = -999;
  event.NNscatTarget = -999;
  event.PiNscatFlag = -999;
  event.PiNscatTarget = -999;

  event.vertex1[0] = -999.9;
  event.vertex1[1] = -999.9;
  event.vertex1[2] = -999.9;
  event.theta1 = -999.9;
  event.cdist1 = -999.9;
  event.u0BcOut = -999.9;
  event.v0BcOut = -999.9;
  event.u0SdcIn = -999.9;
  event.v0SdcIn = -999.9;

  event.vertex2[0] = -999.9;
  event.vertex2[1] = -999.9;
  event.vertex2[2] = -999.9;
  event.theta2 = -999.9;
  event.cdist2 = -999.9;

  event.vertex3[0] = -999.9;
  event.vertex3[1] = -999.9;
  event.vertex3[2] = -999.9;

  event.vertex4[0] = -999.9;
  event.vertex4[1] = -999.9;
  event.vertex4[2] = -999.9;
  event.theta4 = -999.9;
  event.cdist4 = -999.9;

  event.priMomVec[0] = -999.9;
  event.priMomVec[1] = -999.9;
  event.priMomVec[2] = -999.9;
  event.priPosVec[0] = -999.9;
  event.priPosVec[1] = -999.9;
  event.priPosVec[2] = -999.9;

  event.priMomSksInVec[0] = -999.9;
  event.priMomSksInVec[1] = -999.9;
  event.priMomSksInVec[2] = -999.9;
  event.priPosSksInVec[0] = -999.9;
  event.priPosSksInVec[1] = -999.9;
  event.priPosSksInVec[2] = -999.9;

  event.vertex5[0] = -999.9;
  event.vertex5[1] = -999.9;
  event.vertex5[2] = -999.9;
  event.theta5 = -999.9;
  event.cdist5 = -999.9;

  event.chisqrSksIn = -999.9;

  event.vertexScat[0] = -999.9;
  event.vertexScat[1] = -999.9;
  event.vertexScat[2] = -999.9;
  
  event.scatVec[0] = -999.9;
  event.scatVec[1] = -999.9;
  event.scatVec[2] = -999.9;
  
  event.scatPos[0] = -999.9;
  event.scatPos[1] = -999.9;
  event.scatPos[2] = -999.9;
  
  event.hypBeamVec[0] = -999.9;
  event.hypBeamVec[1] = -999.9;
  event.hypBeamVec[2] = -999.9;
  
  event.cdistScat = -999.9;
  event.thetaScat = -999.9;
  event.chisqrXY = -999.9;
  event.chisqrZ = -999.9;
  event.scatMomCal = -999.9;
  event.scatEkinCal = -999.9;
  event.thetaScatCM = -999.9;


  event.vertexDecayPi[0] = -999.9;
  event.vertexDecayPi[1] = -999.9;
  event.vertexDecayPi[2] = -999.9;
  
  event.decayPiVec[0] = -999.9;
  event.decayPiVec[1] = -999.9;
  event.decayPiVec[2] = -999.9;
  
  event.decayPiPos[0] = -999.9;
  event.decayPiPos[1] = -999.9;
  event.decayPiPos[2] = -999.9;
  
  event.cdistDecayPi = -999.9;
  event.thetaDecayPi = -999.9;
  event.chisqrDecayPiXY = -999.9;
  event.chisqrDecayPiZ = -999.9;

  event.decayPiMomCal = -999.9;
  event.decayPiMomVecCal[0] = -999.9;
  event.decayPiMomVecCal[1] = -999.9;
  event.decayPiMomVecCal[2] = -999.9;

  event.decayNMomCal = -999.9;
  event.decayNMomVecCal[0] = -999.9;
  event.decayNMomVecCal[1] = -999.9;
  event.decayNMomVecCal[2] = -999.9;

  event.vertexNpScat[0] = -999.9;
  event.vertexNpScat[1] = -999.9;
  event.vertexNpScat[2] = -999.9;

  event.cdistNpScat = -999.9;
  event.thetaNpScat = -999.9;
  event.scatNpMomCal = -999.9;
  event.scatNpEkinCal = -999.9;
  event.thetaNpScatCM = -999.9;

  //SksTrack
  event.Hits     = -1;

  event.chisqr= -999.0;
  event.P     = -999.0;
  event.Theta = -999.0;
  event.Phi   = -999.0;
  event.X0    = -999.0;
  event.Y0    = -999.0;
  event.U0    = -999.0;
  event.V0    = -999.0;
  event.path  = -999.0;
  event.Pb    = -999.0;
  event.Ub    = -999.0;
  event.Vb    = -999.0;

  event.dP    = -999.0;
  event.dTheta= -999.0;
  event.dPhi  = -999.0;
  event.dX    = -999.0;
  event.dY    = -999.0;
  event.dU    = -999.0;
  event.dV    = -999.0;

  event.m2        = -999.0;
  event.RecoilMom = -999.0;
  event.beta      = -999.0;
  event.Ekin      = -999.0;
  event.MissMass  = -999.0;
  event.BE        = -999.0;
  event.BE_SMF    = -999.0;

  event.DC1Hits  = -1;
  event.DC2Hits  = -1;
  event.DC3Hits  = -1;
  event.DC4Hits  = -1;
  for( int it=0; it<MaxHits; it++){
    event.DC1Posl[it] = -999.0; 
    event.DC1Resd[it] = -999.0;
    event.DC2Posl[it] = -999.0; 
    event.DC2Resd[it] = -999.0;
    event.DC3Posl[it] = -999.0; 
    event.DC3Resd[it] = -999.0;
    event.DC4Posl[it] = -999.0; 
    event.DC4Resd[it] = -999.0;
  }

  event.tofHits  = -1;
  event.smfHits  = -1;
  for( int it=0; it<MaxHits; it++){
    event.tofTime[it]  = -999.0;
    event.tofEdep[it]  = -999.0;
    event.smfTime[it]  = -999.0;
    event.smfEdep[it]  = -999.0;
  }

  for (int i=0; i<FiberLayerNum+1; i++) 
    event.FiberHits[i] = 0;
  for (int i=0; i<FiberLayerNum+1; i++) {
    for (int j=0; j<MaxHits2; j++) {
      event.FiberSeg[i][j] = -1;
      event.FiberTime[i][j] = -999.0;
      event.FiberEdep[i][j] = -999.0;
      event.FiberPID[i][j] = -1;
    }
  }

  for (int i=0; i<CrystalLayerNum+1; i++) 
    event.CrystalHits[i]= 0;
  for (int i=0; i<CrystalLayerNum+1; i++) {
    for (int j=0; j<MaxHits2; j++) {
      event.CrystalSeg [i][j] = -1;
      event.CrystalTime[i][j] = -999.0;
      event.CrystalEdep[i][j] = -999.0;
      event.CrystalPID [i][j] = -1;
    }
  }


  for (int i=0; i<FiberLayerNum+1; i++) {
    event.fClEdepForScat[i] = -999.0;
    event.fClEdepNormForScat[i] = -999.0;
    event.fClPathLengthForScat[i] = -999.0;

    event.fClEdepForPi[i] = -999.0;
    event.fClEdepNormForPi[i] = -999.0;
    event.fClPathLengthForPi[i] = -999.0;
  }

  event.cClEdepForScat = -999.0;
  event.cClEdepForPi = -999.0;
  event.cRcEdepForScat = 0.0;
  event.cRcEdepForPi = 0.0;
}


const int NBin1HPDcIn[13] = {
  100, 
  200, 200,
  200, 200, 200, 200, 200, 200, 
  300, 300, 300, 300
};

const int NBin2HPDcIn[13] = {
  100, 
  50, 50, 
  50, 50, 50, 50, 
  100, 100, 100, 100, 100, 100
};

const double MinHPDcIn[13] = {
  0.0, 
  // SDC0
  -200.0, -200.0,
  // SDC1
  -200.0, -200.0, -200.0, -200.0, -200.0, -200.0,
  // SDC2
  -300.0, -300.0, -300.0, -300.0
};

const double MaxHPDcIn[13] = {
  0.0,
  // SDC0
  200.0, 200.0,
  // SDC1
  200.0, 200.0, 200.0, 200.0, 200.0, 200.0,
  // SDC2
  300.0, 300.0, 300.0, 300.0
};

const int NBin1HPDcOut[13] = {
  100, 
  500, 500, 500, 500, 500, 500,
  500, 500, 500, 500, 500, 500 
};

const int NBin2HPDcOut[13] = {
  100, 
  100, 100, 100, 100, 100, 100,
  100, 100, 100, 100, 100, 100
};

const double MinHPDcOut[13] = {
  0.0,
  // SDC3
  -1200.0, -1200.0, -1200.0, -1200.0, -1200.0, -1200.0,
  // SDC4X
  -1200.0, -1200.0, -1200.0, -1200.0, -1200.0, -1200.0
};

const double MaxHPDcOut[13] = {
  0.0,
  // SDC3
  1200.0, 1200.0, 1200.0, 1200.0, 1200.0, 1200.0,
  // SDC4X 
  1200.0, 1200.0, 1200.0, 1200.0, 1200.0, 1200.0
};

void DefineHistograms( const char *filename )
{
  new TFile( filename, "recreate" );

  //Tree
  HBTree("tree","tree of Sks");
  TTree *tree = dynamic_cast<TTree *>(gFile->Get("tree"));


  char buf1[100], buf2[100];
  sprintf(buf1, "fClEdepForScat[%d]/D", FiberLayerNum+1);
  tree->Branch("fClEdepForScat", event.fClEdepForScat, buf1);
  sprintf(buf1, "fClEdepNormForScat[%d]/D", FiberLayerNum+1);
  tree->Branch("fClEdepNormForScat", event.fClEdepNormForScat, buf1);
  sprintf(buf1, "fClPathLengthForScat[%d]/D", FiberLayerNum+1);
  tree->Branch("fClPathLengthForScat", event.fClPathLengthForScat, buf1);

  sprintf(buf1, "fClEdepForPi[%d]/D", FiberLayerNum+1);
  tree->Branch("fClEdepForPi", event.fClEdepForPi, buf1);
  sprintf(buf1, "fClEdepNormForPi[%d]/D", FiberLayerNum+1);
  tree->Branch("fClEdepNormForPi", event.fClEdepNormForPi, buf1);
  sprintf(buf1, "fClPathLengthForPi[%d]/D", FiberLayerNum+1);
  tree->Branch("fClPathLengthForPi", event.fClPathLengthForPi, buf1);

  tree->Branch("cClEdepForScat", &event.cClEdepForScat, "cClEdepForScat/D");
  tree->Branch("cClEdepForPi", &event.cClEdepForPi, "cClEdepForPi/D");
  tree->Branch("cRcEdepForScat", &event.cRcEdepForScat, "cRcEdepForScat/D");
  tree->Branch("cRcEdepForPi", &event.cRcEdepForPi, "cRcEdepForPi/D");

}
