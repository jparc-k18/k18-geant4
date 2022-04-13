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
#include "EnergyCorrection.hh" 
#include "EvDisp.hh" 

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

  double thetaDecayPi2;
  double decayPiMomCal2;	     

  double vertexNpScat[3];
  double cdistNpScat;
  double thetaNpScat;
  double scatNpMomCal;
  double scatNpEkinCal;
  double thetaNpScatCM;

  double vertexPiPScat[3];
  double cdistPiPScat;
  double thetaDecayPi3;
  double decayPiMomCal3;
  double decayPiMom3;

  double cdistLambdaDecay;
  double vertexLambdaDecay[3];
  double thetaLambdaDecay;
  double momPiFromLambda;
  double momProtonFromLambda;
  double momLambda;
  double momVecPiFromLambda[3];
  double momVecProtonFromLambda[3];
  double momVecLambda[3];
  double cdistLambdaNConv;
  double thetaLambdaNConv;
  double vertexLambdaNConv[3];
  double momCalLambda;
  double momCalLambda2;
  double thetaCMLambdaNConv;
  double invmass_pip;

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

};
static Event event;

const int MinHitForPi = 6;
const int MinHitForP = 6;

static EvDisp evDisp;

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
  if( argc<4 ) {
    std::cerr << "Usage: " << argv[0] 
	      << " confFile fieldMap dataFile"
	      << std::endl;
    exit(-1);
  }

  DCConfMan::GetInstance().Initialize(argv[1]);
  FieldMan::GetInstance().Initialize(argv[2]);

  signal(SIGINT,closeFile);

  evDisp.Initialize();

  inFile = argv[3];
  std::ifstream InputData( inFile );
  while( ProcessOneEvent( InputData ) ){
#ifdef Interactive 
    int it; 
    std::cout << "# ";
    std::cin >> it;
    if( it<0 ) break;
#endif
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

  const int maxMultiSdcIn = 3;
  const int maxMultiSdcOut = 3;
  const int maxMultiBcOut = 3;

  static int events=0;
  events++;
  //if (events%100==0)
    std::cout<<"Events:"<< events<<std::endl;

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

  int multiSdcIn[11];
  for (int i=0; i<11; i++) 
    multiSdcIn[i] = 0;
  int multiSdcOut[13];
  for (int i=0; i<13; i++) 
    multiSdcOut[i] = 0;
  int multiBcOut[13];
  for (int i=0; i<13; i++) 
    multiBcOut[i] = 0;

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
      if (lnumold>=1 && lnumold<=10)
	multiSdcIn[lnumold]++;
      if (lnumold>=41 && lnumold<=52)
	multiSdcOut[lnumold-40]++;

    }
    if( lnumold>112 && lnumold <150){
      DCAna.SetDCHits(lnumold, pos);
      if (lnumold>=113 && lnumold<=124)
	multiBcOut[lnumold-112]++;
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
	if (flayer == 1 || flayer == 3 || flayer == 5 || flayer == 7 || flayer == 9) {
	  evDisp.DrawHitFiber(flayer, segment, pid, edep);
	  evDisp.DrawHitRFiber(flayer, segment, pid, edep);
	} else {
	  evDisp.DrawHitUVFiber(flayer, segment, pid, edep);
	}
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
	/*
	  std::cout << "Crystal "
	  << ", Layer = " << clayer
	  << ", Segment = " << segment
	  << ", Edep = " << edep
	  << ", PID = " << pid << std::endl;
	*/
	evDisp.DrawHitCalori(clayer, segment, pid, edep);

	if (pid == 10 || pid == 11 || pid == 14)
	  CylAna.SetCrystalDataForScatPart(segment, edep);
	else if (pid == 1 || pid == 2)
	  CylAnaPi.SetCrystalDataForScatPart(segment, edep);


      } else {
	std::cerr << "Too many crystal hits : " << event.CrystalHits[clayer]
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


  /*
  if (event.scatFlag != 2)
    return true;
  */
  std::cout << "###" << std::endl;


  for (int i=1; i<=10; i++) {
    //std::cout << "MultiIn[" << i << "] = " << multiSdcIn[i] << std::endl;
    if (multiSdcIn[i] > maxMultiSdcIn) {
      std::cerr << "Too many multi: SdcIn[" << i << "], " << multiSdcIn[i]
		<< std::endl;
      evDisp.EndOfEvent();
      return true;
    }
  }
  for (int i=1; i<=12; i++) {
    //std::cout << "MultiOut[" << i << "] = " << multiSdcOut[i] << std::endl;
    if (multiSdcOut[i] > maxMultiSdcOut) {
      std::cerr << "Too many multi: SdcOut[" << i << "], " << multiSdcOut[i]
		<< std::endl;
      evDisp.EndOfEvent();
      return true;
    }

  }
  for (int i=1; i<=12; i++) {
    //std::cout << "MultiBcOut[" << i << "] = " << multiBcOut[i] << std::endl;
    if (multiBcOut[i] > maxMultiBcOut) {
      std::cerr << "Too many multi: BcOut[" << i << "], " << multiBcOut[i]
		<< std::endl;
      evDisp.EndOfEvent();
      return true;
    }
  }

  DCAna.TrackSearchSdcIn();
  int ntSdcIn = DCAna.GetNtracksSdcIn();
  std::cout << "ntSdcIn = " << ntSdcIn << std::endl;
  std::vector <ThreeVector> sMom1Cont, sPos1Cont;
  for (int i=0; i<ntSdcIn; i++) {
    DCLocalTrack *trIn = DCAna.GetTrackSdcIn(i);
    bool stIn=trIn->GetStatus();
    int nhin=trIn->GetNHit();

    double chisqrIn=trIn->GetChiSquare();
    double x0In=trIn->GetX0(), y0In=trIn->GetY0();
    double u0In=trIn->GetU0(), v0In=trIn->GetV0();
    double thetaIn=acos(1./sqrt(1.+u0In*u0In+v0In*v0In))*Rad2Deg;
    double phiIn=atan2(v0In,u0In)*Rad2Deg;
    ThreeVector sMom1(u0In/sqrt(1.+u0In*u0In+v0In*v0In),
		      v0In/sqrt(1.+u0In*u0In+v0In*v0In),
		      1./sqrt(1.+u0In*u0In+v0In*v0In));
    ThreeVector sPos1(x0In, y0In, 0);
    sMom1Cont.push_back(sMom1);
    sPos1Cont.push_back(sPos1);

    if( stIn ){
      event.u0SdcIn = u0In;
      event.v0SdcIn = v0In;
    }
  }

  DCAna.TrackSearchSdcOut();
  int ntSdcOut = DCAna.GetNtracksSdcOut();

  DCAna.TrackSearchBcOut();
  int ntBcOut = DCAna.GetNtracksBcOut();
  std::cout << "ntBcOut = " << ntBcOut << std::endl;
  std::vector <ThreeVector> bMom1Cont, bPos1Cont;
  for (int i=0; i<ntBcOut; i++) {
    DCLocalTrack *trBcOut = DCAna.GetTrackBcOut(i);
    bool stBcOut=trBcOut->GetStatus();
    int nhBcOut=trBcOut->GetNHit();
    double chisqrBcOut=trBcOut->GetChiSquare();
    double x0BcOut=trBcOut->GetX0(), y0BcOut=trBcOut->GetY0();
    double u0BcOut=trBcOut->GetU0(), v0BcOut=trBcOut->GetV0();
    double thetaBcOut=acos(1./sqrt(1.+u0BcOut*u0BcOut+v0BcOut*v0BcOut))*Rad2Deg;
    double phiBcOut=atan2(v0BcOut,u0BcOut)*Rad2Deg;
    ThreeVector bMom1(u0BcOut/sqrt(1.+u0BcOut*u0BcOut+v0BcOut*v0BcOut),
		      v0BcOut/sqrt(1.+u0BcOut*u0BcOut+v0BcOut*v0BcOut),
		      1./sqrt(1.+u0BcOut*u0BcOut+v0BcOut*v0BcOut));
    double beammomentum=1.3;
    bMom1 *= beammomentum;
    ThreeVector bPos1(x0BcOut, y0BcOut, 0);

    bMom1Cont.push_back(bMom1);
    bPos1Cont.push_back(bPos1);
    
    if( stBcOut ){
      event.u0BcOut = u0BcOut;
      event.v0BcOut = v0BcOut;
    }
  }

  /*
  if( ntSdcIn == 0 || ntSdcOut == 0 || ntBcOut == 0){
    return true; 
  }
  */

  DCAna.TrackSearchSks();
  int ntSks = DCAna.GetNTracksSks();
  std::cout << "ntSks = " << ntSks << std::endl;
  std::vector <ThreeVector> priPosCont, priMomCont;
  std::vector <ThreeVector> priPosSksInCont, priMomSksInCont;
  for (int i=0; i<ntSks; i++) {
    SksTrack *track = DCAna.GetSksTrack(i);
    bool status=track->Status();
    double chisqr=track->chisqr();
    
    ThreeVector priPos = track->PrimaryPosition();
    ThreeVector priMom = track->PrimaryMomentum();
    priPosCont.push_back(priPos);
    priMomCont.push_back(priMom);
    if (i==0) {
      event.priPosVec[0] = priPos.x();
      event.priPosVec[1] = priPos.y();
      event.priPosVec[2] = priPos.z();
      event.priMomVec[0] = priMom.x();
      event.priMomVec[1] = priMom.y();
      event.priMomVec[2] = priMom.z();
    }

    double rx=priPos.x(), ry=priPos.y(), rz=priPos.z();
    double rp=priMom.mag();
    double rt=priMom.theta()*Rad2Deg, rf=priMom.phi()*Rad2Deg;
    double lengTof = track->PathLengthToTOF();
    
    int nhdc = track->GetNHits();
    
    if (status) {
      bool statusSksIn = track->doFitSksIn();
      if (statusSksIn) {
	double chisqrSksIn=track->chisqrSksIn();
	ThreeVector priPosSksIn = track->PrimaryPositionSksIn();
	ThreeVector priMomSksIn = track->PrimaryMomentumSksIn();
	priPosSksInCont.push_back(priPosSksIn);
	priMomSksInCont.push_back(priMomSksIn);
	
	if (i==0) {
	  event.priPosSksInVec[0] = priPosSksIn.x();
	  event.priPosSksInVec[1] = priPosSksIn.y();
	  event.priPosSksInVec[2] = priPosSksIn.z();
	  event.priMomSksInVec[0] = priMomSksIn.x();
	  event.priMomSksInVec[1] = priMomSksIn.y();
	  event.priMomSksInVec[2] = priMomSksIn.z();
	}
      }
    }
    
    
    double gu = tan(gt*Deg2Rad)*cos(gf*Deg2Rad);
    double gv = tan(gt*Deg2Rad)*sin(gf*Deg2Rad);
    double ru = tan(rt*Deg2Rad)*cos(rf*Deg2Rad);
    double rv = tan(rt*Deg2Rad)*sin(rf*Deg2Rad);
    
    event.Hits  = nhdc; 
    event.chisqr= chisqr; 
    
    event.P     = rp;
    event.Theta = rt;//*Deg2Rad;
    event.Phi   = rf;//*Deg2Rad;
    event.X0    = rx;
    event.Y0    = ry;
    event.U0    = ru;
    event.V0    = rv;
    event.path  = lengTof;
    event.Pb    = gpb;
    event.Ub    = gub;
    event.Vb    = gvb;
    
    event.dP    = rp-gp;
    event.dTheta= rt-gt;
    event.dPhi  = rf-gf;
    event.dX    = rx-gx;
    event.dY    = ry-gy;
    event.dU    = ru-gu;
    event.dV    = rv-gv;
    
    //Mass Square
    double beta = lengTof/(event.tofTime[0]*299.792);
    double M2   = ((rp*rp)/(beta*beta))*(1-beta*beta);
    event.m2    = M2;
  }

  bool tofflag = false;
  bool smfflag = false;
  
  if( tofedep>4 && 12<toftime && toftime<24 ) tofflag = true;
  if( smfedep>6 && 12<smftime && smftime<30 ) smfflag = true;
  
  // Vertex2
  if( ntSdcIn>=1 && ntBcOut>=1 && ntSks>=1 ){
    ThreeVector bMom1 = bMom1Cont[0];
    ThreeVector bPos1 = bPos1Cont[0];
    ThreeVector priMom = priMomCont[0];
    ThreeVector priPos = priPosCont[0];

    if (priPosSksInCont.size()>0 && priMomSksInCont.size()>0 ) {
      ThreeVector priPosSksIn = priPosSksInCont[0];
      ThreeVector priMomSksIn = priMomSksInCont[0];
      
      double cost5=bMom1*priMomSksIn/(bMom1.mag()*priMomSksIn.mag());
      ThreeVector Vert5=VertexPoint( bPos1, priPosSksIn, bMom1, priMomSksIn );
      double cdist5=closeDist( bPos1, priPosSksIn, bMom1, priMomSksIn );
      double theta5=acos(cost5)*Rad2Deg;
      
      event.vertex5[0] = Vert5.x();
      event.vertex5[1] = Vert5.y();
      event.vertex5[2] = Vert5.z();
      event.theta5 = theta5;
      event.cdist5 = cdist5;
    }
  }

  if (ntSks<=0) {
    evDisp.EndOfEvent();
    return true;
  }

  //Missing Mass
  ////Beam Mom

  double Pbz = gpb/sqrt(1.0+gub*gub+gvb*gvb);
  double Pbx = Pbz*gub;
  double Pby = Pbz*gvb;
  ThreeVector bMom( Pbx, Pby ,Pbz );

  SksTrack *track = DCAna.GetSksTrack(0);
  bool status=track->Status();
  double chisqr=track->chisqr();
  int nhdc = track->GetNHits();
  ThreeVector priPos = track->PrimaryPosition();
  ThreeVector priMom = track->PrimaryMomentum();
  
  double rx=priPos.x(), ry=priPos.y(), rz=priPos.z();
  double rp=priMom.mag();
  double rt=priMom.theta()*Rad2Deg, rf=priMom.phi()*Rad2Deg;
  double lengTof = track->PathLengthToTOF();
  double gu = tan(gt*Deg2Rad)*cos(gf*Deg2Rad);
  double gv = tan(gt*Deg2Rad)*sin(gf*Deg2Rad);
  double ru = tan(rt*Deg2Rad)*cos(rf*Deg2Rad);
  double rv = tan(rt*Deg2Rad)*sin(rf*Deg2Rad);
  ////Scattered Mom
  double Pz = rp/sqrt(1.0+ru*ru+rv*rv);
  double Px = Pz*ru;
  double Py = Pz*rv;
  ThreeVector sMom( Px, Py ,Pz );

  ThreeVector bMom1 = bMom1Cont[0];
  ThreeVector bPos1 = bPos1Cont[0];

  LorentzVector LvBeam(  bMom1, sqrt(PionMass*PionMass+bMom1.mag2()) );
  LorentzVector LvScat( priMom, sqrt(KaonMass*KaonMass+priMom.mag2()) );
  LorentzVector LvTgt( 0., 0., 0., ProtonMass );
  LorentzVector LvRc = LvBeam+LvTgt-LvScat;
  double mismass = LvRc.mag();//-LvCore.mag();

  event.MissMass = mismass;

  ThreeVector HypBeamVec = LvRc.vect();
  event.hypBeamVec[0] = HypBeamVec.x();
  event.hypBeamVec[1] = HypBeamVec.y();
  event.hypBeamVec[2] = HypBeamVec.z();

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

    flagFitPiXY = CylAnaPi.FitXYPlane();
    if (flagFitPiXY)
      flagFitPiZ = CylAnaPi.FitZPlane();

    double minPos = 999.;
    double maxPos = -999.;
    for (int i=2; i<=FiberLayerNum; i += 2) {
      double pos[3];
      CylAnaPi.GetHitPosition(i, pos);
      evDisp.DrawHitXZPos(pos, 2);
      evDisp.DrawHitYZPos(pos, 2);

      if (minPos > pos[2])
	minPos = pos[2];

      if (maxPos < pos[2])
	maxPos = pos[2];
    }

    if (flagFitPiXY && flagFitPiZ) {

      //std::cout << "CylAnaPi" << std::endl;
      ThreeVector Vert5(event.vertex5[0], event.vertex5[1], event.vertex5[2]);
      ThreeVector DecayPiPos = CylAnaPi.GetScatPosition();
      ThreeVector DecayPiVec = CylAnaPi.GetScatVector();

      double pos1[3], pos2[3];
      pos1[2] = minPos-100.;
      pos1[0] = DecayPiVec.x()/DecayPiVec.z()*(pos1[2]-DecayPiPos.z())+DecayPiPos.x();
      pos1[1] = DecayPiVec.y()/DecayPiVec.z()*(pos1[2]-DecayPiPos.z())+DecayPiPos.y();

      pos2[2] = maxPos+100.;
      pos2[0] = DecayPiVec.x()/DecayPiVec.z()*(pos2[2]-DecayPiPos.z())+DecayPiPos.x();
      pos2[1] = DecayPiVec.y()/DecayPiVec.z()*(pos2[2]-DecayPiPos.z())+DecayPiPos.y();
      evDisp.DrawTrack(pos1, pos2, 2);

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

    flagFitXY = CylAna.FitXYPlane();
    if (flagFitXY)
      flagFitZ = CylAna.FitZPlane();

    double minPos = 999.;
    double maxPos = -999.;

    for (int i=2; i<=FiberLayerNum; i += 2) {
      double pos[3];
      CylAna.GetHitPosition(i, pos);
      evDisp.DrawHitXZPos(pos, 10);
      evDisp.DrawHitYZPos(pos, 10);

      if (minPos > pos[2])
	minPos = pos[2];

      if (maxPos < pos[2])
	maxPos = pos[2];
    }

    if (flagFitXY && flagFitZ) {
      //std::cout << "CylAnaP" << std::endl;
      ThreeVector Vert5(event.vertex5[0], event.vertex5[1], event.vertex5[2]);
      //ThreeVector HypBeamVec = LvRc.vect();
      ThreeVector ScatPos = CylAna.GetScatPosition();
      ThreeVector ScatVec = CylAna.GetScatVector();
      double cdist;
      ThreeVector VertScat = VertexPoint3D( Vert5, ScatPos, HypBeamVec, ScatVec, cdist);

      double pos1[3], pos2[3];
      pos1[2] = minPos-100.;
      pos1[0] = ScatVec.x()/ScatVec.z()*(pos1[2]-ScatPos.z())+ScatPos.x();
      pos1[1] = ScatVec.y()/ScatVec.z()*(pos1[2]-ScatPos.z())+ScatPos.y();

      pos2[2] = maxPos+100.;
      pos2[0] = ScatVec.x()/ScatVec.z()*(pos2[2]-ScatPos.z())+ScatPos.x();
      pos2[1] = ScatVec.y()/ScatVec.z()*(pos2[2]-ScatPos.z())+ScatPos.y();
      evDisp.DrawTrack(pos1, pos2, 10);

      for (int i=1; i<=FiberLayerNum; i++) {
	double edepNorm = CylAna.GetFiberClusterEdepNorm(i);
	event.fClEdepNormForScat[i] = edepNorm;
	double pathLength = CylAna.GetFiberPathLength(i);
	event.fClPathLengthForScat[i] = pathLength;
      }
    }
  }

  evDisp.EndOfEvent();
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

  event.thetaDecayPi2 = -999.9;
  event.decayPiMomCal2 = -999.9;

  event.vertexPiPScat[0] = -999.9;
  event.vertexPiPScat[1] = -999.9;
  event.vertexPiPScat[2] = -999.9;
  event.cdistPiPScat     = -999.9;;
  event.thetaDecayPi3    = -999.9;;
  event.decayPiMomCal3   = -999.9;;
  event.decayPiMom3      = -999.9;;

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


  event.cdistLambdaDecay = -999.9;
  event.thetaLambdaDecay = -999.9;
  event.vertexLambdaDecay[0] = -999.9;
  event.vertexLambdaDecay[1] = -999.9;
  event.vertexLambdaDecay[2] = -999.9;

  event.momPiFromLambda = -999.9;
  event.momProtonFromLambda = -999.9;
  event.momLambda = -999.9;
  event.momVecPiFromLambda[0] = -999.9;
  event.momVecPiFromLambda[1] = -999.9;
  event.momVecPiFromLambda[2] = -999.9;
  event.momVecProtonFromLambda[0] = -999.9;
  event.momVecProtonFromLambda[1] = -999.9;
  event.momVecProtonFromLambda[2] = -999.9;
  event.momVecLambda[0] = -999.9;
  event.momVecLambda[1] = -999.9;
  event.momVecLambda[2] = -999.9;

  event.cdistLambdaNConv = -999.9;
  event.thetaLambdaNConv = -999.9;
  event.vertexLambdaNConv[0] = -999.9;
  event.vertexLambdaNConv[1] = -999.9;
  event.vertexLambdaNConv[2] = -999.9;

  event.momCalLambda = -999.9;
  event.momCalLambda2 = -999.9;
  event.thetaCMLambdaNConv = -999.9;
  event.invmass_pip = -999.9;

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

