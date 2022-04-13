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

  const int maxMultiSdcIn = 3;
  const int maxMultiSdcOut = 3;
  const int maxMultiBcOut = 3;

  static int events=0;
  events++;
  if (events%100==0)
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
  //std::cout << "###" << std::endl;

  HF1( 101, 1. );


  for (int i=1; i<=10; i++) {
    //std::cout << "MultiIn[" << i << "] = " << multiSdcIn[i] << std::endl;
    if (multiSdcIn[i] > maxMultiSdcIn) {
      std::cerr << "Too many multi: SdcIn[" << i << "], " << multiSdcIn[i]
		<< std::endl;
      return true;
    }
  }
  for (int i=1; i<=12; i++) {
    //std::cout << "MultiOut[" << i << "] = " << multiSdcOut[i] << std::endl;
    if (multiSdcOut[i] > maxMultiSdcOut) {
      std::cerr << "Too many multi: SdcOut[" << i << "], " << multiSdcOut[i]
		<< std::endl;
      return true;
    }

  }
  for (int i=1; i<=12; i++) {
    //std::cout << "MultiBcOut[" << i << "] = " << multiBcOut[i] << std::endl;
    if (multiBcOut[i] > maxMultiBcOut) {
      std::cerr << "Too many multi: BcOut[" << i << "], " << multiBcOut[i]
		<< std::endl;
      return true;
    }
  }

  DCAna.TrackSearchSdcIn();
  int ntSdcIn = DCAna.GetNtracksSdcIn();
  //std::cout << "ntSdcIn = " << ntSdcIn << std::endl;
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

    HF1( 3103, double(nhin) );
    if( stIn ){
      event.u0SdcIn = u0In;
      event.v0SdcIn = v0In;
      HF1( 3102, chisqrIn );
      HF1( 3104, x0In );
      HF1( 3105, y0In );
      HF1( 3106, u0In );
      HF1( 3107, v0In );
      HF1( 3108, thetaIn );
      HF1( 3109, phiIn );
      HF2( 3111, x0In, u0In );
      HF2( 3112, y0In, v0In );
      HF2( 3113, x0In, y0In );
      HF2( 3114, thetaIn, phiIn );
      
      double gu = tan(gt*Deg2Rad)*cos(gf*Deg2Rad);
      double gv = tan(gt*Deg2Rad)*sin(gf*Deg2Rad);
      
      HF1( 3115, x0In-gx );
      HF1( 3116, y0In-gy );
      HF1( 3117, u0In-gu );
      HF1( 3118, v0In-gv );
      
      for( std::size_t i=0; i<nhin; ++i ){
	DCHit *hit = trIn->GetHit(i);
	int layer = hit->GetLayer();
	HF1( 3110, layer );
	HF1( 3120+layer, hit->GetLocalHitPos() );
	HF1( 3140+layer, hit->GetResidual() );
	//std::cout << "Residual=" << hit->GetResidual() << std::endl;
	HF1( 3160+layer, hit->GetLocalHitPos(), hit->GetResidual() );
      }
    }
  }

  DCAna.TrackSearchSdcOut();
  int ntSdcOut = DCAna.GetNtracksSdcOut();
  //std::cout << "ntSdcOut = " << ntSdcOut << std::endl;
  for (int i=0; i<ntSdcOut; i++) {
    DCLocalTrack *trOut = DCAna.GetTrackSdcOut(i);
    bool stOut=trOut->GetStatus();
    int nhout=trOut->GetNHit();
    double chisqrOut=trOut->GetChiSquare();
    double z0Out=DCConfMan::GetInstance().GetLocalZ(61);
    double x0Out=trOut->GetX( z0Out ), y0Out=trOut->GetY( z0Out );
    double u0Out=trOut->GetU0(), v0Out=trOut->GetV0();
    double thetaOut=acos(1./sqrt(1.+u0Out*u0Out+v0Out*v0Out))*Rad2Deg;
    double phiOut=atan2(v0Out,u0Out)*Rad2Deg;
    
    HF1( 4101, stOut );
    HF1( 4103, nhout );
    if( stOut ){
      HF1( 4102, chisqrOut );
      HF1( 4104, x0Out );
      HF1( 4105, y0Out );
      HF1( 4106, u0Out );
      HF1( 4107, v0Out );
      HF1( 4108, thetaOut );
      HF1( 4109, phiOut );
      HF2( 4111, x0Out, u0Out );
      HF2( 4112, y0Out, v0Out );
      HF2( 4113, x0Out, y0Out );
      HF2( 4114, thetaOut, phiOut );
      for( std::size_t i=0; i<nhout; ++i ){
	DCHit *hit = trOut->GetHit(i);
	int layer = hit->GetLayer();
	HF1( 4110, layer );
	HF1( 4090+layer, hit->GetLocalHitPos() );
	HF1( 4110+layer, hit->GetResidual() );
	//std::cout << "Layer=" << layer << std::endl;
	//std::cout << "Residual=" << hit->GetResidual() << std::endl;
	HF1( 4130+layer, hit->GetLocalHitPos(), hit->GetResidual() );
      }
    }
  }

  DCAna.TrackSearchBcOut();
  int ntBcOut = DCAna.GetNtracksBcOut();
  //std::cout << "ntBcOut = " << ntBcOut << std::endl;
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
    
    HF1( 5103, double(nhBcOut) );
    if( stBcOut ){
      event.u0BcOut = u0BcOut;
      event.v0BcOut = v0BcOut;
      HF1( 5102, chisqrBcOut );
      HF1( 5104, x0BcOut );
      HF1( 5105, y0BcOut );
      HF1( 5106, u0BcOut );
      HF1( 5107, v0BcOut );
      HF1( 5108, thetaBcOut );
      HF1( 5109, phiBcOut );
      HF2( 5111, x0BcOut, u0BcOut );
      HF2( 5112, y0BcOut, v0BcOut );
      HF2( 5113, x0BcOut, y0BcOut );
      HF2( 5114, thetaBcOut, phiBcOut );
      
      for( std::size_t i=0; i<nhBcOut; ++i ){
	DCHit *hit = trBcOut->GetHit(i);
	int layer = hit->GetLayer()-112;
	HF1( 5110, layer );
	HF1( 5120+layer, hit->GetLocalHitPos() );
	HF1( 5140+layer, hit->GetResidual() );
	//std::cout << "Residual=" << hit->GetResidual() << std::endl;
	HF1( 5160+layer, hit->GetLocalHitPos(), hit->GetResidual() );
      }
    }
    /*
      std::cout << "stIn = " << stIn << ", stOut = " << stOut
      << "stBcOut = " << stBcOut << std::endl;
    */
  }

  /*
  if( ntSdcIn == 0 || ntSdcOut == 0 || ntBcOut == 0){
    return true; 
  }
  */

  HF1( 101, 2. );

  DCAna.TrackSearchSks();
  int ntSks = DCAna.GetNTracksSks();
  //std::cout << "ntSks = " << ntSks << std::endl;
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
	
	HF1(545, chisqrSksIn);
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
    
    //#ifdef Interactive
#if 0
    {
      std::ios::fmtflags oldFlags = std::cout.flags();
      std::size_t oldPrec = std::cout.precision();
      std::cout.setf( std::ios::fixed );
      std::cout.precision(4);
      
      std::cout << "[SksTracking]: St=" << status
		<< " chi2=" << std::setw(10) << chisqr 
		<< " PL=" << std::setprecision(1) 
		<< std::setw(7) << lengTof  << " [mm]" << std::endl;
      std::cout.precision(2);
      std::cout << " Position: ( "
		<< std::setw(7) << rx << ", "
		<< std::setw(7) << ry << " )" << " - ( "
		<< std::setw(7) << gx << ", "
		<< std::setw(7) << gy << " )" << " =====> ( "
		<< std::setw(7) << rx-gx << ", "
		<< std::setw(7) << ry-gy << " ) [mm]" << std::endl;
      std::cout.precision(5);
      std::cout << " P: " << std::setw(8) << rp << " ( "
		<< std::setw(8) << gp << " )[GeV/c] < "
		<< std::setw(8) << (rp-gp)*1000 << " [MeV/c]> ";
      std::cout.precision(1);
      std::cout << "  Angles: ( " 
		<< std::setw(5) << rt << ", " 
		<< std::setw(6) << rf << " ) <( "
		<< std::setw(5) << gt << ", " 
		<< std::setw(5) << gf << " )> [degree]"
		<< std::endl; 
      
      std::cout.flags( oldFlags );
      std::cout.precision( oldPrec );
    }
#endif
    
    double gu = tan(gt*Deg2Rad)*cos(gf*Deg2Rad);
    double gv = tan(gt*Deg2Rad)*sin(gf*Deg2Rad);
    double ru = tan(rt*Deg2Rad)*cos(rf*Deg2Rad);
    double rv = tan(rt*Deg2Rad)*sin(rf*Deg2Rad);
    
    HF1( 101, 3. );
    
    HF1( 102, status );
    HF1( 103, chisqr );
    HF1( 104, nhdc ); 
    HF1( 105, rx ); HF1( 106, ry ); HF1( 107, rp );
    HF1( 108, rt*Deg2Rad ); HF1( 109, rf );
    HF1( 110, lengTof );
    HF1( 111, rx-gx ); HF1( 112, ry-gy ); HF1( 113, rp-gp );
    HF1( 114, rt-gt ); HF1( 115, rf-gf );
    HF1( 116, ru-gu ); HF1( 117, rv-gv );
    HF2( 118, gp, rp-gp );
    HF1( 119, 1000*(rp-gp)/gp );
    
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
  

  // Vertex1
  if( ntSdcIn >= 1 && ntBcOut >= 1 ){
    ThreeVector sMom1 = sMom1Cont[0];
    ThreeVector sPos1 = sPos1Cont[0];
    ThreeVector bMom1 = bMom1Cont[0];
    ThreeVector bPos1 = bPos1Cont[0];

    double dx=sPos1.x()-bPos1.x();
    double dy=sPos1.y()-bPos1.y();
    double us=sMom1.x()/sMom1.z(), vs=sMom1.y()/sMom1.z();
    double ub=bMom1.x()/bMom1.z(), vb=bMom1.y()/bMom1.z();
    
    double cost=bMom1*sMom1/(bMom1.mag()*sMom1.mag());
    ThreeVector Vert=VertexPointByHonly( bPos1, sPos1, bMom1, sMom1 );
    double cdist=closeDist( bPos1, sPos1, bMom1, sMom1 );
    double theta=acos(cost)*Rad2Deg;

    event.vertex1[0] = Vert.x();
    event.vertex1[1] = Vert.y();
    event.vertex1[2] = Vert.z();
    event.theta1 = theta;
    event.cdist1 = cdist;

    HF1(500, Vert.x());
    HF1(501, Vert.y());
    HF1(502, Vert.z());
    HF1(503, theta);
    HF2(504, Vert.z(), theta);
  }

  /*
  std::cout << "bMom1Cont " << bMom1Cont.size()
	    << ", bPos1Cont " << bPos1Cont.size()
	    << ", priMomCont " << priMomCont.size()
	    << ", priPosCont " << priPosCont.size()
	    << ", priMomSksInCont "<< priMomSksInCont.size()
	    << ", priPosSksInCont "<< priPosSksInCont.size()
	    << std::endl;
  */
  // Vertex2
  if( ntSdcIn>=1 && ntBcOut>=1 && ntSks>=1 ){
    ThreeVector bMom1 = bMom1Cont[0];
    ThreeVector bPos1 = bPos1Cont[0];
    ThreeVector priMom = priMomCont[0];
    ThreeVector priPos = priPosCont[0];

    double dx=priPos.x()-bPos1.x();
    double dy=priPos.y()-bPos1.y();
    double us=priMom.x()/priMom.z(), vs=priMom.y()/priMom.z();
    double ub=bMom1.x()/bMom1.z(), vb=bMom1.y()/bMom1.z();
    
    double cost=bMom1*priMom/(bMom1.mag()*priMom.mag());
    ThreeVector Vert=VertexPoint( bPos1, priPos, bMom1, priMom );
    double cdist=closeDist( bPos1, priPos, bMom1, priMom );
    double theta=acos(cost)*Rad2Deg;

    event.vertex2[0] = Vert.x();
    event.vertex2[1] = Vert.y();
    event.vertex2[2] = Vert.z();
    event.theta2 = theta;
    event.cdist2 = cdist;

    HF1(510, Vert.x());
    HF1(511, Vert.y());
    HF1(512, Vert.z());
    HF1(513, theta);
    HF2(514, Vert.z(), theta);

    ThreeVector VertH=VertexPointByHonly( bPos1, priPos, bMom1, priMom );
    event.vertex3[0] = VertH.x();
    event.vertex3[1] = VertH.y();
    event.vertex3[2] = VertH.z();

    HF1(520, VertH.x());
    HF1(521, VertH.y());
    HF1(522, VertH.z());

    ThreeVector VertexSks;
    double      cdistSks;
    ThreeVector sMomSks;

    DCAna.GetSksTrack(0)->VertexSearch(bMom1, bPos1, VertexSks, cdistSks, sMomSks);
    double cost4=bMom1*sMomSks/(bMom1.mag()*sMomSks.mag());
    double theta4=acos(cost4)*Rad2Deg;
    event.vertex4[0] = VertexSks.x();
    event.vertex4[1] = VertexSks.y();
    event.vertex4[2] = VertexSks.z();
    event.cdist4 = cdistSks;
    event.theta4 = theta4;

    HF1(530, VertexSks.x());
    HF1(531, VertexSks.y());
    HF1(532, VertexSks.z());
    HF1(533, theta4);
    HF2(534, VertexSks.z(), theta4);

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
      
      HF1(540, Vert5.x());
      HF1(541, Vert5.y());
      HF1(542, Vert5.z());
      HF1(543, theta5);
      HF2(544, Vert5.z(), theta5);
    }

      /*
	std::cout << "Vertex (x, y, z) = ( " << VertexSks.x()
	<< ", " << VertexSks.y() << ", " << VertexSks.z() 
	<< "), cdist = " << cdistSks << std::endl;
      */
    //getchar();
  }

  if (ntSks<=0)
    return true;


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

  double pSks = track->PrimaryMomentum().mag();
  double pSks_Cor, eSks_Cor;

  ThreeVector Vert5_Sigma;
  if (event.theta5<5. && fabs(event.vertex5[2])>CyLH2TgtZ/2.) 
    Vert5_Sigma = ThreeVector(event.vertex5[0], event.vertex5[1], 0);
  else
    Vert5_Sigma = ThreeVector(event.vertex5[0], event.vertex5[1], event.vertex5[2]);

  CorrElossOut(&pSks_Cor, &eSks_Cor, pSks, KAON, priMom/pSks, Vert5_Sigma);
  priMom *= pSks_Cor/pSks;
  
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

  //LorentzVector LvK(  bMom, sqrt(KaonMass*KaonMass+bMom.mag2()) );
  //LorentzVector LvPi( sMom, sqrt(PionMass*PionMass+sMom.mag2()) );
  LorentzVector LvBeam(  bMom1, sqrt(PionMass*PionMass+bMom1.mag2()) );
  //LorentzVector LvScat( sMom, sqrt(KaonMass*KaonMass+sMom.mag2()) );
  // for (pip,pip)
  LorentzVector LvScat( priMom, sqrt(KaonMass*KaonMass+priMom.mag2()) );
  LorentzVector LvTgt( 0., 0., 0., ProtonMass );
  LorentzVector LvCore( 0., 0., 0., 0.0 );
  //LorentzVector LvTgt( 0., 0., 0., Nucl7LiMass );
  //LorentzVector LvCore( 0., 0., 0., CoreMass6Li );
  //LorentzVector LvTgt( 0., 0., 0., Nucl9BeMass );
  //LorentzVector LvCore( 0., 0., 0., CoreMass8Be );
  //   LorentzVector LvTgt( 0., 0., 0., Nucl16OMass );
  //   LorentzVector LvCore( 0., 0., 0., CoreMass15O );
  //LorentzVector LvTgt( 0., 0., 0., Nucl10BMass );
  //LorentzVector LvCore( 0., 0., 0., CoreMass9Li );
  //LorentzVector LvRc = LvK+LvTgt-LvPi;
  LorentzVector LvRc = LvBeam+LvTgt-LvScat;
  double recoilmom = LvRc.vect().mag();//(LvK.vect()-LvPi.vect()).mag();
  double mismass = LvRc.mag();//-LvCore.mag();
  double BLambda = LvRc.mag()-( LvCore.mag()+LambdaMass );
  double recoilbeta = LvRc.beta();//LvPi.beta()-LvK.beta();
  double recoilEkin = LvRc.e()-LvRc.mag();//LvPi.beta()-LvK.beta();
  event.RecoilMom = recoilmom;
  event.beta = recoilbeta;
  event.Ekin = recoilEkin;
  event.MissMass = mismass;
  event.BE = BLambda;

  ThreeVector HypBeamVec = LvRc.vect();
  event.hypBeamVec[0] = HypBeamVec.x();
  event.hypBeamVec[1] = HypBeamVec.y();
  event.hypBeamVec[2] = HypBeamVec.z();

  if( tofflag && !smfflag ) event.BE_SMF = BLambda;

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

      //std::cout << "CylAnaPi" << std::endl;
      ThreeVector Vert5(event.vertex5[0], event.vertex5[1], event.vertex5[2]);
      ThreeVector DecayPiPos = CylAnaPi.GetScatPosition();
      ThreeVector DecayPiVec = CylAnaPi.GetScatVector();

      double cdistDecayPi;
      VertDecayPi = VertexPoint3D( Vert5, DecayPiPos, HypBeamVec, DecayPiVec, cdistDecayPi);

      double costDecayPi=HypBeamVec*DecayPiVec/(HypBeamVec.mag()*DecayPiVec.mag());
      double thetaDecayPi=acos(costDecayPi)*Rad2Deg;
      double chisqrDecayPiXY = CylAnaPi.GetChisqrXY();
      double chisqrDecayPiZ = CylAnaPi.GetChisqrZ();

      event.vertexDecayPi[0] = VertDecayPi.x();
      event.vertexDecayPi[1] = VertDecayPi.y();
      event.vertexDecayPi[2] = VertDecayPi.z();

      event.decayPiVec[0] = DecayPiVec.x();
      event.decayPiVec[1] = DecayPiVec.y();
      event.decayPiVec[2] = DecayPiVec.z();

      event.decayPiPos[0] = DecayPiPos.x();
      event.decayPiPos[1] = DecayPiPos.y();
      event.decayPiPos[2] = DecayPiPos.z();

      event.cdistDecayPi = cdistDecayPi;
      event.thetaDecayPi = thetaDecayPi;
      event.chisqrDecayPiXY = chisqrDecayPiXY;
      event.chisqrDecayPiZ = chisqrDecayPiZ;

      for (int i=1; i<=FiberLayerNum; i++) {
	double edepNorm = CylAnaPi.GetFiberClusterEdepNorm(i);
	event.fClEdepNormForPi[i] = edepNorm;
	double pathLength = CylAnaPi.GetFiberPathLength(i);
	event.fClPathLengthForPi[i] = pathLength;
      }

      double p1 = HypBeamVec.mag();
      double m1 = SigmaMinusMass;
      double E1 = sqrt(p1*p1+m1*m1);
      double m2 = NeutronMass;
      double m3 = PionMass;
      double A  = (m1*m1+m3*m3-m2*m2)/2.;
      double hanbetu = (A*p1*costDecayPi)*(A*p1*costDecayPi)-(E1*E1-p1*p1*costDecayPi*costDecayPi)*(E1*E1*m3*m3-A*A);
      if (hanbetu >= 0) {
	double ans1 = (A*p1*costDecayPi+sqrt(hanbetu))/(E1*E1-p1*p1*costDecayPi*costDecayPi);
	double ans2 = (A*p1*costDecayPi-sqrt(hanbetu))/(E1*E1-p1*p1*costDecayPi*costDecayPi);
	if (ans1>=0 && ans2<0)
	  event.decayPiMomCal = ans1;	  
	else if (ans2>=0 && ans1<0)
	  event.decayPiMomCal = ans2;	  
	else if (ans1>=0 && ans2>=0) {
	  std::cout << "decayPiMomCal two answers: " << ans1 << ", " << ans2
		    << std::endl;
	  event.decayPiMomCal = -2;	  
	} else if (ans1<0 && ans2<0) {
	  std::cout << "decayPiMomCal two negative answers: " 
		    << ans1 << ", " << ans2
		    << std::endl;
	  event.decayPiMomCal = -1;	  
	}
      }
      
      if (event.decayPiMomCal>0) {
	event.decayPiMomVecCal[0] = event.decayPiMomCal*DecayPiVec.x()/DecayPiVec.mag();
	event.decayPiMomVecCal[1] = event.decayPiMomCal*DecayPiVec.y()/DecayPiVec.mag();
	event.decayPiMomVecCal[2] = event.decayPiMomCal*DecayPiVec.z()/DecayPiVec.mag();
	piDecayMomVec = (event.decayPiMomCal/DecayPiVec.mag())*DecayPiVec;
	nDecayMomVec = HypBeamVec-piDecayMomVec;

	event.decayNMomCal = nDecayMomVec.mag();
	event.decayNMomVecCal[0] = nDecayMomVec.x();
	event.decayNMomVecCal[1] = nDecayMomVec.y();
	event.decayNMomVecCal[2] = nDecayMomVec.z();
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
      //std::cout << "CylAnaP" << std::endl;
      ThreeVector Vert5(event.vertex5[0], event.vertex5[1], event.vertex5[2]);
      //ThreeVector HypBeamVec = LvRc.vect();
      ThreeVector ScatPos = CylAna.GetScatPosition();
      ThreeVector ScatVec = CylAna.GetScatVector();

      //ThreeVector Vert5(event.primaryVertex[0], event.primaryVertex[1], event.primaryVertex[2]);
      //ThreeVector HypBeamVec(event.momVectorHypBeam[0], event.momVectorHypBeam[1], event.momVectorHypBeam[2]);
      //ThreeVector ScatPos(event.scatPos[0], event.scatPos[1], event.scatPos[2]);
      //ThreeVector ScatVec(event.momVectorProtonScat[0], event.momVectorProtonScat[1], event.momVectorProtonScat[2]);

      /*
      std::cout << "Vert5 = ( " << Vert5.x() << ", " << Vert5.y() << ", "
		<< Vert5.z() << ")" << std::endl;
      std::cout << "ScatPos = ( " << ScatPos.x() << ", " << ScatPos.y() << ", "
		<< ScatPos.z() << ")" << std::endl;
      std::cout << "HypBeamVec = ( " << HypBeamVec.x() << ", " << HypBeamVec.y() << ", "
		<< HypBeamVec.z() << ")" << std::endl;
      std::cout << "ScatVec = ( " << ScatVec.x() << ", " << ScatVec.y() << ", "
		<< ScatVec.z() << ")" << std::endl;
      */
      double cdist;
      ThreeVector VertScat = VertexPoint3D( Vert5, ScatPos, HypBeamVec, ScatVec, cdist);
      double costScat=HypBeamVec*ScatVec/(HypBeamVec.mag()*ScatVec.mag());
      double thetaScat=acos(costScat)*Rad2Deg;
      double chisqrXY = CylAna.GetChisqrXY();
      double chisqrZ = CylAna.GetChisqrZ();

      event.vertexScat[0] = VertScat.x();
      event.vertexScat[1] = VertScat.y();
      event.vertexScat[2] = VertScat.z();

      event.scatVec[0] = ScatVec.x();
      event.scatVec[1] = ScatVec.y();
      event.scatVec[2] = ScatVec.z();

      event.scatPos[0] = ScatPos.x();
      event.scatPos[1] = ScatPos.y();
      event.scatPos[2] = ScatPos.z();

      event.cdistScat = cdist;
      event.thetaScat = thetaScat;
      event.chisqrXY = chisqrXY;
      event.chisqrZ = chisqrZ;

      HF1(600, VertScat.x());
      HF1(601, VertScat.y());
      HF1(602, VertScat.z());
      HF1(603, thetaScat);
      HF1(604, chisqrXY);
      HF1(605, chisqrZ);
      HF1(606, cdist);


      double M;
      double p1, E1;
      double phi = event.thetaScat;
      double p4, E4, Ekin4;
      if (event.thetaScat>0. && event.thetaScat<90.) {
	if (event.scatTarget == 0 || event.NNscatTarget == 0 || event.PiNscatTarget == 0) // proton
	  {
	    M = ProtonMass;
	    p1 = HypBeamVec.mag();
	    E1 = sqrt(SigmaMinusMass*SigmaMinusMass+p1*p1);
	    p4 = (2.*M*(E1+M)*p1*cos(phi*Deg2Rad))/
	      ((E1+M)*(E1+M)-p1*p1*cos(phi*Deg2Rad)*cos(phi*Deg2Rad));
	    E4 = sqrt(M*M+p4*p4);
	    Ekin4 = E4 - M;
	    
	    double beta = p1/(E1+M);
	    double gamma = 1./sqrt(1.-beta*beta);

	    double M3 = SigmaMinusMass;
	    double Ecm = sqrt((E1+M)*(E1+M)-p1*p1);
	    double Pcm = 1./(2.*Ecm)*
	      sqrt((Ecm-M3+M)*(Ecm-M3-M)*(Ecm+M3-M)*(Ecm+M3+M));
	    double sintCM = p4*sin(event.thetaScat*Deg2Rad)/Pcm;
	    double thetaCM = asin(sintCM)*Rad2Deg;

	    double P4cm_para = -beta*gamma*E4 + 
	      gamma*p4*cos(event.thetaScat*Deg2Rad);
	    if (P4cm_para>0)
	      thetaCM = 180. - thetaCM;

	    event.scatMomCal = p4;
	    event.scatEkinCal = Ekin4;
	    event.thetaScatCM = thetaCM;

	    // Decay pi analysis when we measure the momntum of pi
	    double fiberEdep = 0.;
	    for (int i=1; i<=FiberLayerNum; i++)
	      if (event.fClEdepForScat[i]>0.)
		fiberEdep += event.fClEdepForScat[1];

	    double totE_p = fiberEdep;
	    if (event.cClEdepForScat>0)
	      totE_p += event.cClEdepForScat;

	    double p_p = sqrt(totE_p*totE_p+2.*totE_p*ProtonMass*1000.);
	    p_p /= 1000.; // GeV/c
	    double p_cor, e_cor;
	    CorrElossOut(&p_cor, &e_cor, p_p, PROTON, ScatVec, VertScat);
	    double Ekin_cor = (e_cor-ProtonMass)*1000.;
	    totE_p = Ekin_cor;
	    ThreeVector ScatProtonMom = ScatVec*p_cor/ScatVec.mag();
	    ThreeVector ScatHypMom = HypBeamVec-ScatProtonMom;
	    if (flagFitPiXY && flagFitPiZ) {
	      //std::cout << "CylAnaPi" << std::endl;
	      ThreeVector DecayPiPos = CylAnaPi.GetScatPosition();
	      ThreeVector DecayPiVec = CylAnaPi.GetScatVector();

	      double costDecayPi2=ScatHypMom*DecayPiVec/(ScatHypMom.mag()*DecayPiVec.mag());
	      double thetaDecayPi2=acos(costDecayPi2)*Rad2Deg;
	      event.thetaDecayPi2 = thetaDecayPi2;

	      double p1 = ScatHypMom.mag();
	      double m1 = SigmaMinusMass;
	      double E1 = sqrt(p1*p1+m1*m1);
	      double m2 = NeutronMass;
	      double m3 = PionMass;
	      double A  = (m1*m1+m3*m3-m2*m2)/2.;
	      double hanbetu = (A*p1*costDecayPi2)*(A*p1*costDecayPi2)-(E1*E1-p1*p1*costDecayPi2*costDecayPi2)*(E1*E1*m3*m3-A*A);
	      if (hanbetu >= 0) {
		double ans1 = (A*p1*costDecayPi2+sqrt(hanbetu))/(E1*E1-p1*p1*costDecayPi2*costDecayPi2);
		double ans2 = (A*p1*costDecayPi2-sqrt(hanbetu))/(E1*E1-p1*p1*costDecayPi2*costDecayPi2);
		if (ans1>=0 && ans2<0)
		  event.decayPiMomCal2 = ans1;	  
		else if (ans2>=0 && ans1<0)
		  event.decayPiMomCal2 = ans2;	  
		else if (ans1>=0 && ans2>=0) {
		  std::cout << "decayPiMomCal two answers: " << ans1 << ", " << ans2
			    << std::endl;
		  event.decayPiMomCal2 = -2;	  
		} else if (ans1<0 && ans2<0) {
		  std::cout << "decayPiMomCal two negative answers: " 
			    << ans1 << ", " << ans2
			    << std::endl;
		  event.decayPiMomCal2 = -1;	  
		}
	      }
	    }
	  }
      }

      // np scattering assumption
      if (event.thetaDecayPi>0) {
	double cdistNpScat;
	ThreeVector VertNpScat = VertexPoint3D( VertDecayPi, ScatPos, nDecayMomVec, ScatVec, cdistNpScat);
	double costNpScat=nDecayMomVec*ScatVec/(nDecayMomVec.mag()*ScatVec.mag());
	double thetaNpScat=acos(costNpScat)*Rad2Deg;
	
	event.vertexNpScat[0] = VertNpScat.x();
	event.vertexNpScat[1] = VertNpScat.y();
	event.vertexNpScat[2] = VertNpScat.z();
	event.cdistNpScat = cdistNpScat;
	event.thetaNpScat = thetaNpScat;
	
	phi = event.thetaNpScat;
	
	if (event.thetaNpScat>0. && event.thetaNpScat<90.) {
	  if (event.scatTarget == 0 || event.NNscatTarget == 0 || event.PiNscatTarget == 0) // proton
	    {
	      M = ProtonMass;
	      p1 = nDecayMomVec.mag();
	      E1 = sqrt(NeutronMass*NeutronMass+p1*p1);
	      p4 = (2.*M*(E1+M)*p1*cos(phi*Deg2Rad))/
		((E1+M)*(E1+M)-p1*p1*cos(phi*Deg2Rad)*cos(phi*Deg2Rad));
	      E4 = sqrt(M*M+p4*p4);
	      Ekin4 = E4 - M;
	      
	      double beta = p1/(E1+M);
	      double gamma = 1./sqrt(1.-beta*beta);
	      
	      double M3 = NeutronMass;
	      double Ecm = sqrt((E1+M)*(E1+M)-p1*p1);
	      double Pcm = 1./(2.*Ecm)*
		sqrt((Ecm-M3+M)*(Ecm-M3-M)*(Ecm+M3-M)*(Ecm+M3+M));
	      double sintCM = p4*sin(event.thetaNpScat*Deg2Rad)/Pcm;
	      double thetaCM = asin(sintCM)*Rad2Deg;
	      
	      double P4cm_para = -beta*gamma*E4 + 
		gamma*p4*cos(event.thetaNpScat*Deg2Rad);
	      if (P4cm_para>0)
		thetaCM = 180. - thetaCM;
	      
	      event.scatNpMomCal = p4;
	      event.scatNpEkinCal = Ekin4;
	      event.thetaNpScatCM = thetaCM;
	      
	    }
	}
      }



      // pip scattering assumption
      if (event.thetaDecayPi>0) {
	ThreeVector ScatPiPos = CylAnaPi.GetScatPosition();
	ThreeVector ScatPiVec = CylAnaPi.GetScatVector();
	ScatPiVec *= event.momDecayPi/ScatPiVec.mag();
	double cdistPiPScat;
	ThreeVector VertPiPScat = VertexPoint3D( ScatPiPos, ScatPos, ScatPiVec, ScatVec, cdistPiPScat);

	double EkinProton=0.; // MeV
	for (int i=1; i<=FiberLayerNum; i++)
	  if (event.fClEdepForScat[i]>0)
	    EkinProton += event.fClEdepForScat[i];
	if (event.cClEdepForScat>0)
	  EkinProton += event.cClEdepForScat;

	EkinProton *= 0.001; // GeV
	double M_p = ProtonMass;
	double pProton = sqrt(EkinProton*EkinProton+2.*EkinProton*M_p) ;
	double pProton_cor, eProton_cor;
	CorrElossOut(&pProton_cor, &eProton_cor, pProton, PROTON, ScatVec, VertPiPScat);

	ThreeVector ScatProtonVec = ScatVec*pProton_cor/ScatVec.mag();

	ThreeVector DecayPiVec=ScatProtonVec+ScatPiVec;
	event.decayPiMom3 = DecayPiVec.mag();

	double costDecayPi3=HypBeamVec*DecayPiVec/(HypBeamVec.mag()*DecayPiVec.mag());
	double thetaDecayPi3=acos(costDecayPi3)*Rad2Deg;

	event.vertexPiPScat[0] = VertPiPScat.x();
	event.vertexPiPScat[1] = VertPiPScat.y();
	event.vertexPiPScat[2] = VertPiPScat.z();
	event.cdistPiPScat  = cdistPiPScat;
	event.thetaDecayPi3 = thetaDecayPi3;

	double p1 = HypBeamVec.mag();
	double m1 = SigmaMinusMass;
	double E1 = sqrt(p1*p1+m1*m1);
	double m2 = NeutronMass;
	double m3 = PionMass;
	double A  = (m1*m1+m3*m3-m2*m2)/2.;
	double hanbetu = (A*p1*costDecayPi3)*(A*p1*costDecayPi3)-(E1*E1-p1*p1*costDecayPi3*costDecayPi3)*(E1*E1*m3*m3-A*A);
	if (hanbetu >= 0) {
	  double ans1 = (A*p1*costDecayPi3+sqrt(hanbetu))/(E1*E1-p1*p1*costDecayPi3*costDecayPi3);
	  double ans2 = (A*p1*costDecayPi3-sqrt(hanbetu))/(E1*E1-p1*p1*costDecayPi3*costDecayPi3);
	  if (ans1>=0 && ans2<0)
	    event.decayPiMomCal3 = ans1;	  
	  else if (ans2>=0 && ans1<0)
	    event.decayPiMomCal3 = ans2;	  
	  else if (ans1>=0 && ans2>=0) {
	    std::cout << "decayPiMomCal3 two answers: " << ans1 << ", " << ans2
		      << std::endl;
	    event.decayPiMomCal3 = -2;	  
	  } else if (ans1<0 && ans2<0) {
	    std::cout << "decayPiMomCal two negative answers: " 
		      << ans1 << ", " << ans2
		      << std::endl;
	    event.decayPiMomCal3 = -1;	  
	  }
	}
      }      


      // LambdaN conversion  assumption
      if (event.thetaDecayPi>0) {
	ThreeVector DecayPiPos = CylAnaPi.GetScatPosition();
	ThreeVector DecayPiVec = CylAnaPi.GetScatVector();
	double cdistLambdaDecay;
	ThreeVector VertLambdaDecay = VertexPoint3D( DecayPiPos, ScatPos, DecayPiVec, ScatVec, cdistLambdaDecay);
	double costLambdaDecay=DecayPiVec*ScatVec/(DecayPiVec.mag()*ScatVec.mag());
	double thetaLambdaDecay=acos(costLambdaDecay)*Rad2Deg;
	
	double theta = thetaLambdaDecay;
	double cost = costLambdaDecay;

	event.cdistLambdaDecay = cdistLambdaDecay;
	event.thetaLambdaDecay = thetaLambdaDecay;
	event.vertexLambdaDecay[0] = VertLambdaDecay.x();
	event.vertexLambdaDecay[1] = VertLambdaDecay.y();
	event.vertexLambdaDecay[2] = VertLambdaDecay.z();

	/*
	std::cout << "cdistLambdaDecay = " << cdistLambdaDecay 
		  << ", thetaLambdaDecay" << thetaLambdaDecay << std::endl;
	*/


	double Ekin1=0.; // MeV
	for (int i=1; i<=FiberLayerNum; i++)
	  if (event.fClEdepForScat[i]>0)
	    Ekin1 += event.fClEdepForScat[i];

	if (event.cClEdepForScat>0)
	  Ekin1 += event.cClEdepForScat;

	Ekin1 *= 0.001; // GeV
	double M1 = ProtonMass;
	double p1 = sqrt(Ekin1*Ekin1+2.*Ekin1*M1) ;
	double p_cor, e_cor;
	CorrElossOut(&p_cor, &e_cor, p1, PROTON, ScatVec, VertLambdaDecay);

	p1 = p_cor;

	// invariant mass of lambda
	LorentzVector LvProton( ScatVec*p1/ScatVec.mag(), sqrt(ProtonMass*ProtonMass+p1*p1) );
	LorentzVector LvPi( DecayPiVec*event.momDecayPi/DecayPiVec.mag(), sqrt(PionMass*PionMass+event.momDecayPi*event.momDecayPi) );
	LorentzVector LvLambda = LvProton + LvPi;
	event.invmass_pip = LvLambda.mag();

	if (event.thetaScat>0. && event.thetaScat<180.) {
	  if (event.scatTarget == 0 || event.NNscatTarget == 0 || event.PiNscatTarget == 0) // proton
	    {
	      // Calculate momentum of pi assume Lambda decay
	      double M0 = LambdaMass;
	      double M2 = PionMass;
	      double E1 = sqrt(p1*p1+M1*M1);
	      
	      double A = (M0*M0-(M1*M1+M2*M2))/2.;
	      //std::cout << "p1 = " << p1 << ", Ekin1 = " << Ekin1 << std::endl;

	      double hanbetu = (A*p1*cost)*(A*p1*cost)-(E1*E1-p1*p1*cost*cost)*(E1*E1*M2*M2-A*A);
	      if (hanbetu>=0) {
		double ans1 = (A*p1*cost+sqrt(hanbetu))/(E1*E1-p1*p1*cost*cost);
		double ans2 = (A*p1*cost-sqrt(hanbetu))/(E1*E1-p1*p1*cost*cost);
		/*		
		std::cout << "--pi momentum" << std::endl;
		if (ans1>=0 && ans2>= 0) {
		  std::cout << "1 : ans1 = " << ans1 << ", ans2 = " << ans2 << std::endl;
		} else if (ans1>=0) {
		  std::cout << "2 : ans1 = " << ans1 << ", ans2 = " << ans2 << std::endl;
		} else if (ans2>=0) {
		  std::cout << "3 : ans1 = " << ans1 << ", ans2 = " << ans2 << std::endl;
		} else {
		  std::cout << "4 : ans1 = " << ans1 << ", ans2 = " << ans2 << std::endl;
		}
		*/
		
		if (ans1>=0) {
		  double p2 = ans1;
		  ThreeVector momVecPi     = DecayPiVec*p2/DecayPiVec.mag();
		  ThreeVector momVecProton = ScatVec*p1/ScatVec.mag();
		  ThreeVector momVecLambda = momVecPi+momVecProton;

		  event.momPiFromLambda = momVecPi.mag();
		  event.momProtonFromLambda = momVecProton.mag();
		  event.momLambda = momVecLambda.mag();
		  event.momVecPiFromLambda[0] = momVecPi.x();
		  event.momVecPiFromLambda[1] = momVecPi.y();
		  event.momVecPiFromLambda[2] = momVecPi.z();
		  event.momVecProtonFromLambda[0] = momVecProton.x();
		  event.momVecProtonFromLambda[1] = momVecProton.y();
		  event.momVecProtonFromLambda[2] = momVecProton.z();
		  event.momVecLambda[0] = momVecLambda.x();
		  event.momVecLambda[1] = momVecLambda.y();
		  event.momVecLambda[2] = momVecLambda.z();
		  /*
		  std::cout << "P_pi = "<< momVecPi.mag() << ", P_proton = " << momVecProton.mag() 
			    << ", P_lambda = " << momVecLambda.mag() << std::endl;
		  */
		  double cdistLambdaNConv;
		  ThreeVector VertLambdaNConv = VertexPoint3D( Vert5, VertLambdaDecay, HypBeamVec, momVecLambda, cdistLambdaNConv);
		  double costLambdaNConv=HypBeamVec*momVecLambda/(HypBeamVec.mag()*momVecLambda.mag());
		  double thetaLambdaNConv=acos(costLambdaNConv)*Rad2Deg;

		  event.cdistLambdaNConv = cdistLambdaNConv;
		  event.thetaLambdaNConv = thetaLambdaNConv;
		  event.vertexLambdaNConv[0] = VertLambdaNConv.x();
		  event.vertexLambdaNConv[1] = VertLambdaNConv.y();
		  event.vertexLambdaNConv[2] = VertLambdaNConv.z();

		  double m1 = SigmaMinusMass;
		  double m2 = ProtonMass;
		  double m3 = LambdaMass;
		  double m4 = NeutronMass;

		  cost = costLambdaNConv;
		  p1 = HypBeamVec.mag();
		  E1 = sqrt(p1*p1+m1*m1);
		  A = m1*m1+m2*m2+m3*m3-m4*m4+2.*E1*m2;

		  hanbetu = 4.*m3*m3*(p1*p1*cost*cost-(E1+m2)*(E1+m2))+A*A;
		  if (hanbetu>=0) {
		    ans1 = (A*p1*cost+(E1+m2)*sqrt(hanbetu))/(2.*((E1+m2)*(E1+m2)-p1*p1*cost*cost));
		    ans2 = (A*p1*cost-(E1+m2)*sqrt(hanbetu))/(2.*((E1+m2)*(E1+m2)-p1*p1*cost*cost));
		    event.momCalLambda = ans1;
		    event.momCalLambda2 = ans2;
		    /*
		    std::cout << "---Calculated Lambda mom" << std::endl;
		    if (ans1>=0 && ans2>= 0) {
		      std::cout << "1 : ans1 = " << ans1 << ", ans2 = " << ans2 << std::endl;
		    } else if (ans1>=0) {
		      std::cout << "2 : ans1 = " << ans1 << ", ans2 = " << ans2 << std::endl;
		    } else if (ans2>=0) {
		      std::cout << "3 : ans1 = " << ans1 << ", ans2 = " << ans2 << std::endl;
		    } else {
		      std::cout << "4 : ans1 = " << ans1 << ", ans2 = " << ans2 << std::endl;
		    }
		    */		    
		    double beta = p1/(E1+m2);
		    double gamma = 1./sqrt(1.-beta*beta);
		    
		    double Ecm = sqrt((E1+m2)*(E1+m2)-p1*p1);
		    double Pcm = 1./(2.*Ecm)*
		      sqrt((Ecm-m3+m4)*(Ecm-m3-m4)*(Ecm+m3-m4)*(Ecm+m3+m4));
		    double p3 = event.momCalLambda;
		    double E3 = sqrt(p3*p3+m3*m3);
		    double sintCM = p3*sin(thetaLambdaNConv*Deg2Rad)/Pcm;
		    double thetaCM = asin(sintCM)*Rad2Deg;
		    
		    double P3cm_para = -beta*gamma*E3 + 
		      gamma*p3*cos(thetaLambdaNConv*Deg2Rad);
		    if (P3cm_para<0)
		      thetaCM = 180. - thetaCM;

		    event.thetaCMLambdaNConv = thetaCM;
		  }
		}
	      }
	    }
	}
      }

      for (int i=1; i<=FiberLayerNum; i++) {
	double res = CylAna.GetResidual(i);
	HF1(610+i, res);
	double edepNorm = CylAna.GetFiberClusterEdepNorm(i);
	event.fClEdepNormForScat[i] = edepNorm;
	double pathLength = CylAna.GetFiberPathLength(i);
	event.fClPathLengthForScat[i] = pathLength;
      }
      /*
      std::cout << "VertScat = ( " << VertScat.x()
		<< ", " << VertScat.y() << ", " << VertScat.z() 
		<<"), cdist = " << cdist << ", TgtFlag = " << event.scatTarget 
		<< std::endl;
      */
    }

  }


  int nhit1=0;
  int nhit2=0;
  int nhit3=0;
  int nhit4=0;
  for (int i=0; i<ntSks; i++) {
    SksTrack *track = DCAna.GetSksTrack(i);
    int nhdc = track->GetNHits();
    for( int j=0; j<nhdc; ++j ){
      TrackHit *thit = track->GetHit(j);
      if( thit ){
	int layer = thit->GetLayer();
	double posl = thit->GetLocalHitPos();
	double resd = thit->GetResidual();
	//std::cout << "Residual=" << resd << std::endl;
	HF1( 200+layer, posl );
	HF1( 300+layer, resd );
	HF2( 400+layer, posl, resd );
	
	if( layer==1 ){
	  if (nhit1<MaxHits) {
	    event.DC1Posl[nhit1] = posl;
	    event.DC1Resd[nhit1] = resd;
	    nhit1++;
	  } else {
	    std::cerr << "Too many nhit1 " << nhit1 << std::endl;
	  }
	}
	if( layer==5 ){
	  if (nhit2<MaxHits) {
	    event.DC2Posl[nhit2] = posl;
	    event.DC2Resd[nhit2] = resd;
	    nhit2++;
	  } else {
	    std::cerr << "Too many nhit2 " << nhit2 << std::endl;
	  }
	}
	if( layer==41 ){
	  if (nhit3<MaxHits) {
	    event.DC3Posl[nhit3] = posl;
	    event.DC3Resd[nhit3] = resd;
	    nhit3++;     
	  } else {
	    std::cerr << "Too many nhit3 " << nhit3 << std::endl;
	  } 
	}
	if( layer==45 ){
	  if (nhit4<MaxHits) {
	    event.DC4Posl[nhit4] = posl;
	    event.DC4Resd[nhit4] = resd;
	    nhit4++;     
	  } else {
	    std::cerr << "Too many nhit4 " << nhit4 << std::endl;
	  }  
	}
      }
    }//std::cout<<"********************"<<std::endl;
  }
  event.DC1Hits = nhit1;
  event.DC2Hits = nhit2;
  event.DC3Hits = nhit3;
  event.DC4Hits = nhit4;

  tree->Fill();

  if( status ){
    HF1( 101, 4. );
    HF1( 123, chisqr );
    HF1( 124, nhdc ); 
    HF1( 125, rx ); HF1( 126, ry ); HF1( 127, rp );
    HF1( 128, rt ); HF1( 129, rf );
    HF1( 130, lengTof );
    HF1( 131, rx-gx ); HF1( 132, ry-gy ); HF1( 133, rp-gp );
    HF1( 134, rt-gt ); HF1( 135, rf-gf );
    HF1( 136, ru-gu ); HF1( 137, rv-gv );
    HF2( 138, gp, rp-gp );

    if( chisqr<MaxChiSqr ){
      HF1( 101, 5. );
      HF1( 163, chisqr );
      HF1( 164, nhdc ); 
      HF1( 165, rx ); HF1( 166, ry ); HF1( 167, rp );
      HF1( 168, rt ); HF1( 169, rf );
      HF1( 170, lengTof );
      HF1( 171, rx-gx ); HF1( 172, ry-gy ); HF1( 173, rp-gp );
      HF1( 174, rt-gt ); HF1( 175, rf-gf );
      HF1( 176, ru-gu ); HF1( 177, rv-gv );
      HF2( 178, gp, rp-gp );

      for( int i=0; i<nhdc; ++i ){
	TrackHit *thit = track->GetHit(i);
	if( thit ){
	  int layer = thit->GetLayer();
	  double posl = thit->GetLocalHitPos();
	  double resd = thit->GetResidual();
	  HF1( 1200+layer, posl );
	  HF1( 1300+layer, resd );
	  HF2( 1400+layer, posl, resd );
	}
      }//std::cout<<"********************"<<std::endl;
      

    }
  }
  else{
    HF1( 143, chisqr );
    HF1( 144, nhdc ); 
    HF1( 145, rx ); HF1( 146, ry ); HF1( 147, rp );
    HF1( 148, rt ); HF1( 149, rf );
    HF1( 150, lengTof );
    HF1( 151, rx-gx ); HF1( 152, ry-gy ); HF1( 153, rp-gp );
    HF1( 154, rt-gt ); HF1( 155, rf-gf );
    HF1( 156, ru-gu ); HF1( 157, rv-gv );
    HF2( 158, gp, rp-gp );
  }


//   std::cout<<"NHitIn="<< nhin <<std::endl;
//   std::cout<<"NHitOut="<< nhout <<std::endl;

  //std::cout<<"********************"<<std::endl;

  //  if( events==1000000 )   return 0;

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
    event.fClEdepForScat[i] = 0.;
    event.fClEdepNormForScat[i] = 0.;
    event.fClPathLengthForScat[i] = -999.0;

    event.fClEdepForPi[i] = 0.;
    event.fClEdepNormForPi[i] = 0.;
    event.fClPathLengthForPi[i] = -999.0;
  }

  event.cClEdepForScat = 0.0;
  event.cClEdepForPi = 0.0;
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

  HB1( 101, "stat",  20, 0., 20. );
  HB1( 102, "status",  10, 0., 10. );
  HB1( 103, "chisqr", 500, 0., 100. );
  HB1( 104, "#Hits",   30, 0., 30. );
  HB1( 105, "X0", 200, -50., 50. );
  HB1( 106, "Y0", 200, -50., 50. );
  HB1( 107, "P [GeV/c]", 100000, 0.5, 2.0 );
  HB1( 108, "Theta [rad]", 1000, 0., 3.14 );
  //HB1( 108, "Theta [degree]", 100, 0., 25. );
  HB1( 109, "Phi [degree]", 100, -200., 200. );
  HB1( 110, "PathLength [mm]", 600, 3500., 6500. );
  HB1( 111, "dX [mm]", 200, -10., 10. ); 
  HB1( 112, "dY [mm]", 200, -10., 10. ); 
  HB1( 113, "dP [GeV/c]", 500, -0.025, 0.025 );
  HB1( 114, "dTheta [deg]", 1000, -5., 5. );
  HB1( 115, "dPhi [deg]", 200, -50., 50. );
  HB1( 116, "dU", 500, -0.05, 0.05 );
  HB1( 117, "DV", 500, -0.05, 0.05 );
  HB2( 118, "dP%P [GeV/c]", 200, 0.5, 1.5, 100, -0.01, 0.01 ); 
  HB1( 119, "dP/P[10^-3]", 500, -5.0, 5.0 );

  HB1( 123, "chisqr <Success>", 500, 0., 100. );
  HB1( 124, "#Hits <Success>",   30, 0., 30. );
  HB1( 125, "X0 <Success>", 200, -50., 50. );
  HB1( 126, "Y0 <Success>", 200, -50., 50. );
  HB1( 127, "P [GeV/c] <Success>", 1000, 0.5, 1.5 );
  HB1( 128, "Theta [degree] <Success>", 100, 0., 25. );
  HB1( 129, "Phi [degree] <Success>", 100, -200., 200. );
  HB1( 130, "PathLength [mm] <Success>", 600, 3500., 6500. );
  HB1( 131, "dX [mm] <Success>", 200, -10., 10. ); 
  HB1( 132, "dY [mm] <Success>", 200, -10., 10. ); 
  HB1( 133, "dP [GeV/c] <Success>", 500, -0.025, 0.025 );
  HB1( 134, "dTheta [deg] <Success>", 1000, -5., 5. );
  HB1( 135, "dPhi [deg] <Success>", 200, -50., 50. );
  HB1( 136, "dU <Success>", 500, -0.05, 0.05 );
  HB1( 137, "DV <Success>", 500, -0.05, 0.05 );
  HB2( 138, "dP%P [GeV/c] <Success>", 200, 0.5, 1.5, 100, -0.01, 0.01 ); 

  HB1( 143, "chisqr <Fail>", 500, 0., 100. );
  HB1( 144, "#Hits <Fail>",   30, 0., 30. );
  HB1( 145, "X0 <Fail>", 200, -50., 50. );
  HB1( 146, "Y0 <Fail>", 200, -50., 50. );
  HB1( 147, "P [GeV/c] <Fail>", 1000, 0.5, 1.5 );
  HB1( 148, "Theta [degree] <Fail>", 100, 0., 25. );
  HB1( 149, "Phi [degree] <Fail>", 100, -200., 200. );
  HB1( 140, "PathLength [mm] <Fail>", 600, 3500., 6500. );
  HB1( 151, "dX [mm] <Fail>", 200, -10., 10. ); 
  HB1( 152, "dY [mm] <Fail>", 200, -10., 10. ); 
  HB1( 153, "dP [GeV/c] <Fail>", 500, -0.025, 0.025 );
  HB1( 154, "dTheta [deg] <Fail>", 1000, -5., 5. );
  HB1( 155, "dPhi [deg] <Fail>", 200, -50., 50. );
  HB1( 156, "dU <Fail>", 500, -0.05, 0.05 );
  HB1( 157, "DV <Fail>", 500, -0.05, 0.05 );
  HB2( 158, "dP%P [GeV/c] <Fail>", 200, 0.5, 1.5, 100, -0.01, 0.01 ); 

  HB1( 164, "#Hits <GoodTrack>",   30, 0., 30. );
  HB1( 165, "X0 <GoodTrack>", 200, -50., 50. );
  HB1( 166, "Y0 <GoodTrack>", 200, -50., 50. );
  HB1( 167, "P [GeV/c] <GoodTrack>", 1000, 0.5, 1.5 );
  HB1( 168, "Theta [degree] <GoodTrack>", 100, 0., 25. );
  HB1( 169, "Phi [degree] <GoodTrack>", 100, -200., 200. );
  HB1( 170, "PathLength [mm] <GoodTrack>", 600, 3500., 6500. );
  HB1( 171, "dX [mm] <GoodTrack>", 200, -10., 10. ); 
  HB1( 172, "dY [mm] <GoodTrack>", 200, -10., 10. ); 
  HB1( 173, "dP [GeV/c] <GoodTrack>", 500, -0.025, 0.025 );
  HB1( 174, "dTheta [deg] <GoodTrack>", 1000, -5., 5. );
  HB1( 175, "dPhi [deg] <GoodTrack>", 200, -50., 50. );
  HB1( 176, "dU <GoodTrack>", 500, -0.05, 0.05 );
  HB1( 177, "DV <GoodTrack>", 500, -0.05, 0.05 );
  HB2( 178, "dP%P [GeV/c] <GoodTrack>", 200, 0.5, 1.5, 100, -0.01, 0.01 ); 

  for( int i=1; i<13; ++i ){
    std::ostringstream title1, title2, title3;
    std::ostringstream title4, title5, title6;

    title1 << "HitPos [mm] Layer#" << i;
    HB1( 200+i, title1.str().c_str(), NBin1HPDcIn[i], 
	 MinHPDcIn[i], MaxHPDcIn[i] );

    title2 << "Residual [mm] Layer#" << i;
    HB1( 300+i, title2.str().c_str(), 250, -2.5, 2.5 );

    title3 << "Resd.%Pos [mm] Layer#" << i;
    HB2( 400+i, title3.str().c_str(), NBin2HPDcIn[i], 
	 MinHPDcIn[i], MaxHPDcIn[i], 100, -2.5, 2.5 );

    title4 << " <GoodTrack> Layer#" << i;
    HB1( 1200+i, title1.str().c_str(), NBin1HPDcIn[i], 
	 MinHPDcIn[i], MaxHPDcIn[i] );

    title5 << " <GoodTrack> Layer#" << i;
    HB1( 1300+i, title2.str().c_str(), 250, -2.5, 2.5 );

    title6 << " <GoodTrack> Layer#" << i;
    HB2( 1400+i, title3.str().c_str(), NBin2HPDcIn[i], 
	 MinHPDcIn[i], MaxHPDcIn[i], 100, -2.5, 2.5 );
  }
  for( int i=11; i<23; ++i ){
    std::ostringstream title1, title2, title3;
    std::ostringstream title4, title5, title6;

    title1 << "HitPos [mm] Layer#" << i;
    HB1( 200+i+30, title1.str().c_str(), NBin1HPDcOut[i-10], 
	 MinHPDcOut[i-10], MaxHPDcOut[i-10] );

    title2 << "Residual [mm] Layer#" << i;
    HB1( 300+i+30, title2.str().c_str(), 250, -2.5, 2.5 );

    title3 << "Resd.%Pos [mm] Layer#" << i;
    HB2( 400+i+30, title3.str().c_str(), NBin2HPDcOut[i-10], 
	 MinHPDcOut[i-10], MaxHPDcOut[i-10], 100, -2.5, 2.5 );

    title4 << " <GoodTrack> Layer#" << i;
    HB1( 1200+i+30, title1.str().c_str(), NBin1HPDcOut[i-10], 
	 MinHPDcOut[i-10], MaxHPDcOut[i-10] );

    title5 << " <GoodTrack> Layer#" << i;
    HB1( 1300+i+30, title2.str().c_str(), 250, -2.5, 2.5 );

    title6 << " <GoodTrack> Layer#" << i;
    HB2( 1400+i+30, title3.str().c_str(), NBin2HPDcOut[i-10], 
	 MinHPDcOut[i-10], MaxHPDcOut[i-10], 100, -2.5, 2.5 );
  }

  // Vertex
  HB1( 500, "Vertex X", 500, -100, 100 );
  HB1( 501, "Vertex Y", 500, -100, 100 );
  HB1( 502, "Vertex Z", 500, -300, 300 );
  HB1( 503, "Theta", 500, 0, 40 );
  HB2( 504, "Theta%Vertex Z", 200, -300, 300, 200, 0, 40 );

  HB1( 510, "Vertex X (w/ SKS Track)", 500, -100, 100 );
  HB1( 511, "Vertex Y (w/ SKS Track)", 500, -100, 100 );
  HB1( 512, "Vertex Z (w/ SKS Track)", 500, -300, 300 );
  HB1( 513, "Theta (w/ SKS Track)", 500, 0, 40 );
  HB2( 514, "Theta%Vertex Z (w/ SKS Track)", 200, -300, 300, 200, 0, 40 );

  HB1( 520, "Vertex X (w/ SKS Track, H-only)", 500, -100, 100 );
  HB1( 521, "Vertex Y (w/ SKS Track, H-only)", 500, -100, 100 );
  HB1( 522, "Vertex Z (w/ SKS Track, H-only)", 500, -300, 300 );
  HB1( 523, "Theta (w/ SKS Track, H-only)", 500, 0, 40 );
  HB2( 524, "Theta%Vertex Z (w/ SKS Track, H-only)", 200, -300, 300, 200, 0, 40 );

  HB1( 530, "Vertex X (w/ SKS Track)", 500, -100, 100 );
  HB1( 531, "Vertex Y (w/ SKS Track)", 500, -100, 100 );
  HB1( 532, "Vertex Z (w/ SKS Track)", 500, -300, 300 );
  HB1( 533, "Theta (w/ SKS Track)", 500, 0, 40 );
  HB2( 534, "Theta%Vertex Z (w/ SKS Track)", 200, -300, 300, 200, 0, 40 );

  HB1( 540, "Vertex X (w/ SksIn Track)", 500, -100, 100 );
  HB1( 541, "Vertex Y (w/ SksIn Track)", 500, -100, 100 );
  HB1( 542, "Vertex Z (w/ SksIn Track)", 500, -300, 300 );
  HB1( 543, "Theta (w/ SksIn Track)", 500, 0, 40 );
  HB2( 544, "Theta%Vertex Z (w/ SksIn Track)", 200, -300, 300, 200, 0, 40 );
  HB1( 545, "chisqrSksIn", 500, 0., 100. );

  HB1( 600, "Scattering Vertex X", 500, -100, 100 );
  HB1( 601, "Scattering Vertex Y", 500, -100, 100 );
  HB1( 602, "Scattering Vertex Z", 500, -300, 300 );
  HB1( 603, "Scattering Theta", 360, 0, 180 );
  HB1( 604, "chisqrXY", 500, 0., 100. );
  HB1( 605, "chisqrZ", 500, 0., 100. );
  HB1( 606, "Closest Distance", 500, 0., 100. );

  for( int i=1; i<=FiberLayerNum; ++i ){
    std::ostringstream title1;
    title1 << "Residual[mm] Layer#" << i;
    HB1( 610+i, title1.str().c_str(), 100, -2.5, 2.5 );
  }  

  // SdcIn
  HB1( 3101, "status", 10, 0., 10. );
  HB1( 3102, "chisquare", 100, 0., 10. );
  HB1( 3103, "#Hits",     20, 0., 20. );
  HB1( 3104, "X0[mm]",   2000, -50., 50. ); 
  HB1( 3105, "Y0[mm]",   2000, -50., 50. ); 
  HB1( 3106, "U0", 4000, -0.4, 0.4 );
  HB1( 3107, "V0", 3000, -0.3, 0.3 );
  HB1( 3108, "Theta[degree]", 1000, 0., 25. );
  HB1( 3109, "Phi[degree]", 1000, -200., 200. );
  HB1( 3110, "LayerDist", 20, 0., 20. );

  HB2( 3111, "U0%X0", 200, -25., 25., 200, -0.4, 0.4 );
  HB2( 3112, "V0%Y0", 200, -25., 25., 200, -0.25, 0.25 );
  HB2( 3113, "Y0%x0", 200, -25., 25., 200, -25., 25. ); 
  HB2( 3114, "Phi%Theta", 100, 0., 25., 100, -200., 200. );

  HB1( 3115, "dX [mm]", 200, -10., 10. ); 
  HB1( 3116, "dY [mm]", 200, -10., 10. ); 
  HB1( 3117, "dU", 5000, -0.05, 0.05 );
  HB1( 3118, "dV", 5000, -0.05, 0.05 );

  for( int i=1; i<=12; ++i ){
    std::ostringstream title1, title2, title3;
    title1 << "HitPos[mm] Layer#" << i;
    HB1( 3120+i, title1.str().c_str(), 500, -250., 250. );
    title2 << "Residual[mm] Layer#" << i;
    HB1( 3140+i, title2.str().c_str(), 100, -2.5, 2.5 );
    title3 << "Resd.%Pos[mm] Layer#" << i;
    HB2( 3160+i, title3.str().c_str(), 100, -250., 250., 100, -2.5, 2.5 );
  }

  // SdcOut
  HB1( 4101, "status",  10, 0., 10. );
  HB1( 4102, "chisquare", 100, 0., 10. );
  HB1( 4103, "#Hits", 20, 0., 20. );
  HB1( 4104, "X0[mm]",  6000, -600., 600. );
  HB1( 4105, "Y0[mm]",  6000, -600., 600. );
  HB1( 4106, "U0", 4000, -0.4, 0.4 );
  HB1( 4107, "V0", 3000, -0.3, 0.3 );
  HB1( 4108, "Theta[degree]", 1000, 0., 25. );
  HB1( 4109, "Phi[degree]", 1000, -200., 200. );
  HB1( 4110, "LayerDist", 20, 30., 50. );
       
  HB2( 4111, "U0%X0", 200, -600., 600., 200, -0.4, 0.4 );
  HB2( 4112, "V0%Y0", 200, -600., 600., 200, -0.25, 0.25 );
  HB2( 4113, "Y0%x0", 200, -600., 600., 200, -600., 600. ); 
  HB2( 4114, "Phi%Theta", 100, 0., 25., 100, -200., 200. );
       
  for( int i=41; i<=52; ++i ){
    std::ostringstream title1, title2, title3;
    title1 << "HitPos[mm] Layer#" << i;
    HB1( 4090+i, title1.str().c_str(), 1200, -1200., 1200. );
    title2 << "Residual[mm] Layer#" << i;
    HB1( 4110+i, title2.str().c_str(), 100, -2.5, 2.5 );
    title3 << "Resd.%Pos[mm] Layer#" << i;
    HB2( 4130+i, title3.str().c_str(), 400, -1200., 1200., 100, -2.5, 2.5 );
  }  

  // BdcOut
  HB1( 5101, "status", 10, 0., 10. );
  HB1( 5102, "chisquare", 100, 0., 10. );
  HB1( 5103, "#Hits",     20, 0., 20. );
  HB1( 5104, "X0[mm]",   2000, -50., 50. ); 
  HB1( 5105, "Y0[mm]",   2000, -50., 50. ); 
  HB1( 5106, "U0", 4000, -0.4, 0.4 );
  HB1( 5107, "V0", 3000, -0.3, 0.3 );
  HB1( 5108, "Theta[degree]", 1000, 0., 25. );
  HB1( 5109, "Phi[degree]", 1000, -200., 200. );
  HB1( 5110, "LayerDist", 20, 0., 20. );

  HB2( 5111, "U0%X0", 200, -25., 25., 200, -0.4, 0.4 );
  HB2( 5112, "V0%Y0", 200, -25., 25., 200, -0.25, 0.25 );
  HB2( 5113, "Y0%x0", 200, -25., 25., 200, -25., 25. ); 
  HB2( 5114, "Phi%Theta", 100, 0., 25., 100, -200., 200. );

  HB1( 5115, "dX [mm]", 200, -10., 10. ); 
  HB1( 5116, "dY [mm]", 200, -10., 10. ); 
  HB1( 5117, "dU", 5000, -0.05, 0.05 );
  HB1( 5118, "dV", 5000, -0.05, 0.05 );

  for( int i=1; i<=14; ++i ){
    std::ostringstream title1, title2, title3;
    title1 << "HitPos[mm] Layer#" << i;
    HB1( 5120+i, title1.str().c_str(), 500, -250., 250. );
    title2 << "Residual[mm] Layer#" << i;
    HB1( 5140+i, title2.str().c_str(), 100, -2.5, 2.5 );
    title3 << "Resd.%Pos[mm] Layer#" << i;
    HB2( 5160+i, title3.str().c_str(), 100, -250., 250., 100, -2.5, 2.5 );
  }


  // SksIn
  HB1( 6101, "status", 10, 0., 10. );
  HB1( 6102, "chisquare", 100, 0., 10. );
  HB1( 6103, "#Hits",     20, 0., 20. );
  HB1( 6104, "X0[mm]",   2000, -50., 50. ); 
  HB1( 6105, "Y0[mm]",   2000, -50., 50. ); 
  HB1( 6106, "U0", 4000, -0.4, 0.4 );
  HB1( 6107, "V0", 3000, -0.3, 0.3 );
  HB1( 6108, "Theta[degree]", 1000, 0., 25. );
  HB1( 6109, "Phi[degree]", 1000, -200., 200. );
  HB1( 6110, "LayerDist", 20, 0., 20. );

  HB2( 6111, "U0%X0", 200, -25., 25., 200, -0.4, 0.4 );
  HB2( 6112, "V0%Y0", 200, -25., 25., 200, -0.25, 0.25 );
  HB2( 6113, "Y0%x0", 200, -25., 25., 200, -25., 25. ); 
  HB2( 6114, "Phi%Theta", 100, 0., 25., 100, -200., 200. );

  HB1( 6115, "dX [mm]", 200, -10., 10. ); 
  HB1( 6116, "dY [mm]", 200, -10., 10. ); 
  HB1( 6117, "dU", 5000, -0.05, 0.05 );
  HB1( 6118, "dV", 5000, -0.05, 0.05 );

  for( int i=1; i<=12; ++i ){
    std::ostringstream title1, title2, title3;
    title1 << "HitPos[mm] Layer#" << i;
    HB1( 6120+i, title1.str().c_str(), 500, -250., 250. );
    title2 << "Residual[mm] Layer#" << i;
    HB1( 6140+i, title2.str().c_str(), 100, -2.5, 2.5 );
    title3 << "Resd.%Pos[mm] Layer#" << i;
    HB2( 6160+i, title3.str().c_str(), 100, -250., 250., 100, -2.5, 2.5 );
  }


  //Tree
  HBTree("tree","tree of Sks");
  TTree *tree = dynamic_cast<TTree *>(gFile->Get("tree"));

  //Primary for Yd scattering
  tree->Branch("thetaMeson",   &event.thetaMeson,  "thetaMeson/D");
  tree->Branch("phiMeson",   &event.phiMeson,  "phiMeson/D");
  tree->Branch("thetaMesonCM",   &event.thetaMesonCM,  "thetaMesonCM/D");
  tree->Branch("phiMesonCM",   &event.phiMesonCM,  "phiMesonCM/D");
  tree->Branch("thetaScatHypCM",   &event.thetaScatHypCM,  "thetaScatHypCM/D");
  tree->Branch("phiScatHypCM",   &event.phiScatHypCM,  "phiScatHypCM/D");
  tree->Branch("momVectorScatMeson",   event.momVectorScatMeson,  "momVectorScatMeson[3]/D");
  tree->Branch("momVectorHypBeam",   event.momVectorHypBeam,  "momVectorHypBeam[3]/D");
  tree->Branch("momVectorHypScat",   event.momVectorHypScat,  "momVectorHypScat[3]/D");
  tree->Branch("momVectorProtonScat",   event.momVectorProtonScat,  "momVectorProtonScat[3]/D");
  tree->Branch("momVectorDecayPi",   event.momVectorDecayPi,  "momVectorDecayPi[3]/D");
  tree->Branch("momVectorDecayNucleon",   event.momVectorDecayNucleon,  "momVectorDecayNucleon[3]/D");
  tree->Branch("momScatMeson",   &event.momScatMeson,  "momScatMeson/D");
  tree->Branch("momHypBeam",   &event.momHypBeam,  "momHypBeam/D");
  tree->Branch("momHypScat",   &event.momHypScat,  "momHypScat/D");
  tree->Branch("momProtonScat",   &event.momProtonScat,  "momProtonScat/D");
  tree->Branch("momDecayPi",   &event.momDecayPi,  "momDecayPi/D");
  tree->Branch("momDecayNucleon",   &event.momDecayNucleon,  "momDecayNucleon/D");
  tree->Branch("primaryVertex",   event.primaryVertex,  "primaryVertex[3]/D");
  tree->Branch("scatPos0",   event.scatPos0,  "scatPos0[3]/D");
  tree->Branch("NNscatPos",   event.NNscatPos,  "NNscatPos[3]/D");
  tree->Branch("PiNscatPos",   event.PiNscatPos,  "PiNscatPos[3]/D");
  tree->Branch("decayPos",   event.decayPos,  "decayPos[3]/D");
  tree->Branch("decayFlag",   &event.decayFlag,  "decayFlag/I");
  tree->Branch("scatFlag",   &event.scatFlag,  "scatFlag/I");
  tree->Branch("scatTarget",   &event.scatTarget,  "scatTarget/I");
  tree->Branch("NNscatFlag",   &event.NNscatFlag,  "NNscatFlag/I");
  tree->Branch("NNscatTarget",   &event.NNscatTarget,  "NNscatTarget/I");
  tree->Branch("PiNscatFlag",   &event.PiNscatFlag,  "PiNscatFlag/I");
  tree->Branch("PiNscatTarget",   &event.PiNscatTarget,  "PiNscatTarget/I");

  tree->Branch("vertex1",   event.vertex1,  "vertex1[3]/D");
  tree->Branch("theta1",   &event.theta1,  "theta1/D");
  tree->Branch("cdist1",   &event.cdist1,  "cdist1/D");

  tree->Branch("u0BcOut",   &event.u0BcOut,  "u0BcOut/D");
  tree->Branch("v0BcOut",   &event.v0BcOut,  "v0BcOut/D");
  tree->Branch("u0SdcIn",   &event.u0SdcIn,  "u0SdcIn/D");
  tree->Branch("v0SdcIn",   &event.v0SdcIn,  "v0SdcIn/D");

  tree->Branch("vertex2",   event.vertex2,  "vertex2[3]/D");
  tree->Branch("theta2",   &event.theta2,  "theta2/D");
  tree->Branch("cdist2",   &event.cdist2,  "cdist2/D");

  tree->Branch("vertex3",   event.vertex3,  "vertex3[3]/D");

  tree->Branch("vertex4",   event.vertex4,  "vertex4[3]/D");
  tree->Branch("theta4",   &event.theta4,  "theta4/D");
  tree->Branch("cdist4",   &event.cdist4,  "cdist4/D");

  tree->Branch("priMomVec",   event.priMomVec,  "priMomVec[3]/D");
  tree->Branch("priPosVec",   event.priPosVec,  "priPosVec[3]/D");

  tree->Branch("priMomSksInVec",   event.priMomSksInVec,  "priMomSksInVec[3]/D");
  tree->Branch("priPosSksInVec",   event.priPosSksInVec,  "priPosSksInVec[3]/D");
  tree->Branch("vertex5",   event.vertex5,  "vertex5[3]/D");
  tree->Branch("theta5",   &event.theta5,  "theta5/D");
  tree->Branch("cdist5",   &event.cdist5,  "cdist5/D");

  tree->Branch("chisqrSksIn",   &event.chisqrSksIn,  "chisqrSksIn/D");

  tree->Branch("vertexScat",   event.vertexScat,  "vertexScat[3]/D");
  tree->Branch("scatVec",   event.scatVec,  "scatVec[3]/D");
  tree->Branch("scatPos",   event.scatPos,  "scatPos[3]/D");
  tree->Branch("hypBeamVec",   event.hypBeamVec,  "hypBeamVec[3]/D");
  tree->Branch("cdistScat",   &event.cdistScat,  "cdistScat/D");
  tree->Branch("thetaScat",   &event.thetaScat,  "thetaScat/D");
  tree->Branch("chisqrXY",   &event.chisqrXY,  "chisqrXY/D");
  tree->Branch("chisqrZ",   &event.chisqrZ,  "chisqrZ/D");
  tree->Branch("scatMomCal",   &event.scatMomCal,  "scatMomCal/D");
  tree->Branch("scatEkinCal",   &event.scatEkinCal,  "scatEkinCal/D");
  tree->Branch("thetaScatCM",   &event.thetaScatCM,  "thetaScatCM/D");

  tree->Branch("vertexDecayPi",   event.vertexDecayPi,  "vertexDecayPi[3]/D");
  tree->Branch("decayPiVec",   event.decayPiVec,  "decayPiVec[3]/D");
  tree->Branch("decayPiPos",   event.decayPiPos,  "decayPiPos[3]/D");
  tree->Branch("cdistDecayPi",   &event.cdistDecayPi,  "cdistDecayPi/D");
  tree->Branch("thetaDecayPi",   &event.thetaDecayPi,  "thetaDecayPi/D");
  tree->Branch("chisqrDecayPiXY",   &event.chisqrDecayPiXY,  "chisqrDecayPiXY/D");
  tree->Branch("chisqrDecayPiZ",   &event.chisqrDecayPiZ,  "chisqrDecayPiZ/D");
  tree->Branch("decayPiMomCal",   &event.decayPiMomCal,  "decayPiMomCal/D");
  tree->Branch("decayPiMomVecCal",   event.decayPiMomVecCal,  "decayPiMomVecCal/D");
  tree->Branch("decayNMomCal",   &event.decayNMomCal,  "decayNMomCal/D");
  tree->Branch("decayNMomVecCal",   event.decayNMomVecCal,  "decayNMomVecCal[3]/D");

  tree->Branch("thetaDecayPi2",   &event.thetaDecayPi2,  "thetaDecayPi2/D");
  tree->Branch("decayPiMomCal2",   &event.decayPiMomCal2,  "decayPiMomCal2/D");

  tree->Branch("vertexPiPScat",   event.vertexPiPScat,  "vertexPiPScat[3]/D");
  tree->Branch("cdistPiPScat",   &event.cdistPiPScat,  "cdistPiPScat/D");
  tree->Branch("thetaDecayPi3",   &event.thetaDecayPi3,  "thetaDecayPi3/D");
  tree->Branch("decayPiMom3",   &event.decayPiMom3,  "decayPiMom3/D");
  tree->Branch("decayPiMomCal3",   &event.decayPiMomCal3,  "decayPiMomCal3/D");

  tree->Branch("vertexNpScat",   event.vertexNpScat,  "vertexNpScat[3]/D");
  tree->Branch("cdistNpScat",   &event.cdistNpScat,  "cdistNpScat/D");
  tree->Branch("thetaNpScat",   &event.thetaNpScat,  "thetaNpScat/D");
  tree->Branch("scatNpMomCal",   &event.scatNpMomCal,  "scatNpMomCal/D");
  tree->Branch("scatNpEkinCal",   &event.scatNpEkinCal,  "scatNpEkinCal/D");
  tree->Branch("thetaNpScatCM",   &event.thetaNpScatCM,  "thetaNpScatCM/D");

  tree->Branch("cdistLambdaDecay",   &event.cdistLambdaDecay,  "cdistLambdaDecay/D");
  tree->Branch("thetaLambdaDecay",   &event.thetaLambdaDecay,  "thetaLambdaDecay/D");
  tree->Branch("vertexLambdaDecay",   event.vertexLambdaDecay,  "vertexLambdaDecay[3]/D");
  tree->Branch("momPiFromLambda",   &event.momPiFromLambda,  "momPiFromLambda/D");
  tree->Branch("momProtonFromLambda",   &event.momProtonFromLambda,  "momProtonFromLambda/D");
  tree->Branch("momLambda",   &event.momLambda,  "momLambda/D");
  tree->Branch("momVecPiFromLambda",   event.momVecPiFromLambda,  "momVecPiFromLambda/D");
  tree->Branch("momVecProtonFromLambda",   event.momVecProtonFromLambda,  "momVecProtonFromLambda/D");
  tree->Branch("momVecLambda",   event.momVecLambda,  "momVecLambda/D");
  tree->Branch("cdistLambdaNConv",   &event.cdistLambdaNConv,  "cdistLambdaNConv/D");
  tree->Branch("thetaLambdaNConv",   &event.thetaLambdaNConv,  "thetaLambdaNConv/D");
  tree->Branch("vertexLambdaNConv",   event.vertexLambdaNConv,  "vertexLambdaNConv[3]/D");
  tree->Branch("momCalLambda",   &event.momCalLambda,  "momCalLambda/D");
  tree->Branch("momCalLambda2",   &event.momCalLambda2,  "momCalLambda2/D");
  tree->Branch("thetaCMLambdaNConv",   &event.thetaCMLambdaNConv,  "thetaCMLambdaNConv/D");
  tree->Branch("invmass_pip",   &event.invmass_pip,  "invmass_pip/D");

  //Primary
  tree->Branch("Hits",   &event.Hits,  "Hits/I");
  tree->Branch("chisqr", &event.chisqr,"chisqr/D");
  tree->Branch("P",      &event.P,     "P/D");
  tree->Branch("Theta",  &event.Theta, "Theta/D");
  tree->Branch("Phi",    &event.Phi,   "Phi/D");
  tree->Branch("X0",     &event.X0,    "X0/D");
  tree->Branch("Y0",     &event.Y0,    "Y0/D");
  tree->Branch("U0",     &event.U0,    "U0/D");
  tree->Branch("V0",     &event.V0,    "V0/D");
  tree->Branch("path",   &event.path,  "path/D");
  tree->Branch("Pb",     &event.Pb,    "Pb/D");
  tree->Branch("Ub",     &event.Ub,    "Ub/D");
  tree->Branch("Vb",     &event.Vb,    "Vb/D");

  tree->Branch("dP",     &event.dP,    "dP/D");
  tree->Branch("dTheta", &event.dTheta,"dTheta/D");
  tree->Branch("dPhi",   &event.dPhi,  "dPhi/D");
  tree->Branch("dX",     &event.dX,    "dX/D");
  tree->Branch("dY",     &event.dY,    "dY/D");
  tree->Branch("dU",     &event.dU,    "dU/D");
  tree->Branch("dV",     &event.dV,    "dV/D");

  tree->Branch("m2",        &event.m2,       "m2/D");
  tree->Branch("RecoilMom", &event.RecoilMom,"RecoilMom/D");
  tree->Branch("beta",      &event.beta,     "beta/D");
  tree->Branch("Ekin",      &event.Ekin,     "Ekin/D");
  tree->Branch("MissmMass", &event.MissMass, "MissMass/D");
  tree->Branch("BE",        &event.BE,       "BE/D");
  tree->Branch("BE_SMF",    &event.BE_SMF,   "BE_SMF/D");

  tree->Branch("DC1Hits", &event.DC1Hits, "DC1Hits/I");
  tree->Branch("DC2Hits", &event.DC2Hits, "DC2Hits/I");
  tree->Branch("DC3Hits", &event.DC3Hits, "DC3Hits/I");
  tree->Branch("DC4Hits", &event.DC4Hits, "DC4Hits/I");
  tree->Branch("DC1Posl", event.DC1Posl,"DC1Posl[DC1Hits]/D");
  tree->Branch("DC1Resd", event.DC1Resd,"DC1Resd[DC1Hits]/D");
  tree->Branch("DC2Posl", event.DC2Posl,"DC2Posl[DC2Hits]/D");
  tree->Branch("DC2Resd", event.DC2Resd,"DC2Resd[DC2Hits]/D");
  tree->Branch("DC3Posl", event.DC3Posl,"DC3Posl[DC3Hits]/D");
  tree->Branch("DC3Resd", event.DC3Resd,"DC3Resd[DC3Hits]/D");
  tree->Branch("DC4Posl", event.DC4Posl,"DC4Posl[DC4Hits]/D");
  tree->Branch("DC4Resd", event.DC4Resd,"DC4Resd[DC4Hits]/D");

  tree->Branch("tofHits", &event.tofHits, "tofHits/I");
  tree->Branch("smfHits", &event.smfHits, "smfHits/I");
  tree->Branch("tofTime", event.tofTime,"tofTime[tofHits]/D");
  tree->Branch("tofEdep", event.tofEdep,"tofEdep[tofHits]/D");
  tree->Branch("smfTime", event.smfTime,"smfTime[smfHits]/D");
  tree->Branch("smfEdep", event.smfEdep,"smfEdep[smfHits]/D");

  char buf1[100], buf2[100];
  for (int i=1; i<=FiberLayerNum; i++) {
    sprintf(buf1, "Fiber%dHits", i);
    sprintf(buf2, "Fiber%dHits/I", i);
    tree->Branch(buf1, &event.FiberHits[i],  buf2);
  }

  for (int i=1; i<=FiberLayerNum; i++) {
    sprintf(buf1, "Fiber%dSeg", i);
    sprintf(buf2, "Fiber%dSeg[Fiber%dHits]/I", i, i);
    tree->Branch(buf1,  event.FiberSeg[i],  buf2);
    sprintf(buf1, "Fiber%dTime", i);
    sprintf(buf2, "Fiber%dTime[Fiber%dHits]/D", i, i);
    tree->Branch(buf1, event.FiberTime[i], buf2);
    sprintf(buf1, "Fiber%dEdep", i);
    sprintf(buf2, "Fiber%dEdep[Fiber%dHits]/D", i, i);
    tree->Branch(buf1, event.FiberEdep[i], buf2);
    sprintf(buf1, "Fiber%dPID", i);
    sprintf(buf2, "Fiber%dPID[Fiber%dHits]/I", i, i);
    tree->Branch(buf1, event.FiberPID[i], buf2);
  }

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


  //Crystal
  for (int i=1; i<=CrystalLayerNum; i++) {
    sprintf(buf1, "Crystal%dHits", i);
    sprintf(buf2, "Crystal%dHits/I", i);
    tree->Branch(buf1, &event.CrystalHits[i],  buf2);
  }

  for (int i=1; i<=CrystalLayerNum; i++) {
    sprintf(buf1, "Crystal%dSeg", i);
    sprintf(buf2, "Crystal%dSeg[Crystal%dHits]/I", i, i);
    tree->Branch(buf1,  event.CrystalSeg[i],  buf2);
    sprintf(buf1, "Crystal%dTime", i);
    sprintf(buf2, "Crystal%dTime[Crystal%dHits]/D", i, i);
    tree->Branch(buf1, event.CrystalTime[i], buf2);
    sprintf(buf1, "Crystal%dEdep", i);
    sprintf(buf2, "Crystal%dEdep[Crystal%dHits]/D", i, i);
    tree->Branch(buf1, event.CrystalEdep[i], buf2);
    sprintf(buf1, "Crystal%dPID", i);
    sprintf(buf2, "Crystal%dPID[Crystal%dHits]/I", i, i);
    tree->Branch(buf1, event.CrystalPID[i], buf2);
  }

  tree->Branch("cClEdepForScat", &event.cClEdepForScat, "cClEdepForScat/D");
  tree->Branch("cClEdepForPi", &event.cClEdepForPi, "cClEdepForPi/D");
  tree->Branch("cRcEdepForScat", &event.cRcEdepForScat, "cRcEdepForScat/D");
  tree->Branch("cRcEdepForPi", &event.cRcEdepForPi, "cRcEdepForPi/D");

}
