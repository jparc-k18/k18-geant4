/*
  SksTracking.cc
*/

#include "DCLocalTrack.hh"
#include "DCHit.hh"
#include "DCConfMan.hh"
#include "K18Track.hh"
#include "K18TransMatrix.hh"
#include "TrackHit.hh"
#include "ThreeVector.hh"
#include "LorentzVector.hh"
#include "Kinematics.hh"
#include "RadDeg.hh"
#include "DCAnalyzer.hh" 

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
  double x0In0;
  double y0In0;
  double u0In0;
  double v0In0;
  double theta0;
  double p0;
  double x0In;
  double y0In;
  double u0In;
  double v0In;
  double p;
  double thetaIn;
  double x0;
  double y0;
  double u0;
  double v0;
  double thetat;
  double chisqr;
  double mm1;
  double mm2;
  double bind;
};
static Event event;

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
	      << " geomFile matrixFile RootFile dataFile"
	      << std::endl;
    exit(-1);
  }

  DCConfMan::GetInstance().Initialize(argv[1]);
  K18TransMatrix *K18Matrix = new K18TransMatrix(argv[2]);
  K18Matrix->Initialize();
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

  const int maxMultiBcIn = 3;
  const int maxMultiBcOut = 3;

  static int events=0;
  events++;
  if (events%100==0)
    std::cout<<"Events:"<< events<<std::endl;


  //Tree
  TTree *tree = dynamic_cast<TTree *>(gFile->Get("tree"));
  InitializeEvent();

  if( !In ) return false;

  In >> event.x0In0 >> event.y0In0 >> event.u0In0 >> event.v0In0 >> event.p0;
  event.theta0=atan(sqrt(event.u0In0*event.u0In0 + event.v0In0*event.v0In0));

  DCAnalyzer DCAna;

  int lnum;
  In >> lnum;
  volatile int lnumold=lnum;

  int multiBcIn[13];
  for (int i=0; i<13; i++) 
    multiBcIn[i] = 0;
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
    else {
      std::cout << "Strange lnum : " << lnum << std::endl;
      return true;
    }

    if(( lnumold>=101 && lnumold <=124) || lnumold ==131 || lnumold ==132){
      DCAna.SetDCHits(lnumold, pos);
      if (lnumold>=101 && lnumold<=112)
	multiBcIn[lnumold-100]++;
      else if (lnumold>=113 && lnumold<=124)
	multiBcOut[lnumold-112]++;
    }
    lnumold = lnum;
  }

  HF1( 1, 1. );

  for (int i=1; i<=12; i++) {
    //std::cout << "MultiBcIn[" << i << "] = " << multiBcIn[i] << std::endl;
    if (multiBcIn[i] > maxMultiBcIn) {
      std::cerr << "Too many multi: BcIn[" << i << "], " << multiBcIn[i]
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

  DCAna.TrackSearchBcIn();
  int ntBcIn = DCAna.GetNtracksBcIn();
  HF1( 10, ntBcIn );
  //std::cout << "ntBcIn = " << ntBcIn << std::endl;
  for (int i=0; i<ntBcIn; i++) {
    DCLocalTrack *trBcIn = DCAna.GetTrackBcIn(i);
    bool stBcIn=trBcIn->GetStatus();
    int nhBcIn=trBcIn->GetNHit();
    double chisqrBcIn=trBcIn->GetChiSquare();
    double x0BcIn=trBcIn->GetX0(), y0BcIn=trBcIn->GetY0();
    double u0BcIn=trBcIn->GetU0(), v0BcIn=trBcIn->GetV0();
    //std::cout << "chisqrBcIn = " << chisqrBcIn << std::endl;
    HF1( 11, nhBcIn );
    HF1( 12, chisqrBcIn );
    HF1( 13, x0BcIn );
    HF1( 14, y0BcIn );
    HF1( 15, u0BcIn );
    HF1( 16, v0BcIn );
    event.x0In = x0BcIn;
    event.y0In = y0BcIn;
    event.u0In = u0BcIn;
    event.v0In = v0BcIn;
    event.thetaIn=atan(sqrt(event.u0In*event.u0In + event.v0In*event.v0In));
  }

  DCAna.TrackSearchBcOut();
  int ntBcOut = DCAna.GetNtracksBcOut();
  HF1( 30, ntBcOut );
  //std::cout << "ntBcOut = " << ntBcOut << std::endl;
  for (int i=0; i<ntBcOut; i++) {
    DCLocalTrack *trBcOut = DCAna.GetTrackBcOut(i);
    bool stBcOut=trBcOut->GetStatus();
    int nhBcOut=trBcOut->GetNHit();
    double chisqrBcOut=trBcOut->GetChiSquare();
    double x0BcOut=trBcOut->GetX0(), y0BcOut=trBcOut->GetY0();
    double u0BcOut=trBcOut->GetU0(), v0BcOut=trBcOut->GetV0();
    //std::cout << "chisqrBcOut = " << chisqrBcOut << std::endl;
    HF1( 31, nhBcOut );
    HF1( 32, chisqrBcOut );
    HF1( 33, x0BcOut );
    HF1( 34, y0BcOut );
    HF1( 35, u0BcOut );
    HF1( 36, v0BcOut );
  }

  if (ntBcIn>=1 && ntBcOut>=1) {
    DCAna.TrackSearchK18();

    double mhy=11.576;//11B + xi
    double mt=0.93146*12;//carbon target
    double mb=0.493677;
    double ms=0.493677;
    double pb=1.8;
    double ps=1.3;
    double eb=sqrt(mb*mb + pb*pb);
    double es=sqrt(ms*ms + ps*ps);

    int ntK18 = DCAna.GetNTracksK18();
    //std::cout << "ntK18 = " << ntK18 << std::endl;
    HF1( 50, ntK18 );
    for( int i=0; i<ntK18; ++i ){
      K18Track *tp=DCAna.GetK18Track(i);
      if(!tp) continue;
      int nh=tp->GetNHitsTotal();
      double chisqr=tp->chisquare();
      //std::cout << "chisqrK18=" << chisqr << std::endl;
      double xi=tp->Xin(), yi=tp->Yin();
      double ui=tp->Uin(), vi=tp->Vin();
      double xo=tp->Xout(), yo=tp->Yout();
      double uo=tp->Uout(), vo=tp->Vout();
      double xt=tp->Xtgt(), yt=tp->Ytgt();
      double ut=tp->Utgt(), vt=tp->Vtgt();
      double p=tp->P();
      double p2=tp->PD2U();

      event.p = p;
      ps = p;
      es=sqrt(ms*ms + ps*ps);

      event.thetat=atan(sqrt(ui*ui + vi*vi));

      event.mm1 = sqrt(mt*mt + mb*mb + ms*ms + 2*(mt*eb - eb*es - mt*es + pb*ps*cos(event.thetaIn)));
      event.mm2 = sqrt(mt*mt + mb*mb + ms*ms + 2*(mt*eb - eb*es - mt*es + pb*ps*cos(event.thetat)));
      event.bind = mhy - event.mm1;

      //      std::cout<<"delta="<<tp->Delta()<<std::endl;
      //      std::cout<<"p="<<p<<std::endl;
      //      std::cout<<"p0="<<tp->P()/(1+0.01*tp->Delta())<<std::endl;

      HF1( 51, nh );      
      HF1( 52, chisqr );      
      HF1( 53, xi );      
      HF1( 54, yi );      
      HF1( 55, ui );      
      HF1( 56, vi );      
      HF1( 57, xo );      
      HF1( 58, yo );      
      HF1( 59, uo );      
      HF1( 60, vo );
      HF1( 61, xt );      
      HF1( 62, yt );      
      HF1( 63, ut );      
      HF1( 64, vt );
      event.x0 = xi;
      event.y0 = yi;
      event.u0 = ui;
      event.v0 = vi;
      event.chisqr = chisqr;
      HF1( 65, p );      
      HF1( 66, event.p0-p );      
      HF1( 67, p2 );      
      HF1( 68, event.p0-p2 );      
      HF2( 69, event.p0, event.p0-p );      
      HF2( 70, event.p0, event.p0-p2 );      
    }
  }
  tree->Fill();

  return true;

}

void InitializeEvent( void )
{

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

  HB1( 1, "stat",  20, 0., 20. );

  HB1(10, "ntDcIn", 10, 0., 10);
  HB1(11, "nhDcIn", 15, 0., 15);
  HB1(12, "chisqrDcIn", 100, 0., 10);
  HB1(13, "x0DcIn", 100, -50., 50);
  HB1(14, "y0DcIn", 60, -30., 30);
  HB1(15, "u0DcIn", 200, -1., 1);
  HB1(16, "v0DcIn", 200, -1., 1);

  HB1(30, "ntDcOut", 10, 0., 10);
  HB1(31, "nhDcOut", 15, 0., 15);
  HB1(32, "chisqrDcOut", 100, 0., 10);
  HB1(33, "x0DcOut", 1200, -600., 600);
  HB1(34, "y0DcOut", 600, -300., 300);
  HB1(35, "u0DcOut", 200, -1., 1);
  HB1(36, "v0DcOut", 200, -1., 1);

  HB1(50, "ntS2s", 10, 0., 10);
  HB1(51, "nhS2s", 15, 0., 15);
  HB1(52, "chisqrS2s", 100, 0., 10);
  HB1(53, "xiS2s", 100, -50., 50);
  HB1(54, "yiS2s", 60, -30., 30);
  HB1(55, "uiS2s", 200, -1., 1);
  HB1(56, "viS2s", 200, -1., 1);
  HB1(57, "xoS2s", 1200, -600., 600);
  HB1(58, "yoS2s", 600, -300., 300);
  HB1(59, "uoS2s", 200, -1., 1);
  HB1(60, "voS2s", 200, -1., 1);
  HB1(61, "xtS2s", 200, -100., 100);
  HB1(62, "ytS2s", 200, -100., 100);
  HB1(63, "utS2s", 200, -1., 1);
  HB1(64, "vtS2s", 200, -1., 1);
  HB1(65, "p", 3000, 1.15, 1.45);
  HB1(66, "#Delta p", 200, -0.01, 0.01);
  HB1(67, "p (Calc Xin information)", 3000, 1.15, 1.45);
  HB1(68, "#Delta p (Calc Xin information)", 1100, -10.5, 10.5);
  HB2(69, "#Delta p:p0", 600, 1.15, 1.45, 400, -0.1, 0.1);
  HB2(70, "#Delta p2:p0", 600, 1.15, 1.45, 400, -0.1, 0.1);
  //Tree
  HBTree("tree","tree of Sks");
  TTree *tree = dynamic_cast<TTree *>(gFile->Get("tree"));
  tree->Branch("x0In0",   &event.x0In0,  "x0In0/D");
  tree->Branch("y0In0",   &event.y0In0,  "y0In0/D");
  tree->Branch("u0In0",   &event.u0In0,  "u0In0/D");
  tree->Branch("v0In0",   &event.v0In0,  "v0In0/D");
  tree->Branch("theta0",  &event.theta0, "theta0/D");
  tree->Branch("p0",   &event.p0,  "p0/D");
  tree->Branch("x0In",   &event.x0In,  "x0In/D");
  tree->Branch("y0In",   &event.y0In,  "y0In/D");
  tree->Branch("u0In",   &event.u0In,  "u0In/D");
  tree->Branch("v0In",   &event.v0In,  "v0In/D");
  tree->Branch("thetaIn",  &event.thetaIn, "thetaIn/D");
  tree->Branch("p",   &event.p,  "p/D");
  tree->Branch("x0",   &event.x0,  "x0/D");
  tree->Branch("y0",   &event.y0,  "y0/D");
  tree->Branch("u0",   &event.u0,  "u0/D");
  tree->Branch("v0",   &event.v0,  "v0/D");
  tree->Branch("thetat",  &event.thetat, "thetat/D");
  tree->Branch("chisqr",  &event.chisqr, "chisqr/D");
  tree->Branch("mm1",   &event.mm1,  "mm1/D");
  tree->Branch("mm2",   &event.mm2,  "mm2/D");
  tree->Branch("bind",   &event.bind,  "bind/D");

  //Primary for Yd scattering
  //tree->Branch("thetaMeson",   &event.thetaMeson,  "thetaMeson/D");

}
