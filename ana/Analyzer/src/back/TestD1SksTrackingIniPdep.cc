/*
  TestD1SksTrackingIniPdep.cc

  2004/7/12   T.Takahashi

*/

#include "FieldMan.hh"
#include "DCLocalTrack.hh"
#include "DCHit.hh"
#include "DCConfMan.hh"
#include "SksTrack.hh"
#include "TrackHit.hh"
#include "ThreeVector.hh"
#include "DMag1.hh"
 
#include "RootHelper.hh"

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <string>

#define Interactive 1

const double Deg2Rad = acos(-1.)/180.;
const double Rad2Deg = 180./acos(-1.);

bool ProcessOneEvent( std::ifstream & );
void DefineHistograms( const char * );

const double MaxChiSqr = 20.;

int main( int argc, char **argv )
{
  if( argc<5 ) {
    std::cerr << "Usage: " << argv[0] 
	      << " confFile fieldMap dataFile RootFile"
	      << std::endl;
    exit(-1);
  }

  DCConfMan::GetInstance().Initialize(argv[1]);
  FieldMan::GetInstance().Initialize(argv[2]);

  // D1

  //  DMag1 D1( "D1", ThreeVector( 1698.30, 44.53, 0.00 ),
  //	    0.0, 90.0, -125.0, 150.0, -1.5 ); 
  DMag1 D1( "D1", ThreeVector( 1732.50, 138.50, 0.00 ),
	    0.0, 90.0, -125.0, 150.0, -1.5 ); 
 
  FieldMan::GetInstance().AddElement( &D1 );

  std::ifstream InputData(argv[3]);

  DefineHistograms(argv[4]);

  while( ProcessOneEvent( InputData ) ){
#ifdef Interactive 
    int it; 
    std::cout << "# "; 
    std::cin >> it;
    if( it<0 ) break;
#endif
   }

  gFile->Write();
  gFile->Close();

  return 0;
}

const int NBin1HPDcIn[11] = {
  100, 250, 250, 250, 250, 500, 500, 700, 700, 700, 700
};

const int NBin2HPDcIn[11] = {
  100, 50, 50, 50, 50, 100, 100, 140, 140, 140, 140
};

const double MinHPDcIn[11] = {
  0.0, 
  // DC1
  -125.0, -125.0, -125.0, -125.0,
  // DC2 DC3
  -250.0, -250.0, -350.0, -350.0, -350.0, -350.0
};

const double MaxHPDcIn[11] = {
  0.0,
  // DC1
  125.0, 125.0, 125.0, 125.0,
  // DC2
  250.0, 250.0, 350.0, 350.0, 350.0, 350.0
};

const int NBin1HPDcOut[13] = {
  100, 450, 450, 450, 450, 500, 500, 500, 500, 500, 500, 500, 500
};

const int NBin2HPDcOut[13] = {
  100, 90, 90, 90, 90, 100, 100, 100, 100, 100, 100, 100, 100
};

const double MinHPDcOut[13] = {
  0.0,
  // DC4
  -450.0, -450.0, -450.0, -450.0, 
  // DC5 DC6
  -500.0, -500.0, -500.0, -500.0, -500.0, -500.0, -500.0, -500.0
};

const double MaxHPDcOut[13] = {
  0.0,
  // DC4
  450.0, 450.0, 450.0, 450.0, 
  // DC5 DC6
  500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0
};



void DefineHistograms( const char *filename )
{
  new TFile( filename, "recreate" );

  HB1( 101, "stat (1)",  20, 0., 20. );
  HB1( 102, "status (1)",  10, 0., 10. );
  HB1( 103, "chisqr (1)", 500, 0., 100. );
  HB1( 104, "#Hits (1)",   30, 0., 30. );
  HB1( 105, "X0 (1)", 200, -50., 50. );
  HB1( 106, "Y0 (1)", 200, -50., 50. );
  HB1( 107, "P [GeV/c] (1)", 1000, 0.5, 1.5 );
  HB1( 108, "Theta [degree] (1)", 100, 0., 25. );
  HB1( 109, "Phi [degree] (1)", 100, -200., 200. );
  HB1( 110, "PathLength [mm] (1)", 600, 4500., 7500. );
  HB1( 111, "dX [mm] (1)", 200, -10., 10. ); 
  HB1( 112, "dY [mm] (1)", 200, -10., 10. ); 
  HB1( 113, "dP [GeV/c] (1)", 500, -0.025, 0.025 );
  HB1( 114, "dTheta [deg] (1)", 100, -5., 5. );
  HB1( 115, "dPhi [deg] (1)", 200, -50., 50. );
  HB1( 116, "dU (1)", 500, -0.05, 0.05 );
  HB1( 117, "DV (1)", 500, -0.05, 0.05 );
  HB2( 118, "dP%P [GeV/c] (1)", 200, 0.5, 1.5, 100, -0.01, 0.01 );

  HB1( 123, "chisqr <Success> (1)", 500, 0., 100. );
  HB1( 124, "#Hits <Success> (1)",   30, 0., 30. );
  HB1( 125, "X0 <Success> (1)", 200, -50., 50. );
  HB1( 126, "Y0 <Success> (1)", 200, -50., 50. );
  HB1( 127, "P [GeV/c] <Success> (1)", 1000, 0.5, 1.5 );
  HB1( 128, "Theta [degree] <Success> (1)", 100, 0., 25. );
  HB1( 129, "Phi [degree] <Success> (1)", 100, -200., 200. );
  HB1( 130, "PathLength [mm] <Success> (1)", 600, 4500., 7500. );
  HB1( 131, "dX [mm] <Success> (1)", 200, -10., 10. ); 
  HB1( 132, "dY [mm] <Success> (1)", 200, -10., 10. ); 
  HB1( 133, "dP [GeV/c] <Success> (1)", 500, -0.025, 0.025 );
  HB1( 134, "dTheta [deg] <Success> (1)", 100, -5., 5. );
  HB1( 135, "dPhi [deg] <Success> (1)", 200, -50., 50. );
  HB1( 136, "dU <Success> (1)", 500, -0.05, 0.05 );
  HB1( 137, "DV <Success> (1)", 500, -0.05, 0.05 );
  HB2( 138, "dP%P [GeV/c] <Success> (1)", 200, 0.5, 1.5, 100, -0.01, 0.01 );

  HB1( 143, "chisqr <Fail> (1)", 500, 0., 100. );
  HB1( 144, "#Hits <Fail> (1)",   30, 0., 30. );
  HB1( 145, "X0 <Fail> (1)", 200, -50., 50. );
  HB1( 146, "Y0 <Fail> (1)", 200, -50., 50. );
  HB1( 147, "P [GeV/c] <Fail> (1)", 1000, 0.5, 1.5 );
  HB1( 148, "Theta [degree] <Fail> (1)", 100, 0., 25. );
  HB1( 149, "Phi [degree] <Fail> (1)", 100, -200., 200. );
  HB1( 150, "PathLength [mm] <Fail> (1)", 600, 4500., 7500. );
  HB1( 151, "dX [mm] <Fail> (1)", 200, -10., 10. ); 
  HB1( 152, "dY [mm] <Fail> (1)", 200, -10., 10. ); 
  HB1( 153, "dP [GeV/c] <Fail> (1)", 500, -0.025, 0.025 );
  HB1( 154, "dTheta [deg] <Fail> (1)", 100, -5., 5. );
  HB1( 155, "dPhi [deg] <Fail> (1)", 200, -50., 50. );
  HB1( 156, "dU <Fail> (1)", 500, -0.05, 0.05 );
  HB1( 157, "DV <Fail> (1)", 500, -0.05, 0.05 );
  HB2( 158, "dP%P [GeV/c] <Fail> (1)", 200, 0.5, 1.5, 100, -0.01, 0.01 );

  HB1( 164, "#Hits <GoodTrack> (1)",   30, 0., 30. );
  HB1( 165, "X0 <GoodTrack> (1)", 200, -50., 50. );
  HB1( 166, "Y0 <GoodTrack> (1)", 200, -50., 50. );
  HB1( 167, "P [GeV/c] <GoodTrack> (1)", 1000, 0.5, 1.5 );
  HB1( 168, "Theta [degree] <GoodTrack> (1)", 100, 0., 25. );
  HB1( 169, "Phi [degree] <GoodTrack> (1)", 100, -200., 200. );
  HB1( 170, "PathLength [mm] <GoodTrack> (1)", 600, 4500., 7500. );
  HB1( 171, "dX [mm] <GoodTrack> (1)", 200, -10., 10. ); 
  HB1( 172, "dY [mm] <GoodTrack> (1)", 200, -10., 10. ); 
  HB1( 173, "dP [GeV/c] <GoodTrack> (1)", 500, -0.025, 0.025 );
  HB1( 174, "dTheta [deg] <GoodTrack> (1)", 100, -5., 5. );
  HB1( 175, "dPhi [deg] <GoodTrack> (1)", 200, -50., 50. );
  HB1( 176, "dU <GoodTrack> (1)", 500, -0.05, 0.05 );
  HB1( 177, "DV <GoodTrack> (1)", 500, -0.05, 0.05 );
  HB2( 178, "dP%P [GeV/c] <GoodTrack> (1)", 200, 0.5, 1.5, 100, -0.01, 0.01 );

  HB1( 201, "stat (2)",  20, 0., 20. );
  HB1( 202, "status (2)",  10, 0., 10. );
  HB1( 203, "chisqr (2)", 500, 0., 100. );
  HB1( 204, "#Hits (2)",   30, 0., 30. );
  HB1( 205, "X0 (2)", 200, -50., 50. );
  HB1( 206, "Y0 (2)", 200, -50., 50. );
  HB1( 207, "P [GeV/c] (2)", 1000, 0.5, 1.5 );
  HB1( 208, "Theta [degree] (2)", 100, 0., 25. );
  HB1( 209, "Phi [degree] (2)", 100, -200., 200. );
  HB1( 210, "PathLength [mm] (2)", 600, 4500., 7500. );
  HB1( 211, "dX [mm] (2)", 200, -10., 10. ); 
  HB1( 212, "dY [mm] (2)", 200, -10., 10. ); 
  HB1( 213, "dP [GeV/c] (2)", 500, -0.025, 0.025 );
  HB1( 214, "dTheta [deg] (2)", 100, -5., 5. );
  HB1( 215, "dPhi [deg] (2)", 200, -50., 50. );
  HB1( 216, "dU (2)", 500, -0.05, 0.05 );
  HB1( 217, "DV (2)", 500, -0.05, 0.05 );
  HB2( 218, "dP%P [GeV/c] (2)", 200, 0.5, 1.5, 100, -0.01, 0.01 );

  HB1( 223, "chisqr <Success> (2)", 500, 0., 100. );
  HB1( 224, "#Hits <Success> (2)",   30, 0., 30. );
  HB1( 225, "X0 <Success> (2)", 200, -50., 50. );
  HB1( 226, "Y0 <Success> (2)", 200, -50., 50. );
  HB1( 227, "P [GeV/c] <Success> (2)", 1000, 0.5, 1.5 );
  HB1( 228, "Theta [degree] <Success> (2)", 100, 0., 25. );
  HB1( 229, "Phi [degree] <Success> (2)", 100, -200., 200. );
  HB1( 230, "PathLength [mm] <Success> (2)", 600, 4500., 7500. );
  HB1( 231, "dX [mm] <Success> (2)", 200, -10., 10. ); 
  HB1( 232, "dY [mm] <Success> (2)", 200, -10., 10. ); 
  HB1( 233, "dP [GeV/c] <Success> (2)", 500, -0.025, 0.025 );
  HB1( 234, "dTheta [deg] <Success> (2)", 100, -5., 5. );
  HB1( 235, "dPhi [deg] <Success> (2)", 200, -50., 50. );
  HB1( 236, "dU <Success> (2)", 500, -0.05, 0.05 );
  HB1( 237, "DV <Success> (2)", 500, -0.05, 0.05 );
  HB2( 238, "dP%P [GeV/c] <Success> (2)", 200, 0.5, 1.5, 100, -0.01, 0.01 );

  HB1( 243, "chisqr <Fail> (2)", 500, 0., 100. );
  HB1( 244, "#Hits <Fail> (2)",   30, 0., 30. );
  HB1( 245, "X0 <Fail> (2)", 200, -50., 50. );
  HB1( 246, "Y0 <Fail> (2)", 200, -50., 50. );
  HB1( 247, "P [GeV/c] <Fail> (2)", 1000, 0.5, 1.5 );
  HB1( 248, "Theta [degree] <Fail> (2)", 100, 0., 25. );
  HB1( 249, "Phi [degree] <Fail> (2)", 100, -200., 200. );
  HB1( 250, "PathLength [mm] <Fail> (2)", 600, 4500., 7500. );
  HB1( 251, "dX [mm] <Fail> (2)", 200, -10., 10. ); 
  HB1( 252, "dY [mm] <Fail> (2)", 200, -10., 10. ); 
  HB1( 253, "dP [GeV/c] <Fail> (2)", 500, -0.025, 0.025 );
  HB1( 254, "dTheta [deg] <Fail> (2)", 100, -5., 5. );
  HB1( 255, "dPhi [deg] <Fail> (2)", 200, -50., 50. );
  HB1( 256, "dU <Fail> (2)", 500, -0.05, 0.05 );
  HB1( 257, "DV <Fail> (2)", 500, -0.05, 0.05 );
  HB2( 258, "dP%P [GeV/c] <Fail> (2)", 200, 0.5, 1.5, 100, -0.01, 0.01 );

  HB1( 264, "#Hits <GoodTrack> (2)",   30, 0., 30. );
  HB1( 265, "X0 <GoodTrack> (2)", 200, -50., 50. );
  HB1( 266, "Y0 <GoodTrack> (2)", 200, -50., 50. );
  HB1( 267, "P [GeV/c] <GoodTrack> (2)", 1000, 0.5, 1.5 );
  HB1( 268, "Theta [degree] <GoodTrack> (2)", 100, 0., 25. );
  HB1( 269, "Phi [degree] <GoodTrack> (2)", 100, -200., 200. );
  HB1( 270, "PathLength [mm] <GoodTrack> (2)", 600, 4500., 7500. );
  HB1( 271, "dX [mm] <GoodTrack> (2)", 200, -10., 10. ); 
  HB1( 272, "dY [mm] <GoodTrack> (2)", 200, -10., 10. ); 
  HB1( 273, "dP [GeV/c] <GoodTrack> (2)", 500, -0.025, 0.025 );
  HB1( 274, "dTheta [deg] <GoodTrack> (2)", 100, -5., 5. );
  HB1( 275, "dPhi [deg] <GoodTrack> (2)", 200, -50., 50. );
  HB1( 276, "dU <GoodTrack> (2)", 500, -0.05, 0.05 );
  HB1( 277, "DV <GoodTrack> (2)", 500, -0.05, 0.05 );
  HB2( 278, "dP%P [GeV/c] <GoodTrack> (2)", 200, 0.5, 1.5, 100, -0.01, 0.01 );


  HB2( 301, "ChiSqr(2)%ChiSqr(1)", 100, 0., 50., 100, 0., 50. );
  HB2( 302, "dX(2)%dX(1)", 100, -5., 5., 100, -5., 5. );
  HB2( 303, "dY(2)%dY(1)", 100, -5., 5., 100, -5., 5. );
  HB2( 304, "dP(2)%dP(1)", 100, -0.02, 0.02, 100, -0.02, 0.02 );
  HB2( 305, "dTheta(2)%dTheta(1)", 100, -5., 5., 100, -5., 5. );
  HB2( 306, "dPhi(2)%dPhi(1)", 100, -25., 25., 100, -25., 25. );
  HB2( 307, "dU(2)%dU(1)", 100, -0.05, 0.05, 100, -0.05, 0.05 );
  HB2( 308, "dV(2)%dV(1)", 100, -0.05, 0.05, 100, -0.05, 0.05 );

  HB2( 321, "ChiSqr(2)%ChiSqr(1) <Success>", 100, 0., 50., 100, 0., 50. );
  HB2( 322, "dX(2)%dX(1) <Success>", 100, -5., 5., 100, -5., 5. );
  HB2( 323, "dY(2)%dY(1) <Success>", 100, -5., 5., 100, -5., 5. );
  HB2( 324, "dP(2)%dP(1) <Success>", 100, -0.02, 0.02, 100, -0.02, 0.02 );
  HB2( 325, "dTheta(2)%dTheta(1) <Success>", 100, -5., 5., 100, -5., 5. );
  HB2( 326, "dPhi(2)%dPhi(1) <Success>", 100, -25., 25., 100, -25., 25. );
  HB2( 327, "dU(2)%dU(1) <Success>", 100, -0.05, 0.05, 100, -0.05, 0.05 );
  HB2( 328, "dV(2)%dV(1) <Success>", 100, -0.05, 0.05, 100, -0.05, 0.05 );

  HB2( 341, "ChiSqr(2)%ChiSqr(1) <Fail>", 100, 0., 50., 100, 0., 50. );
  HB2( 342, "dX(2)%dX(1) <Fail>", 100, -5., 5., 100, -5., 5. );
  HB2( 343, "dY(2)%dY(1) <Fail>", 100, -5., 5., 100, -5., 5. );
  HB2( 344, "dP(2)%dP(1) <Fail>", 100, -0.02, 0.02, 100, -0.02, 0.02 );
  HB2( 345, "dTheta(2)%dTheta(1) <Fail>", 100, -5., 5., 100, -5., 5. );
  HB2( 346, "dPhi(2)%dPhi(1) <Fail>", 100, -25., 25., 100, -25., 25. );
  HB2( 347, "dU(2)%dU(1) <Fail>", 100, -0.05, 0.05, 100, -0.05, 0.05 );
  HB2( 348, "dV(2)%dV(1) <Fail>", 100, -0.05, 0.05, 100, -0.05, 0.05 );

  HB2( 361, "ChiSqr(2)%ChiSqr(1) <GoodTrack>", 100, 0., 50., 100, 0., 50. );
  HB2( 362, "dX(2)%dX(1) <GoodTrack>", 100, -5., 5., 100, -5., 5. );
  HB2( 363, "dY(2)%dY(1) <GoodTrack>", 100, -5., 5., 100, -5., 5. );
  HB2( 364, "dP(2)%dP(1) <GoodTrack>", 100, -0.02, 0.02, 100, -0.02, 0.02 );
  HB2( 365, "dTheta(2)%dTheta(1) <GoodTrack>", 100, -5., 5., 100, -5., 5. );
  HB2( 366, "dPhi(2)%dPhi(1) <GoodTrack>", 100, -25., 25., 100, -25., 25. );
  HB2( 367, "dU(2)%dU(1) <GoodTrack>", 100, -0.05, 0.05, 100, -0.05, 0.05 );
  HB2( 368, "dV(2)%dV(1) <GoodTrack>", 100, -0.05, 0.05, 100, -0.05, 0.05 );

}

bool ProcessOneEvent( std::ifstream &In )
{
  const std::string funcname = "ProcessEvent";

  if( !In ) return false;
  double gx,gy,gz,gp,gt,gf;
  In >> gx >> gy >> gz >> gp >> gt >> gf;
  DCLocalTrack trIn, trOut;
  int lnum;
  In >> lnum;
  volatile int lnumold=lnum;

  while( lnum!=-1 ){
    //    std::cout << lnumold << " "; 
    double pos=0.;
    In >> pos >> lnum;
    DCHit *hit = new DCHit( lnumold, 0, 0 );
    if( hit ){
      hit->SetLocalHitPos( pos );
      hit->SetTiltAngle( DCConfMan::GetInstance().GetTiltAngle(lnumold) );
      if( lnumold>0 && lnumold<=30 )
	trIn.AddHit( hit );
      else if( lnumold>40 && lnumold<=52 )
	trOut.AddHit( hit );
      else
	delete hit;
    }
    lnumold = lnum;
  }
  //  std::cout << std::endl;

  HF1( 101, 0. );
  HF1( 201, 0. );

  //  trIn.DoFit(); 
  trOut.DoFit();

  //  bool stIn=trIn.GetStatus();
  bool stOut=trOut.GetStatus();

  std::size_t nhin=trIn.GetNHit(), nhout=trOut.GetNHit();

  //  if( !stIn || !stOut || nhin!=10 || nhout !=12 ){
  if( !stOut ){
    for( std::size_t i=0; i<nhin; ++i )
      delete trIn.GetHit(i);
    for( std::size_t i=0; i<nhout; ++i )
      delete trOut.GetHit(i);

    return true; 
  }


  SksTrack track1( &trIn, &trOut );
  SksTrack track2( &trIn, &trOut );
  //  track.SetInitialMomentum(gp);
  track1.SetInitialMomentum( 1.2 );
  track2.SetInitialMomentum( 1.3 );
  track1.doFit();
  track2.doFit();
  bool status1=track1.Status();
  bool status2=track2.Status();
  double chisqr1=track1.chisqr();
  double chisqr2=track2.chisqr();

  ThreeVector priPos1 = track1.PrimaryPosition();
  ThreeVector priPos2 = track2.PrimaryPosition();
  ThreeVector priMom1 = track1.PrimaryMomentum();
  ThreeVector priMom2 = track2.PrimaryMomentum();

  double rx1=priPos1.x(), ry1=priPos1.y(), rz1=priPos1.z();
  double rp1=priMom1.mag();
  double rt1=priMom1.theta()*Rad2Deg, rf1=priMom1.phi()*Rad2Deg;
  double rx2=priPos2.x(), ry2=priPos2.y(), rz2=priPos2.z();
  double rp2=priMom2.mag();
  double rt2=priMom2.theta()*Rad2Deg, rf2=priMom2.phi()*Rad2Deg;
  double lengTof1 = track1.PathLengthToTOF();
  double lengTof2 = track2.PathLengthToTOF();

  int nhdc1 = track1.GetNHits();
  int nhdc2 = track2.GetNHits();

#ifdef Interactive
  {
    std::ios::fmtflags oldFlags = std::cout.flags();
    std::size_t oldPrec = std::cout.precision();
    std::cout.setf( std::ios::fixed );
    std::cout.precision(4);

    std::cout << "[TestQSksTracking]: Track1 St=" << status1
	      << " chi2=" << std::setw(10) << chisqr1 
	      << " PL=" << std::setprecision(1) 
	      << std::setw(7) << lengTof1  << " [mm]" << std::endl;
    std::cout.precision(2);
    std::cout << " Position: ( "
	      << std::setw(7) << rx1 << ", "
	      << std::setw(7) << ry1 << " )" << " - ( "
	      << std::setw(7) << gx << ", "
	      << std::setw(7) << gy << " )" << " =====> ( "
	      << std::setw(7) << rx1-gx << ", "
	      << std::setw(7) << ry1-gy << " ) [mm]" << std::endl;
    std::cout.precision(5);
    std::cout << " P: " << std::setw(8) << rp1 << " ( "
	      << std::setw(8) << gp << " ) < "
	      << std::setw(8) << rp1-gp << " > [GeV/c]";
    std::cout.precision(1);
    std::cout << "  Angles: ( " 
	      << std::setw(5) << rt1 << ", " 
	      << std::setw(6) << rf1 << " ) <( "
	      << std::setw(5) << gt << ", " 
	      << std::setw(5) << gf << " )> [degree]"
	      << std::endl; 

    std::cout << "[TestQSksTracking]: Track2 St=" << status2
	      << " chi2=" << std::setw(10) << chisqr2 
	      << " PL=" << std::setprecision(1) 
	      << std::setw(7) << lengTof2  << " [mm]" << std::endl;
    std::cout.precision(2);
    std::cout << " Position: ( "
	      << std::setw(7) << rx2 << ", "
	      << std::setw(7) << ry2 << " )" << " - ( "
	      << std::setw(7) << gx << ", "
	      << std::setw(7) << gy << " )" << " =====> ( "
	      << std::setw(7) << rx2-gx << ", "
	      << std::setw(7) << ry2-gy << " ) [mm]" << std::endl;
    std::cout.precision(5);
    std::cout << " P: " << std::setw(8) << rp2 << " ( "
	      << std::setw(8) << gp << " ) < "
	      << std::setw(8) << rp2-gp << " > [GeV/c]";
    std::cout.precision(1);
    std::cout << "  Angles: ( " 
	      << std::setw(5) << rt2 << ", " 
	      << std::setw(6) << rf2 << " ) <( "
	      << std::setw(5) << gt << ", " 
	      << std::setw(5) << gf << " )> [degree]"
	      << std::endl; 
      
    std::cout.flags( oldFlags );
    std::cout.precision( oldPrec );
  }
#endif
  double gu = tan(gt*Deg2Rad)*cos(gf*Deg2Rad);
  double gv = tan(gt*Deg2Rad)*sin(gf*Deg2Rad);
  double ru1 = tan(rt1*Deg2Rad)*cos(rf1*Deg2Rad);
  double rv1 = tan(rt1*Deg2Rad)*sin(rf1*Deg2Rad);
  double ru2 = tan(rt2*Deg2Rad)*cos(rf2*Deg2Rad);
  double rv2 = tan(rt2*Deg2Rad)*sin(rf2*Deg2Rad);

  HF1( 101, 1. ); 
  HF1( 102, status1 );
  HF1( 103, chisqr1 );
  HF1( 104, nhdc1 );
  HF1( 105, rx1 ); HF1( 106, ry1 ); HF1( 107, rp1 );
  HF1( 108, rt1 ); HF1( 109, rf1 );
  HF1( 110, lengTof1 );
  HF1( 111, rx1-gx ); HF1( 112, ry1-gy ); HF1( 113, rp1-gp );
  HF1( 114, rt1-gt ); HF1( 115, rf1-gf );
  HF1( 116, ru1-gu ); HF1( 117, rv1-gv );
  HF2( 118, gp, rp1-gp );
 
  HF1( 201, 1. );
  HF1( 202, status2 );
  HF1( 203, chisqr2 );
  HF1( 204, nhdc2 );
  HF1( 205, rx2 ); HF1( 206, ry2 ); HF1( 207, rp2 );
  HF1( 208, rt2 ); HF1( 209, rf2 );
  HF1( 210, lengTof2 );
  HF1( 211, rx2-gx ); HF1( 212, ry2-gy ); HF1( 213, rp2-gp );
  HF1( 214, rt2-gt ); HF1( 215, rf2-gf );
  HF1( 216, ru2-gu ); HF1( 217, rv2-gv );
  HF2( 218, gp, rp2-gp );

  HF2( 301, chisqr1, chisqr2 );
  HF2( 302, rx1-gx, rx2-gx );
  HF2( 303, ry1-gy, ry2-gy );
  HF2( 304, rp1-gp, rp2-gp );
  HF2( 305, rt1-gt, rt2-gt );
  HF2( 306, rf1-gf, rf2-gf );
  HF2( 307, ru1-gu, ru2-gu );
  HF2( 308, rv1-gv, rv2-gv );

  if( status1 ){
    HF1( 101, 2. ); 
    HF1( 123, chisqr1 );
    HF1( 124, nhdc1 );
    HF1( 125, rx1 ); HF1( 126, ry1 ); HF1( 127, rp1 );
    HF1( 128, rt1 ); HF1( 129, rf1 );
    HF1( 130, lengTof1 );
    HF1( 131, rx1-gx ); HF1( 132, ry1-gy ); HF1( 133, rp1-gp );
    HF1( 134, rt1-gt ); HF1( 135, rf1-gf );
    HF1( 136, ru1-gu ); HF1( 137, rv1-gv );
    HF2( 138, gp, rp1-gp );

    if( chisqr1<MaxChiSqr ){
      HF1( 101, 3. ); 
      HF1( 164, nhdc1 );
      HF1( 165, rx1 ); HF1( 166, ry1 ); HF1( 167, rp1 );
      HF1( 168, rt1 ); HF1( 169, rf1 );
      HF1( 170, lengTof1 );
      HF1( 171, rx1-gx ); HF1( 172, ry1-gy ); HF1( 173, rp1-gp );
      HF1( 174, rt1-gt ); HF1( 175, rf1-gf );
      HF1( 176, ru1-gu ); HF1( 177, rv1-gv );
      HF2( 178, gp, rp1-gp );
    }
  }
  else {
    HF1( 144, nhdc1 );
    HF1( 145, rx1 ); HF1( 146, ry1 ); HF1( 147, rp1 );
    HF1( 148, rt1 ); HF1( 149, rf1 );
    HF1( 150, lengTof1 );
    HF1( 151, rx1-gx ); HF1( 152, ry1-gy ); HF1( 153, rp1-gp );
    HF1( 154, rt1-gt ); HF1( 155, rf1-gf );
    HF1( 156, ru1-gu ); HF1( 157, rv1-gv );
    HF2( 158, gp, rp1-gp );
  }

  if( status2 ){
    HF1( 201, 2. );
    HF1( 223, chisqr2 );
    HF1( 224, nhdc2 );
    HF1( 225, rx2 ); HF1( 226, ry2 ); HF1( 227, rp2 );
    HF1( 228, rt2 ); HF1( 229, rf2 );
    HF1( 230, lengTof2 );
    HF1( 231, rx2-gx ); HF1( 232, ry2-gy ); HF1( 233, rp2-gp );
    HF1( 234, rt2-gt ); HF1( 235, rf2-gf );
    HF1( 236, ru2-gu ); HF1( 237, rv2-gv );
    HF2( 238, gp, rp2-gp );

    if( chisqr2<MaxChiSqr ){
      HF1( 201, 3. );
      HF1( 264, nhdc2 );
      HF1( 265, rx2 ); HF1( 266, ry2 ); HF1( 267, rp2 );
      HF1( 268, rt2 ); HF1( 269, rf2 );
      HF1( 270, lengTof2 );
      HF1( 271, rx2-gx ); HF1( 272, ry2-gy ); HF1( 273, rp2-gp );
      HF1( 274, rt2-gt ); HF1( 275, rf2-gf );
      HF1( 276, ru2-gu ); HF1( 277, rv2-gv );
      HF2( 278, gp, rp2-gp );
    }
  }
  else {
    HF1( 244, nhdc2 );
    HF1( 245, rx2 ); HF1( 246, ry2 ); HF1( 247, rp2 );
    HF1( 248, rt2 ); HF1( 249, rf2 );
    HF1( 250, lengTof2 );
    HF1( 251, rx2-gx ); HF1( 252, ry2-gy ); HF1( 253, rp2-gp );
    HF1( 254, rt2-gt ); HF1( 255, rf2-gf );
    HF1( 256, ru2-gu ); HF1( 257, rv2-gv );
    HF2( 258, gp, rp2-gp );
  }

  if( status1 && status2 ){
    HF2( 321, chisqr1, chisqr2 );
    HF2( 322, rx1-gx, rx2-gx );
    HF2( 323, ry1-gy, ry2-gy );
    HF2( 324, rp1-gp, rp2-gp );
    HF2( 325, rt1-gt, rt2-gt );
    HF2( 326, rf1-gf, rf2-gf );
    HF2( 327, ru1-gu, ru2-gu );
    HF2( 328, rv1-gv, rv2-gv );
    if( chisqr1<MaxChiSqr && chisqr2<MaxChiSqr ){
      HF2( 361, chisqr1, chisqr2 );
      HF2( 362, rx1-gx, rx2-gx );
      HF2( 363, ry1-gy, ry2-gy );
      HF2( 364, rp1-gp, rp2-gp );
      HF2( 365, rt1-gt, rt2-gt );
      HF2( 366, rf1-gf, rf2-gf );
      HF2( 367, ru1-gu, ru2-gu );
      HF2( 368, rv1-gv, rv2-gv );
    }
  }
  else {
    HF2( 341, chisqr1, chisqr2 );
    HF2( 342, rx1-gx, rx2-gx );
    HF2( 343, ry1-gy, ry2-gy );
    HF2( 344, rp1-gp, rp2-gp );
    HF2( 345, rt1-gt, rt2-gt );
    HF2( 346, rf1-gf, rf2-gf );
    HF2( 347, ru1-gu, ru2-gu );
    HF2( 348, rv1-gv, rv2-gv );
  }


  for( std::size_t i=0; i<nhin; ++i )
    delete trIn.GetHit(i);
  for( std::size_t i=0; i<nhout; ++i )
    delete trOut.GetHit(i);

  return true;
}




















