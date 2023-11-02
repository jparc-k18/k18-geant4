/*
  LocalTrackingIn.cc
*/

#include "DCLocalTrack.hh"
#include "DCHit.hh"
#include "DCConfMan.hh"

#include "RootHelper.hh"

#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>

//#define Interactive 1

const double Deg2Rad = acos(-1.)/180.;
const double Rad2Deg = 180./acos(-1.);

bool ProcessOneEvent( std::istream & );
void DefineHistograms( const char * );


int main( int argc, char **argv )
{
  if( argc<4 ){
    std::cerr << "Usage: " << argv[0] << " confFile dataFile RootFile" 
	      << std::endl;
    exit(-1);
  }

  DCConfMan::GetInstance().Initialize(argv[1]);

  std::ifstream InputData(argv[2]);

  DefineHistograms(argv[3]);

  while( ProcessOneEvent( InputData ) ){
#ifdef Interactive 
    int it; 
    std::cout << "# "; 
    std::cin >> it;
#endif
  }

  gFile->Write();
  gFile->Close();

  return 0;
}

void DefineHistograms( const char *filename )
{
  new TFile( filename, "recreate" );

  HB1( 101, "status",    10, 0., 10. );
  HB1( 102, "chisquare", 100, 0., 10. );
  HB1( 103, "#Hits",     20, 0., 20. );
  HB1( 104, "X0[mm]",   2000, -50., 50. ); 
  HB1( 105, "Y0[mm]",   2000, -50., 50. ); 
  HB1( 106, "U0", 4000, -0.4, 0.4 );
  HB1( 107, "V0", 3000, -0.3, 0.3 );
  HB1( 108, "Theta[degree]", 1000, 0., 25. );
  HB1( 109, "Phi[degree]", 1000, -200., 200. );
  HB1( 110, "LayerDist", 20, 0., 20. );

  HB2( 111, "U0%X0", 200, -25., 25., 200, -0.4, 0.4 );
  HB2( 112, "V0%Y0", 200, -25., 25., 200, -0.25, 0.25 );
  HB2( 113, "Y0%x0", 200, -25., 25., 200, -25., 25. ); 
  HB2( 114, "Phi%Theta", 100, 0., 25., 100, -200., 200. );

  HB1( 115, "dX [mm]", 200, -10., 10. ); 
  HB1( 116, "dY [mm]", 200, -10., 10. ); 
  HB1( 117, "dU", 5000, -0.05, 0.05 );
  HB1( 118, "dV", 5000, -0.05, 0.05 );

  for( int i=1; i<=10; ++i ){
    std::ostringstream title1, title2, title3;
    title1 << "HitPos[mm] Layer#" << i;
    HB1( 120+i, title1.str().c_str(), 500, -250., 250. );
    title2 << "Residual[mm] Layer#" << i;
    HB1( 130+i, title2.str().c_str(), 100, -2.5, 2.5 );
    title3 << "Resd.%Pos[mm] Layer#" << i;
    HB2( 140+i, title3.str().c_str(), 100, -250., 250., 100, -2.5, 2.5 );
  }
}


bool ProcessOneEvent( std::istream &In )
{
  const std::string funcname = "ProcessEvent";

  if( !In ) return false;
  double gx,gy,gz,gp,gt,gf;
  In >> gx >> gy >> gz >> gp >> gt >> gf;

  DCLocalTrack track;
  while( In ){
    int lnum;
    double pos;
    In >> lnum;
    if( lnum<0 ) break;
    In >> pos;
    DCHit *hit = new DCHit( lnum, 0, 0 );
    if( hit ){
      hit->SetLocalHitPos( pos );
      hit->SetTiltAngle( DCConfMan::GetInstance().GetTiltAngle(lnum) );
      track.AddHit( hit );
    }
  }
  track.DoFit(); 

  bool status=track.GetStatus();
  double chisqr=track.GetChiSquare();
  int nhit=track.GetNHit();
  double x0=track.GetX0(), y0=track.GetY0();
  double u0=track.GetU0(), v0=track.GetV0();
  double theta=acos(1./sqrt(1.+u0*u0+v0*v0))*Rad2Deg;
  double phi=atan2(v0,u0)*Rad2Deg;

  HF1( 101, status ); 
  HF1( 103, double(nhit) );
  if( status ){
    HF1( 102, chisqr );
    HF1( 104, x0 );
    HF1( 105, y0 );
    HF1( 106, u0 );
    HF1( 107, v0 );
    HF1( 108, theta );
    HF1( 109, phi );
    HF2( 111, x0, u0 );
    HF2( 112, y0, v0 );
    HF2( 113, x0, y0 );
    HF2( 114, theta, phi );

    double gu = tan(gt*Deg2Rad)*cos(gf*Deg2Rad);
    double gv = tan(gt*Deg2Rad)*sin(gf*Deg2Rad);

    HF1( 115, x0-gx );
    HF1( 116, y0-gy );
    HF1( 117, u0-gu );
    HF1( 118, v0-gv );

    for( std::size_t i=0; i<nhit; ++i ){
      DCHit *hit = track.GetHit(i);
      int layer = hit->GetLayer();
      HF1( 110, layer );
      HF1( 120+layer, hit->GetLocalHitPos() );
      HF1( 130+layer, hit->GetResidual() );
      std::cout << "Residual=" << hit->GetResidual() << std::endl;
      HF1( 140+layer, hit->GetLocalHitPos(), hit->GetResidual() );
    }
  }


  //#ifdef Interactive 
#if 1
  std::cout << "TrackStatus=" << status 
	    << " ChiSqr=" << chisqr
	    << " #Hits=" << nhit << std::endl;
  std::cout << "      " 
	    << " X0:" << x0 << " Y0:" << y0
	    << " U0:" << u0 << " V0:" << v0
	    << " Theta=" << theta << "[degree]"
	    << " Phi=" << phi << "[degree]"
	    << std::endl;
#endif


  for( std::size_t i=0; i<nhit; ++i ){
    delete track.GetHit(i);
  }

  return In;
}
