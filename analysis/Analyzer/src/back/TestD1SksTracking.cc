/*
  TestD1SksTracking.cc

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

//#define Interactive 1

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

  HB1( 101, "stat",  20, 0., 20. );
  HB1( 102, "status",  10, 0., 10. );
  HB1( 103, "chisqr", 500, 0., 100. );
  HB1( 104, "#Hits",   30, 0., 30. );
  HB1( 105, "X0", 200, -50., 50. );
  HB1( 106, "Y0", 200, -50., 50. );
  HB1( 107, "P [GeV/c]", 1000, 0.5, 1.5 );
  HB1( 108, "Theta [degree]", 100, 0., 25. );
  HB1( 109, "Phi [degree]", 100, -200., 200. );
  HB1( 110, "PathLength [mm]", 600, 4500., 7500. );
  HB1( 111, "dX [mm]", 200, -10., 10. ); 
  HB1( 112, "dY [mm]", 200, -10., 10. ); 
  HB1( 113, "dP [GeV/c]", 500, -0.025, 0.025 );
  HB1( 114, "dTheta [deg]", 100, -5., 5. );
  HB1( 115, "dPhi [deg]", 200, -50., 50. );
  HB1( 116, "dU", 500, -0.05, 0.05 );
  HB1( 117, "DV", 500, -0.05, 0.05 );
  HB2( 118, "dP%P [GeV/c]", 200, 0.5, 1.5, 100, -0.01, 0.01 );

  HB1( 123, "chisqr <Success>", 500, 0., 100. );
  HB1( 124, "#Hits <Success>",   30, 0., 30. );
  HB1( 125, "X0 <Success>", 200, -50., 50. );
  HB1( 126, "Y0 <Success>", 200, -50., 50. );
  HB1( 127, "P [GeV/c] <Success>", 1000, 0.5, 1.5 );
  HB1( 128, "Theta [degree] <Success>", 100, 0., 25. );
  HB1( 129, "Phi [degree] <Success>", 100, -200., 200. );
  HB1( 130, "PathLength [mm] <Success>", 600, 4500., 7500. );
  HB1( 131, "dX [mm] <Success>", 200, -10., 10. ); 
  HB1( 132, "dY [mm] <Success>", 200, -10., 10. ); 
  HB1( 133, "dP [GeV/c] <Success>", 500, -0.025, 0.025 );
  HB1( 134, "dTheta [deg] <Success>", 100, -5., 5. );
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
  HB1( 140, "PathLength [mm] <Fail>", 600, 4500., 7500. );
  HB1( 151, "dX [mm] <Fail>", 200, -10., 10. ); 
  HB1( 152, "dY [mm] <Fail>", 200, -10., 10. ); 
  HB1( 153, "dP [GeV/c] <Fail>", 500, -0.025, 0.025 );
  HB1( 154, "dTheta [deg] <Fail>", 100, -5., 5. );
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
  HB1( 170, "PathLength [mm] <GoodTrack>", 600, 4500., 7500. );
  HB1( 171, "dX [mm] <GoodTrack>", 200, -10., 10. ); 
  HB1( 172, "dY [mm] <GoodTrack>", 200, -10., 10. ); 
  HB1( 173, "dP [GeV/c] <GoodTrack>", 500, -0.025, 0.025 );
  HB1( 174, "dTheta [deg] <GoodTrack>", 100, -5., 5. );
  HB1( 175, "dPhi [deg] <GoodTrack>", 200, -50., 50. );
  HB1( 176, "dU <GoodTrack>", 500, -0.05, 0.05 );
  HB1( 177, "DV <GoodTrack>", 500, -0.05, 0.05 );
  HB2( 178, "dP%P [GeV/c] <GoodTrack>", 200, 0.5, 1.5, 100, -0.01, 0.01 );

  for( int i=1; i<=6; ++i ){
    std::ostringstream title1, title2, title3;
    title1 << "HitPos [mm] Layer#" << i;
    HB1( 200+i, title1.str().c_str(), NBin1HPDcIn[i], 
	 MinHPDcIn[i], MaxHPDcIn[i] );
    title2 << "Residual [mm] Layer#" << i;
    HB1( 300+i, title2.str().c_str(), 250, -2.5, 2.5 );
    title3 << "Resd.%Pos [mm] Layer#" << i;
    HB2( 400+i, title3.str().c_str(), NBin2HPDcIn[i], 
	 MinHPDcIn[i], MaxHPDcIn[i], 100, -2.5, 2.5 );

    title1 << " <GoodTrack>";
    HB1( 1200+i, title1.str().c_str(), NBin1HPDcIn[i], 
	 MinHPDcIn[i], MaxHPDcIn[i] );
    title2 << " <GoodTrack>";
    HB1( 1300+i, title2.str().c_str(), 250, -2.5, 2.5 );
    title3 << " <GoodTrack>";
    HB2( 1400+i, title3.str().c_str(), NBin2HPDcIn[i], 
	 MinHPDcIn[i], MaxHPDcIn[i], 100, -2.5, 2.5 );
  }

  for( int i=21; i<=24; ++i ){
    std::ostringstream title1, title2, title3;
    title1 << "HitPos [mm] Layer#" << i;
    HB1( 200+i, title1.str().c_str(), NBin1HPDcIn[i-14], 
	 MinHPDcIn[i-14], MaxHPDcIn[i-14] );
    title2 << "Residual [mm] Layer#" << i;
    HB1( 300+i, title2.str().c_str(), 250, -2.5, 2.5 );
    title3 << "Resd.%Pos [mm] Layer#" << i;
    HB2( 400+i, title3.str().c_str(), NBin2HPDcIn[i-14], 
	 MinHPDcIn[i-14], MaxHPDcIn[i-14], 100, -2.5, 2.5 );

    title1 << " <GoodTrack>";
    HB1( 1200+i, title1.str().c_str(), NBin1HPDcIn[i-14], 
	 MinHPDcIn[i-14], MaxHPDcIn[i-14] );
    title2 << " <GoodTrack>";
    HB1( 1300+i, title2.str().c_str(), 250, -2.5, 2.5 );
    title3 << " <GoodTrack>";
    HB2( 1400+i, title3.str().c_str(), NBin2HPDcIn[i-14], 
	 MinHPDcIn[i-14], MaxHPDcIn[i-14], 100, -2.5, 2.5 );
  }

  for( int i=41; i<=52; ++i ){
    std::ostringstream title1, title2, title3;
    title1 << "HitPos [mm] Layer#" << i;
    HB1( 200+i, title1.str().c_str(), NBin1HPDcOut[i-40], 
	 MinHPDcOut[i-40], MaxHPDcOut[i-40] );
    title2 << "Residual [mm] Layer#" << i;
    HB1( 300+i, title2.str().c_str(), 250, -2.5, 2.5 );
    title3 << "Resd.%Pos [mm] Layer#" << i;
    HB2( 400+i, title3.str().c_str(), NBin2HPDcOut[i-40], 
	 MinHPDcOut[i-40], MaxHPDcOut[i-40], 100, -2.5, 2.5 );

    title1 << " <GoodTrack>";
    HB1( 1200+i, title1.str().c_str(), NBin1HPDcOut[i-40], 
	 MinHPDcOut[i-40], MaxHPDcOut[i-40] );
    title2 << " <GoodTrack>";
    HB1( 1300+i, title2.str().c_str(), 250, -2.5, 2.5 );
    title3 << " <GoodTrack>";
    HB2( 1400+i, title3.str().c_str(), NBin2HPDcOut[i-40], 
	 MinHPDcOut[i-40], MaxHPDcOut[i-40], 100, -2.5, 2.5 );
  }

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


  SksTrack track( &trIn, &trOut );
  //  track.SetInitialMomentum(gp);
  track.SetInitialMomentum( 1.2 );
  //  track.SetInitialMomentum( 1.3 );
  track.doFit();
  bool status=track.Status();
  double chisqr=track.chisqr();

  ThreeVector priPos = track.PrimaryPosition();
  ThreeVector priMom = track.PrimaryMomentum();

  double rx=priPos.x(), ry=priPos.y(), rz=priPos.z();
  double rp=priMom.mag();
  double rt=priMom.theta()*Rad2Deg, rf=priMom.phi()*Rad2Deg;
  double lengTof = track.PathLengthToTOF();

  int nhdc = track.GetNHits();

#ifdef Interactive
  {
    std::ios::fmtflags oldFlags = std::cout.flags();
    std::size_t oldPrec = std::cout.precision();
    std::cout.setf( std::ios::fixed );
    std::cout.precision(4);

    std::cout << "[TestQSksTracking]: St=" << status
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
	      << std::setw(8) << gp << " ) < "
	      << std::setw(8) << rp-gp << " > [GeV/c]";
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


  HF1( 101, 1. );

  HF1( 102, status );
  HF1( 103, chisqr );
  HF1( 104, nhdc ); 
  HF1( 105, rx ); HF1( 106, ry ); HF1( 107, rp );
  HF1( 108, rt ); HF1( 109, rf );
  HF1( 110, lengTof );
  HF1( 111, rx-gx ); HF1( 112, ry-gy ); HF1( 113, rp-gp );
  HF1( 114, rt-gt ); HF1( 115, rf-gf );
  HF1( 116, ru-gu ); HF1( 117, rv-gv );
  HF2( 118, gp, rp-gp );

  for( std::size_t i=0; i<nhdc; ++i ){
    TrackHit *thit = track.GetHit(i);
    if( thit ){
      int layer = thit->GetLayer();
      double posl = thit->GetLocalHitPos();
      double resd = thit->GetResidual();
      HF1( 200+layer, posl );
      HF1( 300+layer, resd );
      HF2( 400+layer, posl, resd );
    }
  }

  if( status ){
    HF1( 101, 2. );
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
      HF1( 101, 3. );
      HF1( 163, chisqr );
      HF1( 164, nhdc ); 
      HF1( 165, rx ); HF1( 166, ry ); HF1( 167, rp );
      HF1( 168, rt ); HF1( 169, rf );
      HF1( 170, lengTof );
      HF1( 171, rx-gx ); HF1( 172, ry-gy ); HF1( 173, rp-gp );
      HF1( 174, rt-gt ); HF1( 175, rf-gf );
      HF1( 176, ru-gu ); HF1( 177, rv-gv );
      HF2( 178, gp, rp-gp );

      for( std::size_t i=0; i<nhdc; ++i ){
	TrackHit *thit = track.GetHit(i);
	if( thit ){
	  int layer = thit->GetLayer();
	  double posl = thit->GetLocalHitPos();
	  double resd = thit->GetResidual();
	  HF1( 1200+layer, posl );
	  HF1( 1300+layer, resd );
	  HF2( 1400+layer, posl, resd );
	}
      }
    }
  }
  else {
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

  for( std::size_t i=0; i<nhin; ++i )
    delete trIn.GetHit(i);
  for( std::size_t i=0; i<nhout; ++i )
    delete trOut.GetHit(i);

  return true;
}




















