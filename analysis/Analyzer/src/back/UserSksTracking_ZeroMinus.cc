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
 
#include "RootHelper.hh"

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <string>
#include <vector>

//#define Interactive 1

////////////////Physics Parameters///////////////////
const double AtomicMassUnit  = 0.93149432;
const double PionMass        = 0.1395701;
const double KaonMass        = 0.493677;
const double ProtonMass      = 0.93827200;
const double NeutronMass     = 0.93956563;
const double LambdaMass      = 1.115648;
//2H
const double DeutronMass     = 2.*AtomicMassUnit+0.01313672;
//4He
const double Nucl4HeMass     = 4.0*AtomicMassUnit+0.0024249;
const double DeltaHe3        = 0.0149312;
const double He3CoreMass     = 3.*AtomicMassUnit+DeltaHe3;
//7Li
const double Nucl7LiMass     = 7.0*AtomicMassUnit+0.014908;
const double DeltaLi6        = 0.014086;
const double Li6CoreMass     = 6.*AtomicMassUnit+DeltaLi6;
//12C
const double Nucl12CMass      = 12.*AtomicMassUnit;
const double DeltaC11        = 0.0106502;
const double C11CoreMass     = 11.*AtomicMassUnit+DeltaC11;
////////////////////////////////////////////////////

const double Deg2Rad = acos(-1.)/180.;
const double Rad2Deg = 180./acos(-1.);

bool ProcessOneEvent( std::ifstream & );
void DefineHistograms( const char * );
void InitializeEvent( void );

const double MaxChiSqr = 100.;

#ifndef MaxHits 
#define MaxHits 30
#endif

struct Event{
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
};
static Event event;

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

  //Tree
  TTree *tree = dynamic_cast<TTree *>(gFile->Get("tree"));
  InitializeEvent();

  if( !In ) return false;
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

  DCLocalTrack trIn, trOut;
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
    In >> pos >> lnum;
#if 0
      std::cout<< "*********************************************" <<std::endl;
  std::cout<<" "<< pos 
	   <<" "<< lnum 
	   <<" "<< lnumold 
	   <<std::endl;
  std::cout<< "*********************************************" <<std::endl;
#endif
    if( lnum<60 && lnumold<60 ){
      DCHit *hit = new DCHit( lnumold, 0, 0 );
      if( hit ){
	hit->SetLocalHitPos( pos );
	hit->SetTiltAngle( DCConfMan::GetInstance().GetTiltAngle(lnumold) );
	if( lnumold>0 && lnumold<11 )
	  trIn.AddHit( hit );
	else if( lnumold>40 && lnumold<55 )
	  trOut.AddHit( hit );
	else
	  delete hit;
      }
    }

    //TOF
    if( lnum == 62 && lnumold == 61 ){
      toftime = pos;
      event.tofTime[tofhits] = toftime; 
      tofhits++;
    }
    if( (lnum == 61 && lnumold == 62) ||
	(lnum == -1 && lnumold == 62) ){
      tofedep = pos;
      event.tofEdep[tofhits-1] = tofedep;
    }
    //SMF
    if( lnum == 72 && lnumold == 71 ){
      smftime = pos;
      event.smfTime[smfhits] = smftime; 
      smfhits++;
    }
    if( (lnum == 71 && lnumold == 72) ||
	(lnum == -1 && lnumold == 72) ){
      smfedep = pos;
      event.smfEdep[smfhits-1] = smfedep;
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

  HF1( 101, 1. );

  trIn.DoFit(); trOut.DoFit();
  bool stIn=trIn.GetStatus(), stOut=trOut.GetStatus();
  int nhin=trIn.GetNHit(), nhout=trOut.GetNHit();

  if( !stIn || !stOut ){
    for( int i=0; i<nhin; ++i )
      delete trIn.GetHit(i);
    for( int i=0; i<nhout; ++i )
      delete trOut.GetHit(i);
    return true; 
  }

  HF1( 101, 2. );

  SksTrack track( &trIn, &trOut );
  track.SetInitialMomentum( 0.9 );
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

  //#ifdef Interactive
#if 1
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
  
  bool tofflag = false;
  bool smfflag = false;

  if( tofedep>4 && 12<toftime && toftime<24 ) tofflag = true;
  if( smfedep>6 && 12<smftime && smftime<30 ) smfflag = true;

  //Missing Mass
  ////Beam Mom
  double Pbz = gpb/sqrt(1.0+gub*gub+gvb*gvb);
  double Pbx = Pbz*gub;
  double Pby = Pbz*gvb;
  ThreeVector bMom( Pbx, Pby ,Pbz );
  ////Scattered Mom
  double Pz = rp/sqrt(1.0+ru*ru+rv*rv);
  double Px = Pz*ru;
  double Py = Pz*rv;
  ThreeVector sMom( Px, Py ,Pz );

  LorentzVector LvK(  bMom, sqrt(KaonMass*KaonMass+bMom.mag2()) );
  LorentzVector LvPi( sMom, sqrt(PionMass*PionMass+sMom.mag2()) );
//   LorentzVector LvTgt( 0., 0., 0., ProtonMass );
//   LorentzVector LvCore( 0., 0., 0., 0.0 );
  LorentzVector LvTgt( 0., 0., 0., Nucl7LiMass );
  LorentzVector LvCore( 0., 0., 0., Li6CoreMass );
  LorentzVector LvRc = LvK+LvTgt-LvPi;
  double recoilmom = (LvK.vect()-LvPi.vect()).mag();
  double mismass = LvRc.mag()-LvCore.mag();
  double BLambda = LvRc.mag()-( LvCore.mag()+LambdaMass );
  event.RecoilMom = recoilmom;
  event.MissMass = mismass;
  event.BE = BLambda;

  if( tofflag && !smfflag ) event.BE_SMF = BLambda;

  int nhit1=0;
  int nhit2=0;
  int nhit3=0;
  int nhit4=0;
  for( int i=0; i<nhdc; ++i ){
    TrackHit *thit = track.GetHit(i);
    if( thit ){
      int layer = thit->GetLayer();
      double posl = thit->GetLocalHitPos();
      double resd = thit->GetResidual();
      //std::cout << "Residual=" << resd << std::endl;
      HF1( 200+layer, posl );
      HF1( 300+layer, resd );
      HF2( 400+layer, posl, resd );

      if( layer==1 ){
	event.DC1Posl[nhit1] = posl;
	event.DC1Resd[nhit1] = resd;
	nhit1++;
      }
      if( layer==5 ){
	event.DC2Posl[nhit2] = posl;
	event.DC2Resd[nhit2] = resd;
	nhit2++;
      }
      if( layer==41 ){
	event.DC3Posl[nhit3] = posl;
	event.DC3Resd[nhit3] = resd;
	nhit3++;     
      }
      if( layer==45 ){
	event.DC4Posl[nhit4] = posl;
	event.DC4Resd[nhit4] = resd;
	nhit4++;      
      }
    }
  }//std::cout<<"********************"<<std::endl;
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
	TrackHit *thit = track.GetHit(i);
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

  for( int i=0; i<nhin; ++i ){
    delete trIn.GetHit(i);
  }
  for( int i=0; i<nhout; ++i ){
    delete trOut.GetHit(i);
  }
  //std::cout<<"********************"<<std::endl;

  return true;
}

void InitializeEvent( void )
{
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
}


const int NBin1HPDcIn[11] = {
  100, 
  200, 200, 200, 200, 200, 200, 
  300, 300, 300, 300
};

const int NBin2HPDcIn[11] = {
  100, 
  50, 50, 50, 50, 
  100, 100, 100, 100, 100, 100
};

const double MinHPDcIn[11] = {
  0.0, 
  // SDC1
  -200.0, -200.0, -200.0, -200.0, -200.0, -200.0,
  // SDC2
  -300.0, -300.0, -300.0, -300.0
};

const double MaxHPDcIn[11] = {
  0.0,
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

  for( int i=1; i<11; ++i ){
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
    HB1( 200+i+20, title1.str().c_str(), NBin1HPDcOut[i-10], 
	 MinHPDcOut[i-10], MaxHPDcOut[i-10] );

    title2 << "Residual [mm] Layer#" << i;
    HB1( 300+i+20, title2.str().c_str(), 250, -2.5, 2.5 );

    title3 << "Resd.%Pos [mm] Layer#" << i;
    HB2( 400+i+20, title3.str().c_str(), NBin2HPDcOut[i-10], 
	 MinHPDcOut[i-10], MaxHPDcOut[i-10], 100, -2.5, 2.5 );

    title4 << " <GoodTrack> Layer#" << i;
    HB1( 1200+i+20, title1.str().c_str(), NBin1HPDcOut[i-10], 
	 MinHPDcOut[i-10], MaxHPDcOut[i-10] );

    title5 << " <GoodTrack> Layer#" << i;
    HB1( 1300+i+20, title2.str().c_str(), 250, -2.5, 2.5 );

    title6 << " <GoodTrack> Layer#" << i;
    HB2( 1400+i+20, title3.str().c_str(), NBin2HPDcOut[i-10], 
	 MinHPDcOut[i-10], MaxHPDcOut[i-10], 100, -2.5, 2.5 );
  }

  //Tree
  HBTree("tree","tree of Sks");
  TTree *tree = dynamic_cast<TTree *>(gFile->Get("tree"));

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
}
