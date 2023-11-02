/*
  K18TrackD2UFCN.cc
*/

#include "K18TrackD2UFCN.hh"
#include "K18Track.hh"
#include "K18TransMatrix.hh"
#include "TrackHit.hh"
//#include "DCGeomMan.hh"
#include "DCConfMan.hh"

#include <string>
#include <cmath>

const double Deg2Rad = acos(-1.)/180.;
const double Rad2Deg = 180./acos(-1.);

K18TrackD2UFCN::K18TrackD2UFCN( K18Track *track, K18TransMatrix *trMatrix )
  : Tr_(track), trM_(trMatrix)
{}

K18TrackD2UFCN::~K18TrackD2UFCN()
{}

double K18TrackD2UFCN::operator()( int np, double *g, double *u, int flag )
{
  static const std::string funcname = "[K18TrackD2UFCN::operator()]";
  const DCConfMan & geomMan=DCConfMan::GetInstance();

  double xo=-u[0], yo=u[1], uo=u[2], vo=-u[3], delta=u[4];//convert to TRANS
  double xi, yi, ui, vi; 


  trM_->Transport( xo, yo, uo, vo, delta, xi, yi, ui, vi );
  /*
  xi *= -1.;
  yi *= -1.;
  ui *= -1.;
  vi *= -1.;
  */
  //in x v
  xo *= -1.;
  vo *= -1.;
  //out x v
  xi *= -1.;
  vi *= -1.;

  /*
  std::cout << "xo = " << xo << ", yo = " << yo << ", uo = " << uo
	    << ", vo = " << vo << ", delta = " << delta 
	    << "xi = " << xi << ", yi = " << yi 
	    << ", ui " << ui << ", vi = " << vi << std::endl;
  */
  int nin=Tr_->GetNHitsIn(), nout=Tr_->GetNHitsOut();

  double chi=0.0;
  int nh=0;

  for( int i=0; i<nout; ++i ){
    TrackHit *hit=Tr_->GetK18HitOut(i);
    if(!hit) continue;
    int lnum=hit->GetLayer();
    //std::cout<<"out lnum="<<lnum<<std::endl;
    double z= geomMan.GetLocalZ( lnum );
    double dd=geomMan.GetResolution( lnum );
    double pos= hit->GetLocalHitPos();
    double aa=hit->GetDCHit()->GetTiltAngle()*Deg2Rad;
    // x --> reverse dehanaku sonomama, y --> sonomama
    double cpx=xo+uo*z, cpy=yo+vo*z;
    double c=cpx*cos(aa)+cpy*sin(aa);
    chi += ((pos-c)*(pos-c)/(dd*dd)); ++nh;
    if( flag==3 ) hit->SetCalLPos(c);
  }
  
  for( int i=0; i<nin; ++i ){
    TrackHit *hit=Tr_->GetK18HitIn(i);
    if(!hit) continue;
    int lnum=hit->GetLayer();
    //if(1){
    if (lnum == 131 || lnum == 132) {
      //std::cout<<"lnum131"<<std::endl;
      double z= geomMan.GetLocalZ( lnum );
      double dd=geomMan.GetResolution( lnum );
      double pos= hit->GetLocalHitPos();
      double aa=hit->GetDCHit()->GetTiltAngle()*Deg2Rad;
      // x --> reverse dehanaku sonomama, y --> sonomama
      double cpx=xi+ui*z, cpy=yi+vi*z;
      double c=cpx*cos(aa)+cpy*sin(aa);
      chi += ((pos-c)*(pos-c)/(dd*dd)); ++nh;
      //std::cout<<"lnum="<<lnum<<" chi="<<chi<<std::endl;
      if( flag==3 ) hit->SetCalLPos(c);
    }
  }

  if( nh>5 ) chi /= double(nh-5);

  return chi;
}
