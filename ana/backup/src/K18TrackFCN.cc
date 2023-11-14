/*
  K18TrackFCN.cc
*/

#include "K18TrackFCN.hh"
#include "K18Track.hh"
#include "K18TransMatrix.hh"
#include "TrackHit.hh"
//#include "DCGeomMan.hh"
#include "DCConfMan.hh"

#include <string>
#include <cmath>

const double Deg2Rad = acos(-1.)/180.;
const double Rad2Deg = 180./acos(-1.);

K18TrackFCN::K18TrackFCN( K18Track *track, K18TransMatrix *trMatrix )
  : Tr_(track), trM_(trMatrix)
{}

K18TrackFCN::~K18TrackFCN()
{}

double K18TrackFCN::operator()( int np, double *g, double *u, int flag )
{
  static const std::string funcname = "[K18TrackFCN::operator()]";
  const DCConfMan & geomMan=DCConfMan::GetInstance();

  double xi=u[0], yi=u[1], ui=u[2], vi=u[3], delta=u[4];
  double xo, yo, uo, vo; 

  trM_->Transport( xi, yi, ui, vi, delta, xo, yo, uo, vo );
  /*
  std::cout << "xi = " << xi << ", yi = " << yi << ", ui = " << ui
	    << ", vi = " << vi << ", delta = " << delta 
	    << "xo = " << xo << ", yo = " << yo 
	    << ", uo = " << uo << ", vo = " << vo << std::endl;
  */
  int nin=Tr_->GetNHitsIn(), nout=Tr_->GetNHitsOut();

  double chi=0.0;
  int nh=0;

  for( int i=0; i<nin; ++i ){
    TrackHit *hit=Tr_->GetK18HitIn(i);
    if(!hit) continue;
    int lnum=hit->GetLayer();
    //if(lnum==131||lnum==132) std::cout<<"lnum="<<lnum<<std::endl;
    //std::cout<<"try"<<std::endl;
    if(lnum==103||lnum==104/*||lnum==110||lnum==111*/){
    //if(lnum==131||lnum==132){
      //std::cout<<"hit"<<std::endl;
      double z=geomMan.GetLocalZ( lnum );
      double dd=geomMan.GetResolution( lnum );
      double pos=hit->GetLocalHitPos();
      double aa=hit->GetDCHit()->GetTiltAngle()*Deg2Rad;
      double cpx=xi+ui*z, cpy=yi+vi*z;
      double c=cpx*cos(aa)+cpy*sin(aa);
      chi += ((pos-c)*(pos-c)/(dd*dd)); ++nh;
      if( flag==3 ) hit->SetCalLPos(c);
    }
  }

  for( int i=0; i<nout; ++i ){
    TrackHit *hit=Tr_->GetK18HitOut(i);
    if(!hit) continue;
    int lnum=hit->GetLayer();
    //if(lnum==115||lnum==116||lnum==122||lnum==123){
    if(1){
      double z=geomMan.GetLocalZ( lnum );
      double dd=geomMan.GetResolution( lnum );
      double pos=hit->GetLocalHitPos();
      double aa=hit->GetDCHit()->GetTiltAngle()*Deg2Rad;
      double cpx=xo+uo*z, cpy=yo+vo*z;
      double c=cpx*cos(aa)+cpy*sin(aa);
      chi += ((pos-c)*(pos-c)/(dd*dd)); ++nh;
      if( flag==3 ) hit->SetCalLPos(c);
    }
  }

  if( nh>5 ) chi /= double(nh-5);

  //std::cout<<"return chi="<<chi<<std::endl;
  return chi;
}




//if 0


#if 0
double K18TrackFCN::operator()( int np, double *g, double *u, int flag )
{
  static const std::string funcname = "[K18TrackFCN::operator()]";
  const DCConfMan & geomMan=DCConfMan::GetInstance();

  double xi=-u[0], yi=-u[1], ui=u[2], vi=u[3], delta=u[4];
  double xo, yo, uo, vo; 

  trM_->Transport( xi, yi, ui, vi, delta, xo, yo, uo, vo );
  /*
  std::cout << "xi = " << xi << ", yi = " << yi << ", ui = " << ui
	    << ", vi = " << vi << ", delta = " << delta 
	    << "xo = " << xo << ", yo = " << yo 
	    << ", uo = " << uo << ", vo = " << vo << std::endl;
  */
  int nin=Tr_->GetNHitsIn(), nout=Tr_->GetNHitsOut();

  double chi=0.0;
  int nh=0;

  for( int i=0; i<nin; ++i ){
    TrackHit *hit=Tr_->GetK18HitIn(i);
    if(!hit) continue;
    int lnum=hit->GetLayer();
    double z= -geomMan.GetLocalZ( lnum );
    double dd=geomMan.GetResolution( lnum );
    double pos= -hit->GetLocalHitPos();
    double aa=hit->GetDCHit()->GetTiltAngle()*Deg2Rad;
    double cpx=xi+ui*z, cpy=yi+vi*z;
    double c=cpx*cos(aa)+cpy*sin(aa);
    chi += ((pos-c)*(pos-c)/(dd*dd)); ++nh;
    if( flag==3 ) hit->SetCalLPos(c);
  }

  for( int i=0; i<nout; ++i ){
    TrackHit *hit=Tr_->GetK18HitOut(i);
    if(!hit) continue;
    int lnum=hit->GetLayer();
    double z= -geomMan.GetLocalZ( lnum );
    double dd=geomMan.GetResolution( lnum );
    double pos= -hit->GetLocalHitPos();
    double aa=hit->GetDCHit()->GetTiltAngle()*Deg2Rad;
    double cpx=xo+uo*z, cpy=yo+vo*z;
    double c=cpx*cos(aa)+cpy*sin(aa);
    chi += ((pos-c)*(pos-c)/(dd*dd)); ++nh;
    if( flag==3 ) hit->SetCalLPos(c);
  }

  if( nh>5 ) chi /= double(nh-5);

  return chi;
}
#endif
