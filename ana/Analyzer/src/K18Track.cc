/*
  K18Track.cc
*/

#include "K18Track.hh"
#include "DCLocalTrack.hh"
#include "K18TransMatrix.hh"
#include "DCConfMan.hh"
#include "DetectorId.hh"
#include "TrackHit.hh"
#include "TemplateLib.hh"
#include "K18TrackFCN.hh"
#include "K18TrackD2UFCN.hh"
#include "Minuit.hh"
//#include "ConfMan.hh"
#include "K18Parameters.hh"
//#include "SksObjectId.hh"
#include "DCAnalyzer.hh"

#include <string>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <cstring>
#include <stdexcept>
#include <sstream>

const double LowBand[5] = 
  { MinK18InX, MinK18InY, MinK18InU, MinK18InV, MinK18Delta };

const double UpperBand[5] = 
  { MaxK18InX, MaxK18InY, MaxK18InU, MaxK18InV, MaxK18Delta };

const double LowBandOut[5] = 
  { MinK18OutX, MinK18OutY, MinK18OutU, MinK18OutV, MinK18Delta };

const double UpperBandOut[5] = 
  { MaxK18OutX, MaxK18OutY, MaxK18OutU, MaxK18OutV, MaxK18Delta };

const double IniError[5] =
  { 0.3, 0.3, 1.E-4, 1.E-4, 1.E-3 };
//  { 0.4, 0.4, 1.E-5, 1.E-4, 1.E-3 };

const int MaxFCNCall = 300;
const double EPS = 1.;  


K18Track::K18Track( DCLocalTrack *tin, DCLocalTrack *tout, double P0, DCHit *SF1, DCHit *SF2)
  : TrIn_(tin), TrOut_(tout), P0_(P0), sf1_(SF1), sf2_(SF2),
    Status_(false), StatusD2U_(false), gfastatus_(true)
{}

K18Track::~K18Track()
{
  deleteHits();
}

bool K18Track::doFit( void )
{
  static const std::string funcname = "[K18Track::doFit]";
  Status_=false;
  deleteHits();
  addHits();

  //K18TransMatrix *mat=ConfMan::GetConfManager()->GetK18Matrix();
  K18TransMatrix *mat=K18TransMatrix::GetTransMatrix();

  K18TrackFCN FCN( this, mat );

  double param[5], error[5];
  param[0]=TrIn_->GetX0(); 
  param[1]=TrIn_->GetY0();
  param[2]=TrIn_->GetU0();
  param[3]=TrIn_->GetV0();
  param[4]=0.0;
  error[0]=IniError[0];
  error[1]=IniError[1];
  error[2]=IniError[2];


  error[3]=IniError[3];
  error[4]=IniError[4];

  double LowBand_[5], UpperBand_[5];
  for( int i=0; i<5; ++i ){
    LowBand_[i]=LowBand[i]; 
    UpperBand_[i]=UpperBand[i];
    if( param[i]<LowBand_[i] )   LowBand_[i]=param[i];
    if( param[i]>UpperBand_[i] ) UpperBand_[i]=param[i];
  }

#if 0
  std::cout << funcname << ": before fitting." << std::endl;
#endif
  
  Status_ = Minuit(&FCN).Fit( 5, param, error, LowBand_, UpperBand_,
			      MaxFCNCall, EPS, chisqr_ );

  //std::cout<<"final chi="<<chisqr_<<std::endl;
#if 0
  std::cout << funcname << ": after fitting. " 
	    << " Status=" << Status_  << std::endl;
#endif

  if( Status_ ){
    Xi_=param[0]; Yi_=param[1]; Ui_=param[2]; Vi_=param[3];
    Delta_=param[4];
    mat->Transport( Xi_, Yi_, Ui_, Vi_, Delta_,
		    Xo_, Yo_, Uo_, Vo_ );
  }

  return Status_;
}

bool K18Track::doFitD2U( void )
{
  static const std::string funcname = "[K18Track::doFitD2U]";
  Status_=false;
  deleteHits();
  addHits();

  //K18TransMatrix *mat=ConfMan::GetConfManager()->GetK18Matrix();
  K18TransMatrix *mat=K18TransMatrix::GetTransMatrix();

  K18TrackD2UFCN FCN( this, mat );

  double param[5], error[5];
  param[0]=TrOut_->GetX0(); 
  param[1]=TrOut_->GetY0();
  param[2]=TrOut_->GetU0();
  param[3]=TrOut_->GetV0();
  param[4]=0.0;
  error[0]=IniError[0];
  error[1]=IniError[1];
  error[2]=IniError[2];

  error[3]=IniError[3];
  error[4]=IniError[4];

  double LowBand_[5], UpperBand_[5];
  for( int i=0; i<5; ++i ){
    LowBand_[i]=LowBandOut[i]; 
    UpperBand_[i]=UpperBandOut[i];
    if( param[i]<LowBand_[i] )   LowBand_[i]=param[i];
    if( param[i]>UpperBand_[i] ) UpperBand_[i]=param[i];
  }

#if 0
  std::cout << funcname << ": before fitting." << std::endl;
#endif
  
  Status_ = Minuit(&FCN).Fit( 5, param, error, LowBand_, UpperBand_,
			      MaxFCNCall, EPS, chisqr_ );

#if 0
  std::cout << funcname << ": after fitting. " 
	    << " Status=" << Status_  << std::endl;
#endif

  if( Status_ ){
    Xo_=param[0]; Yo_=param[1]; Uo_=param[2]; Vo_=param[3];
    Delta_=param[4];
    mat->Transport( -Xo_, Yo_, Uo_, -Vo_, Delta_,
		    Xi_, Yi_, Ui_, Vi_ );
    Xi_ *= -1.;
    Vi_ *= -1.;
  }

  return Status_;
}

bool K18Track::CalcMomentumD2U( void )
{
  static const std::string funcname = "[K18Track::CalcMomentumD2U]";
  StatusD2U_=false;

  //K18TransMatrix *mat=ConfMan::GetConfManager()->GetK18Matrix();
  K18TransMatrix *mat=K18TransMatrix::GetTransMatrix();

  double xi, yi, ui, vi, xo, delta1, delta2;
  xi= -TrOut_->GetX0(); 
  yi=  TrOut_->GetY0();
  ui=  TrOut_->GetU0();
  vi= -TrOut_->GetV0();
  xo= -TrIn_->GetX0(); 

  /*
  DCHit *hit=TrIn_->GetHitOfLayerNumber(112);
  if (hit) {
    xo = -hit->GetLocalHitPos();
  }
  */
  StatusD2U_ = mat->CalcDeltaD2U(xi, yi, ui, vi, xo, delta1, delta2);

#if 0
  std::cout << funcname << ": after calculation. " 
	    << " StatusD2U=" << StatusD2U_  << std::endl;
#endif

  if( StatusD2U_ ){
    //std::cout << "delta1 = " << delta1 << ", delta2 = " << delta2 << std::endl;
    DeltaD2U_=delta1;
  }

  return StatusD2U_;
}

ThreeVector K18Track::BeamMomentum( void ) const
{
  double u=TrOut_->GetU0(), v=TrOut_->GetV0();
  double pz=P()/sqrt(1.+u*u+v*v);

  return ThreeVector( pz*u, pz*v, pz );
}

double K18Track::Xtgt( void ) const
{
  double z=DCConfMan::GetInstance().GetLocalZ( IdK18Target );
  return TrOut_->GetU0()*z+TrOut_->GetX0();
}

double K18Track::Ytgt( void ) const
{
  double z=DCConfMan::GetInstance().GetLocalZ( IdK18Target );
  return TrOut_->GetV0()*z+TrOut_->GetY0();
}

double K18Track::Utgt( void ) const
{
  return TrOut_->GetU0();
}

double K18Track::Vtgt( void ) const
{
  return TrOut_->GetV0();
}

void K18Track::deleteHits( void )
{
  for_each( hitContIn.begin(),  hitContIn.end(),  DeleteObject() );
  for_each( hitContOut.begin(), hitContOut.end(), DeleteObject() );
  hitContIn.clear(); hitContOut.clear();
}

void K18Track::addHits( void )
{
  static const std::string funcname="K18Track::addHits";

  int nhIn=TrIn_->GetNHit(), nhOut=TrOut_->GetNHit();

  for( int i=0; i<nhIn; ++i ){
    DCHit *lhit = TrIn_->GetHit(i);
    if(!lhit) continue;
    TrackHit *hit = new TrackHit(lhit);
    if(hit){
      hitContIn.push_back(hit);
    }
    else{
      std::cerr << funcname << ": new fail" << std::endl;
    }      
  }

  {
    TrackHit *hit = new TrackHit(sf1_);
    if(hit){
      hitContIn.push_back(hit);
    }
    else{
      std::cerr << funcname << ": new fail" << std::endl;
    }      
  }
  {
    TrackHit *hit = new TrackHit(sf2_);
    if(hit){
      hitContIn.push_back(hit);
    }
    else{
      std::cerr << funcname << ": new fail" << std::endl;
    }      
  }

  for( int i=0; i<nhOut; ++i ){
    DCHit *lhit = TrOut_->GetHit(i);
    if(!lhit) continue;
    TrackHit *hit = new TrackHit(lhit);
    //std::cout<<"layer "<<hit->GetLayer()<<std::endl;
    if(hit){
      hitContOut.push_back(hit);
    }
    else{
      std::cerr << funcname << ": new fail" << std::endl;
    }      
  }
}

TrackHit * K18Track::GetK18HitIn( int i )
{
  if( i>=0 && i<hitContIn.size() )
    return hitContIn[i];
  else
    return 0;
}

TrackHit * K18Track::GetK18HitOut( int i )
{
  if( i>=0 && i<hitContOut.size() )
    return hitContOut[i];
  else
    return 0;
}

TrackHit * K18Track::GetK18HitTotal( int i )
{
  int nin=hitContIn.size();
  if( i<0 ) 
    return 0;
  else if( i<nin )
    return hitContIn[i];
  else if( i<nin+hitContOut.size() )
    return hitContOut[i-nin];
  else
    return 0;
} 


TrackHit *K18Track::GetK18HitByPlaneId( int PlId )
{
  int id = PlId-PlOffsBc;
  TrackHit *hit=0;
  if( id>=PlMinBcIn && id<=PlMaxBcIn ){
    for( int i=0; i<hitContIn.size(); ++i ){
      if( hitContIn[i]->GetLayer()==PlId ){
	hit=hitContIn[i]; break;
      }
    }
  }
  else if( id>=PlMinBcOut && id<=PlMaxBcOut ){
    for( int i=0; i<hitContOut.size(); ++i ){
      if( hitContOut[i]->GetLayer()==PlId ){
	hit=hitContOut[i]; break;
      }
    }
  }
  return hit;
}

/*
bool K18Track::ReCalc( bool applyRecursively )
{
  static const std::string funcname = "[K18Track::ReCalc]"; 

  bool ret1=true, ret2=true, ret=false;
  if( applyRecursively ){
    ret1=TrIn_->ReCalc(applyRecursively);
    ret2=TrOut_->ReCalc(applyRecursively);
  }
  if( ret1 && ret2 ){
    ret=doFit();
  }
  if(!ret){
    std::cerr << funcname << ": ReCalculation fails" << std::endl;
  }

  return ret;
}
*/

////////////
// for DS //
////////////

// std::size_t K18Track::DSSize( void ) const
// {
//   static const std::string funcname="[K18Track::DSSize]";

//   std::size_t size=sizeof(DSFormat)/sizeof(int);
//   if( sizeof(DSFormat)%sizeof(int) ) size+=1;

//   int nhIn=hitContIn.size(), nhOut=hitContOut.size();
//   size+=1; 
//   for( int i=0; i<nhIn; ++i ){
//     if( TrackHit *hit=hitContIn[i] )
//       size += hit->DSSize();
//   }
//   size+=1; 
//   for( int i=0; i<nhOut; ++i ){
//     if( TrackHit *hit=hitContOut[i] )
//       size += hit->DSSize();
//   }

// #if 0
//   std::cout << funcname << ": Size=" << size << std::endl;
// #endif

//   return size;
// }

// std::size_t K18Track::DSSave( unsigned int *bufp ) const
// {
//   static const std::string funcname = "[K18Track::DSSave]";
//   DSFormat data;
//   std::size_t size=sizeof(data)/sizeof(int);
//   if( sizeof(data)%sizeof(int) ) size+=1;


//   std::size_t sizeIn=1, sizeOut=1;
//   int nhIn=hitContIn.size(), nhOut=hitContOut.size();
//   for( int i=0; i<nhIn; ++i ){
//     if( TrackHit *hit=hitContIn[i] )
//       sizeIn += hit->DSSize();
//   }
//   for( int i=0; i<nhOut; ++i ){
//     if( TrackHit *hit=hitContOut[i] )
//       sizeOut += hit->DSSize();
//   }
//   std::size_t totalSize = size+sizeIn+sizeOut;
  
//   int headerId=DSMakeHeaderID( SksObjK18Track );
//   data.header_=DSMakeHeader( headerId, totalSize );
//   data.pkeyIn_ =TrIn_->DSGetSaveKey();
//   data.pkeyOut_=TrOut_->DSGetSaveKey();
//   data.xi_=Xi_; data.yi_=Yi_; data.ui_=Ui_; data.vi_=Vi_;
//   data.xo_=Xo_; data.yo_=Yo_; data.uo_=Uo_; data.vo_=Vo_;
//   data.p0_=P0_; data.delta_=Delta_; data.chisqr_=chisqr_;
//   if( Status_ ) data.status_=1;
//   else          data.status_=0;
//   if( gfastatus_ ) data.gfastatus_=1;
//   else             data.gfastatus_=0;

//   std::memcpy( bufp, &data, sizeof(data) );
//   bufp+=size;

//   int headerIn =DSMakeHeaderID(SksObjTrackHitInC,true);
//   int headerOut=DSMakeHeaderID(SksObjTrackHitOutC,true);

//   *bufp=DSMakeHeader( headerIn, sizeIn );
//   bufp+=1;
//   for( int i=0; i<nhIn; ++i ){
//     if( TrackHit *hit=hitContIn[i] ){
//       std::size_t s = hit->DSSave(bufp);
//       hit->DSSetSaveKey( DSMakeKey( headerIn, i ) ); 
//       bufp += s;
//     }
//   }
//   *bufp=DSMakeHeader( headerOut, sizeOut );
//   bufp+=1;
//   for( int i=0; i<nhOut; ++i ){
//     if( TrackHit *hit=hitContOut[i] ){
//       std::size_t s = hit->DSSave(bufp);
//       hit->DSSetSaveKey( DSMakeKey( headerOut, i ) ); 
//       bufp += s;
//     }
//   }

// #if 0
//   std::cout << funcname << ": Size=" << totalSize << std::endl;
// #endif

//   return totalSize;

// } 

// bool K18Track::DSRestore( unsigned int *bufp, DCAnalyzer *DCana )
// {
//   static const std::string funcname = "[K18Track::DSRestore]";

//   DSFormat data;
//   unsigned int *hp=bufp;
//   std::memcpy( &data, bufp, sizeof(data) );
//   std::size_t size=sizeof(data)/sizeof(int);
//   if( sizeof(data)%sizeof(int) ) size+=1;
//   bufp+=size;

//   std::size_t totalSize=DSGetSize( data.header_ );

//   if( totalSize==0 ){
//     std::ostringstream mess;
//     mess << funcname << ": Size=0";
//     std::cerr << mess.str() << std::endl;
//     throw std::invalid_argument(mess.str());
//   }

//   bool ret=true;

//   Xi_=data.xi_; Yi_=data.yi_; Ui_=data.ui_; Vi_=data.vi_;
//   Xo_=data.xo_; Yo_=data.yo_; Uo_=data.uo_; Vi_=data.vo_;
//   P0_=data.p0_; Delta_=data.delta_; chisqr_=data.chisqr_;
//   if( data.status_ ) Status_=true;
//   else               Status_=false;
//   if( data.gfastatus_ ) gfastatus_=true;
//   else                  gfastatus_=false;

//   int HCIdIn =DSGetObjectID( DSGetHeaderID( data.pkeyIn_ ) );
//   int HCIdOut=DSGetObjectID( DSGetHeaderID( data.pkeyOut_ ) );
//   int seqNumIn = DSGetSeqNum( data.pkeyIn_ );
//   int seqNumOut= DSGetSeqNum( data.pkeyOut_ );

//   TrIn_=TrOut_=0;
//   if( HCIdIn==SksObjLocalTrackBIC )  TrIn_=DCana->GetTrackBcIn(seqNumIn);
//   if( HCIdOut==SksObjLocalTrackBOC ) TrOut_=DCana->GetTrackBcOut(seqNumOut);

//   deleteHits();

//   // HitIn
//   unsigned int *hpIn=bufp;
//   std::size_t sizeIn=DSGetSize( *hpIn ); 
//   int objIdIn=DSGetObjectID( DSGetHeaderID( *hpIn ) );
//   bufp+=1;
//   while( bufp<hpIn+sizeIn ){
//     std::size_t s1=DSGetSize( *bufp );
//     int objId=DSGetObjectID( DSGetHeaderID( *bufp ) );
// #if 0    
//     std::cout << funcname << ": Size=" << std::setw(4) << size
//               << " objId=" << std::setw(3) << objId << std::endl;
// #endif
//     if( objId==SksObjTrackHit && objIdIn==SksObjTrackHitInC && TrIn_){
//       TrackHit *hit = new TrackHit( bufp, TrIn_, TrOut_ );
//       if( hit ) hitContIn.push_back(hit);
//       else {
// 	std::cerr << funcname << ": new fail objId=" 
//                   << std::setw(3) << objId << std::endl;
//         ret=false;
//       }
//     }
//     else{
//       std::cerr << funcname << ": invalid object type " 
//                 << std::setw(3) << objId << std::endl;
//       ret=false;
//     }
//     bufp+=s1;
//   }
//   // HitOut
//   unsigned int *hpOut=bufp;
//   std::size_t sizeOut=DSGetSize( *hpOut ); 
//   int objIdOut=DSGetObjectID( DSGetHeaderID( *hpOut ) );
//   bufp+=1;
//   while( bufp<hpOut+sizeOut ){
//     std::size_t s1=DSGetSize( *bufp );
//     int objId=DSGetObjectID( DSGetHeaderID( *bufp ) );
// #if 0    
//     std::cout << funcname << ": Size=" << std::setw(4) << size
//               << " objId=" << std::setw(3) << objId << std::endl;
// #endif
//     if( objId==SksObjTrackHit && objIdOut==SksObjTrackHitOutC && TrOut_){
//       TrackHit *hit = new TrackHit( bufp, TrIn_, TrOut_ );
//       if( hit ) hitContOut.push_back(hit);
//       else {
// 	std::cerr << funcname << ": new fail objId=" 
//                   << std::setw(3) << objId << std::endl;
//         ret=false;
//       }
//     }
//     else{
//       std::ostringstream mess;
//       mess << funcname << ": invalid object type " 
//                 << std::setw(3) << objId;
//       std::cerr << mess.str() << std::endl;
//       ret=false;
//       throw std::invalid_argument(mess.str());
//     }
//     bufp+=s1;
//   }

//   return ret;
// }
