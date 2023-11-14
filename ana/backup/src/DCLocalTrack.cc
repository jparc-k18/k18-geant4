/*
  DCLocalTrack.cc
*/

#include "DCLocalTrack.hh"
#include "DCHit.hh"
#include "DCConfMan.hh"

#include "MathTools.hh"

#include <string>
#include <vector>
#include <cmath>

const double Deg2Rad = acos(-1.)/180.;
const double Rad2Deg = 180./acos(-1.);

const int ReservedNumOfHits = 100;
const int DCLocalMinNHits = 6;

DCLocalTrack::DCLocalTrack()
  : status_(false), x0_(0.0), y0_(0.0), u0_(0.0), v0_(0.0)
{
  hitArray.reserve( ReservedNumOfHits );
}

DCLocalTrack::~DCLocalTrack()
{
}

DCHit *DCLocalTrack::GetHit( std::size_t nth ) const
{
  if( nth<hitArray.size() )
    return hitArray[nth];
  else
    return 0;
}

DCHit *DCLocalTrack::GetHitOfLayerNumber( int lnum ) const
{
  for( std::size_t i=0; i<hitArray.size(); ++i )
    if( hitArray[i]->GetLayer()==lnum )
      return hitArray[i];
  return 0;
}

bool DCLocalTrack::DoFit( void )
{
  const std::string funcname = "[DCLocalTrack::DoFit()]";

  const DCConfMan & confMan=DCConfMan::GetInstance();

  std::size_t n = hitArray.size();
  if( n < DCLocalMinNHits ) return status_ = false;

  //std::cout<<"hitArray="<< n<<std::endl;

  std::vector <double> z, w, s, ct, st;
  z.reserve(n); w.reserve(n); w.reserve(n);
  ct.reserve(n); st.reserve(n);

  for( std::size_t i=0; i<n; ++i ){
    DCHit *hitp = hitArray[i];
    if( hitp ){
      int lnum = hitp->GetLayer();
      //if(lnum==131) std::cout<<"131 Local"<<std::endl;
      double ww = confMan.GetResolution( lnum );
      double zz = confMan.GetLocalZ( lnum );
      double aa = hitp->GetTiltAngle()*Deg2Rad;
      double LocalPos;
      LocalPos =  hitp->GetLocalHitPos();

      z.push_back( zz ); w.push_back( 1./(ww*ww) ); 
      s.push_back( LocalPos );
      ct.push_back( cos(aa) ); st.push_back( sin(aa) );
    }
  }
  std::size_t nn = z.size();

  double matrx[16], *mtp[4], fitp[4];
  mtp[0]=&matrx[0]; mtp[1]=&matrx[4]; mtp[2]=&matrx[8]; mtp[3]=&matrx[12];

  for( int i=0; i<4; ++i ){
    fitp[i]=0.0;
    for( int j=0; j<4; ++j ){
      mtp[i][j]=0.0;
    }
  }

  for( std::size_t i=0; i<nn; ++i ){
    double ww=w[i], zz=z[i], ss=s[i], ctt=ct[i], stt=st[i];
    mtp[0][0] += ww*ctt*ctt;
    mtp[0][1] += ww*zz*ctt*ctt;
    mtp[0][2] += ww*ctt*stt;
    mtp[0][3] += ww*zz*ctt*stt;
    mtp[1][1] += ww*zz*zz*ctt*ctt;
    mtp[1][2] += ww*zz*ctt*stt;
    mtp[1][3] += ww*zz*zz*ctt*stt;
    mtp[2][2] += ww*stt*stt;
    mtp[2][3] += ww*zz*stt*stt;
    mtp[3][3] += ww*zz*zz*stt*stt;

    fitp[0] += ww*ss*ctt;
    fitp[1] += ww*zz*ss*ctt;
    fitp[2] += ww*ss*stt;
    fitp[3] += ww*zz*ss*stt;
  }
  mtp[1][0]=mtp[0][1]; mtp[2][0]=mtp[0][2]; mtp[3][0]=mtp[0][3];
  mtp[2][1]=mtp[1][2]; mtp[3][1]=mtp[1][3]; mtp[3][2]=mtp[2][3];

  std::vector<int> indxc(nn), indxd(nn), ipiv(nn);

  if( MathTools::GaussJordan(mtp,4,fitp,&indxc[0],
			     &indxd[0],&ipiv[0])==false ){
    std::cerr << funcname << ": Fitting fails" << std::endl;
    return status_=false;
  }
  x0_=fitp[0]; y0_=fitp[2]; u0_=fitp[1]; v0_=fitp[3];

  double chisqr=0.0;
  for( std::size_t i=0; i<nn; ++i ){
    double ww=w[i], zz=z[i];
    double scal=GetX(zz)*ct[i]+GetY(zz)*st[i];
    chisqr += ww*(s[i]-scal)*(s[i]-scal);
  }
  chisqr /= nn-4.;
  chisqr_=chisqr;

  /*
  for( std::size_t i=0; i<n; ++i ){
    DCHit *hitp = hitArray[i];
    if( hitp ){
      int lnum = hitp->GetLayer();
      double zz = DCConfMan::GetInstance().GetLocalZ( lnum );
      hitp->SetCalPosition( GetX(zz), GetY(zz) );
    }
  }
  */

#if 0
  std::cout<< "*********************************************" <<std::endl;
  std::cout<< "x0="  << x0_
	   << " y0=" << y0_
	   << " u0=" << u0_
	   << " v0=" << v0_ << std::endl;
  std::cout<< "LocalChisqr=" << chisqr_ << std::endl;
  std::cout<< "*********************************************" <<std::endl;
#endif

  return status_=true;
}
