/*
  "S2SFieldMap.cc"
  
  Modified by Toshi Gogami , 21Nov2014
*/
#include "MagnetConstant.hh"
#include "S2SFieldMap.hh"

#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <cstdlib>

const double Deg2Rad = acos(-1.)/180.;
const double Rad2Deg = 180./acos(-1.);

S2SFieldMap::S2SFieldMap( const char *filename, double ScaleFactorQ1, double ScaleFactorQ2, double ScaleFactor)
  : filename_(filename),
    ScaleFactorQ1_(ScaleFactorQ1),
    ScaleFactorQ2_(ScaleFactorQ2),
    ScaleFactor_(ScaleFactor),
    Nx(0),
    Ny(0),
    Nz(0)
{
}

S2SFieldMap::~S2SFieldMap()
{
  cleanupMap();
}

bool S2SFieldMap::Initialize( void )
{
  static const std::string funcname = "S2SFieldMap::Initialize";

  std::ifstream fsin( filename_.c_str() );
  
  if(!fsin){
    std::cerr << "[" << funcname << "]: file open fail : " << filename_ << std::endl;
    std::exit(-1);
  }
  cleanupMap();

  if( !(fsin >> Ny >> Nz >> Nx >>
	Y0 >> Z0 >> X0 >> dY >> dZ >> dX ) ){
    std::cerr << "[" << funcname << "]: Invalid format " << std::endl;
    std::exit(-1);
  }
  dX = dX*10.; // [cm] -> [mm]
  dY = dY*10.;
  dZ = dZ*10.;
  X0 = X0*10.; // [cm] -> [mm]
  Y0 = Y0*10.;
  Z0 = Z0*10.;

  B.resize(Nx);
  for( int ix=0; ix<Nx; ++ix ){
    B[ix].resize(Ny);
    for( int iy=0; iy<Ny; ++iy ){
      B[ix][iy].resize(Nz);
    }
  }
  

  double xlim_Q1 = -(rhoD*tan(bendAngleD/2.*Deg2Rad) + driftL2 + Q2z + driftL1/2.);
  double xlim_Q2 = -(rhoD*tan(bendAngleD/2.*Deg2Rad) + driftL2/2.);
  
  double x,y,z,bx,by,bz;  
  int npoint=0;
  int readn = 0;
  std::cout << "Now reading Field Map " <<std::endl;
  while( fsin ){
    if (npoint%50000 == 0) {
      //std::cout << " (-_-)p[Wait]q ";
      if(readn<5)std::cout << "(๑ˇεˇ)  ";
      else if (readn==5)std::cout << "(๑ˇεˇ)  " << std::endl;
      else if (readn<10)std::cout << "ヾ(｡>﹏<｡)ﾉﾞ  " ;
      else if (readn==10)std::cout << "ヾ(｡>﹏<｡)ﾉﾞ" << std::endl;
      else if (readn<15)std::cout << "(๑ˇεˇ)  ";
      else if (readn==15)std::cout << "(๑ˇεˇ) " << std::endl;
      else if (readn<20)std::cout << "ヾ(｡>﹏<｡)ﾉﾞ  " ;
      else if (readn==20)std::cout << "ヾ(｡>﹏<｡)ﾉﾞ" << std::endl;
      else std::cout << "." ;
      fflush( stdout );
      readn++;
    }
    npoint++;

    fsin >> y >> z >> x >> by >> bz >> bx;
    x = x*10.; // [cm] -> [mm]
    y = y*10.;
    z = z*10.;
    int ix = int((x-X0+0.1*dX)/dX);
    int iy = int((y-Y0+0.1*dY)/dY);
    int iz = int((z-Z0+0.1*dZ)/dZ);
    if( ix>=0 && ix<Nx && iy>=0 && iy<Ny && iz>=0 && iz<Nz ){
      if(x<xlim_Q1){
	B[ix][iy][iz].x = bx*ScaleFactorQ1_;
	B[ix][iy][iz].y = by*ScaleFactorQ1_;
	B[ix][iy][iz].z = bz*ScaleFactorQ1_;
      }
      else if (x<xlim_Q2){
	B[ix][iy][iz].x = bx*ScaleFactorQ2_;
	B[ix][iy][iz].y = by*ScaleFactorQ2_;
	B[ix][iy][iz].z = bz*ScaleFactorQ2_;
      }
      else{
	B[ix][iy][iz].x = bx*ScaleFactor_;
	B[ix][iy][iz].y = by*ScaleFactor_;
	B[ix][iy][iz].z = bz*ScaleFactor_;
      }
    }
  }
  std::cout << std::endl << "Finished reading Field Map " << std::endl;;
  return true;
}

bool S2SFieldMap::GetFieldValue( const double point[3],
				 double *Bfield ) const
{
  static const std::string funcname = "S2SFieldMap::GetFieldValue";
  double xt=point[0], yt=point[1], zt=point[2];

  int ix1, ix2, iy1, iy2, iz1, iz2;
  ix1=int( (xt-X0)/dX );
  iy1=int( (yt-Y0)/dY );
  iz1=int( (zt-Z0)/dZ );

  double wx1, wx2, wy1, wy2, wz1, wz2;
  if( ix1<0 ) { ix1=ix2=0; wx1=1.; wx2=0.; }
  else if( ix1>=Nx-1 ) { ix1=ix2=Nx-1; wx1=1.; wx2=0.; }
  else { ix2=ix1+1; wx1=(X0+dX*ix2-xt)/dX; wx2=1.-wx1; }

  if( iy1<0 ) { iy1=iy2=0; wy1=1.; wy2=0.; }
  else if( iy1>=Ny-1 ) { iy1=iy2=Ny-1; wy1=1.; wy2=0.; }
  else { iy2=iy1+1; wy1=(Y0+dY*iy2-yt)/dY; wy2=1.-wy1; }

  if( iz1<0 ) { iz1=iz2=0; wz1=1.; wz2=0.; }
  else if( iz1>=Nz-1 ) { iz1=iz2=Nz-1; wz1=1.; wz2=0.; }
  else { iz2=iz1+1; wz1=(Z0+dZ*iz2-zt)/dZ; wz2=1.-wz1; }

  double bx1=wx1*wy1*B[ix1][iy1][iz1].x+wx1*wy2*B[ix1][iy2][iz1].x
    +wx2*wy1*B[ix2][iy1][iz1].x+wx2*wy2*B[ix2][iy2][iz1].x;
  double bx2=wx1*wy1*B[ix1][iy1][iz2].x+wx1*wy2*B[ix1][iy2][iz2].x
    +wx2*wy1*B[ix2][iy1][iz2].x+wx2*wy2*B[ix2][iy2][iz2].x;
  double bx=wz1*bx1+wz2*bx2;
  double by1=wx1*wy1*B[ix1][iy1][iz1].y+wx1*wy2*B[ix1][iy2][iz1].y
    +wx2*wy1*B[ix2][iy1][iz1].y+wx2*wy2*B[ix2][iy2][iz1].y;
  double by2=wx1*wy1*B[ix1][iy1][iz2].y+wx1*wy2*B[ix1][iy2][iz2].y
    +wx2*wy1*B[ix2][iy1][iz2].y+wx2*wy2*B[ix2][iy2][iz2].y;
  double by=wz1*by1+wz2*by2;
  double bz1=wx1*wy1*B[ix1][iy1][iz1].z+wx1*wy2*B[ix1][iy2][iz1].z
    +wx2*wy1*B[ix2][iy1][iz1].z+wx2*wy2*B[ix2][iy2][iz1].z;
  double bz2=wx1*wy1*B[ix1][iy1][iz2].z+wx1*wy2*B[ix1][iy2][iz2].z
    +wx2*wy1*B[ix2][iy1][iz2].z+wx2*wy2*B[ix2][iy2][iz2].z;
  double bz=wz1*bz1+wz2*bz2;

  Bfield[0]=bx; Bfield[1]=by; Bfield[2]=bz; 
  return true;
}

void S2SFieldMap::cleanupMap( void )
{
  for( int ix=0; ix<Nx; ++ix ){
    for( int iy=0; iy<Ny; ++iy ){
      B[ix][iy].clear();
    }
    B[ix].clear();
  }
  B.clear();
}
