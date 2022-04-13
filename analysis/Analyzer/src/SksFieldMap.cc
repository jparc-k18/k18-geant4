/*
  SksFieldMap.cc
*/

#include "SksFieldMap.hh"

#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>

const double Deg2Rad = acos(-1.)/180.;
const double Rad2Deg = 180./acos(-1.);

SksFieldMap::SksFieldMap( const char *filename )
  : filename_(filename), Nx(0), Ny(0), Nz(0)
{
}

SksFieldMap::~SksFieldMap()
{
  cleanupMap();
}

bool SksFieldMap::Initialize( void )
{
  static const std::string funcname = "SksFieldMap::Initialize";

  std::ifstream fsin(filename_.c_str());

  if(!fsin){
    std::cerr << "[" << funcname << "]: file open fail" << std::endl;
    exit(-1);
  }
  cleanupMap();

  if( !(fsin >> Nx >> Ny >> Nz >> X0 >> Y0 >> Z0 >> dX >> dY >> dZ ) ){
    std::cerr << "[" << funcname << "]: Invalid format " << std::endl;
    exit(-1);
  }

  B.resize(Nx);
  for( int ix=0; ix<Nx; ++ix ){
    B[ix].resize(Ny);
    for( int iy=0; iy<Ny; ++iy ){
      B[ix][iy].resize(Nz);
    }
  }

  double x,y,z,bx,by,bz;
  int npoint=0;
  //double scale = 0.924;//QM4 to sksmap395a
  //double scale = 0.9298;//QM4 to QM4Full
  //double scale = 0.6;
  double scale = 1.;
  std::cout << "Scale of Field Map is " << scale << ", OK? (press return)"<<std::endl;;
  //getchar();
  std::cout << "Now reading Field Map " <<std::endl;;
  while( fsin ){
    if (npoint%50000 == 0) {
      std::cout << " (-_-)p[Wait]q ";
      fflush( stdout );
    }
    npoint++;

    fsin >> x >> y >> z >> bx >> by >> bz;
    int ix = int((x-X0+0.1*dX)/dX);
    int iy = int((y-Y0+0.1*dY)/dY);
    int iz = int((z-Z0+0.1*dZ)/dZ);
    if( ix>=0 && ix<Nx && iy>=0 && iy<Ny && iz>=0 && iz<Nz ){
//       B[ix][iy][iz].x=bx;
//       B[ix][iy][iz].y=by;
//       B[ix][iy][iz].z=bz;
      //       B[ix][iy][iz].x=0.88888888*bx;
      //       B[ix][iy][iz].y=0.88888888*by;
      //       B[ix][iy][iz].z=0.88888888*bz;
            B[ix][iy][iz].x = scale*bx;
            B[ix][iy][iz].y = scale*by;
            B[ix][iy][iz].z = scale*bz;
    }
  }

#if 1
  std::cout << "[" << funcname << "]: Finish Initialization " 
	    << std::endl;
#endif

  return true;
}

bool SksFieldMap::GetFieldValue( const double pointCM[3],
				 double *BfieldTesla ) const
{
  static const std::string funcname = "SksFieldMap::GetFieldValue";
  double xt=pointCM[0], yt=pointCM[1], zt=pointCM[2];

  //std::cout<<"[GetFieldValue] pos("<<xt<<","<<yt<<","<<zt<<")"<<std::endl;

  int ix1, ix2, iy1, iy2, iz1, iz2;
  ix1=int( (xt-X0)/dX );
  iy1=int( (yt-Y0)/dY );
  iz1=int( (zt-Z0)/dZ );

  // std::cout<<"[GetFieldValue] id1("<<ix1<<","<<iy1<<","<<iz1<<")"<<std::endl;
//   std::cout<<"[GetFieldValue] B("<<B[ix1][iy1][iz1].x<<","<<B[ix1][iy1][iz1].y<<","<<B[ix1][iy1][iz1].z<<")"<<std::endl;

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

  BfieldTesla[0]=bx; BfieldTesla[1]=by; BfieldTesla[2]=bz; 
  return true;
}

void SksFieldMap::cleanupMap( void )
{
  for( int ix=0; ix<Nx; ++ix ){
    for( int iy=0; iy<Ny; ++iy ){
      B[ix][iy].clear();
    }
    B[ix].clear();
  }
  B.clear();
}
