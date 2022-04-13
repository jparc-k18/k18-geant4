/*
  "S2SFieldMap.hh"
  
  Modified by Toshi Gogami , 21Nov2014
*/

#ifndef S2SFieldMap_h
#define S2SFieldMap_h 1

#include <string>
#include <vector>

class S2SFieldMap
{
public:
  S2SFieldMap( const char *filename=0, double ScaleFactor=1.0 );
  ~S2SFieldMap();

private:
  S2SFieldMap( const S2SFieldMap & );
  S2SFieldMap & operator = ( const S2SFieldMap & );
  
private:
  std::string filename_;
  double ScaleFactor_;
  
public:
  bool Initialize( void );
  bool GetFieldValue( const double point[3], double *Bfield ) const;
private:
  struct FD {
    float x, y, z;
  };
  
  typedef std::vector < std::vector < std::vector < FD > > > FDContainer;
  FDContainer B;
  double X0, Y0, Z0, dX, dY, dZ;
  int Nx, Ny, Nz;
  
  void cleanupMap( void );
};

#endif
