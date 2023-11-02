/*
  DCHit.cc
*/

#include "DCHit.hh"
#include <cmath>

const double Deg2Rad = acos(-1)/180.;
const double Rad2Deg = 180./acos(-1);

DCHit::DCHit()
  : layer_(-1), wire_(-1), tdc_(-1),
    belongTrack_(false)
{
}

DCHit::DCHit( int layer, int wire, int tdc )
  : layer_(layer), wire_(wire), tdc_(tdc),
    belongTrack_(false)
{
}

DCHit::~DCHit()
{
}

double DCHit::GetResidual( void ) const
{
  double xlcal=xcal_*cos(angle_*Deg2Rad)+ycal_*sin(angle_*Deg2Rad);
  return xl_-xlcal;
}

