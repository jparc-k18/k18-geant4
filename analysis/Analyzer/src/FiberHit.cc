/*
  FiberHit.cc
*/

#include "FiberHit.hh"
#include <cmath>

const double Deg2Rad = acos(-1)/180.;
const double Rad2Deg = 180./acos(-1);

FiberHit::FiberHit()
  : layer_(-1), segment_(-1), x_(-999.9), y_(-999.9), z_(-999.9),
    edep_(-999.9), time_(-999.9), pid_(-1)
{
}

FiberHit::FiberHit( int layer, int segment, double edep, double time, int pid )
  : layer_(layer), segment_(segment), edep_(edep), time_(time), pid_(pid),
    x_(-999.9), y_(-999.9), z_(-999.9)
{
}

FiberHit::~FiberHit()
{
}

