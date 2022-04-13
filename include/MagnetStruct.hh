#ifndef __MAGPAR_H__
#define __MAGPAR_H__

#include "magnetcommon.h"

struct Quad 
{
  double B0;  /* Magnetic Field at Pole [kGauss] */
  double a0;  /* aperture [m] */
  double L;   /* Length of Magnet */
  double EffL; /* Effective Length*/
  double Pos[XYZ]; /* Position of Center of Magnet */
};

struct Bend
{
  double B0;  /* Magnetic Field [kGauss] */
  double X;   /* half length of X [m] */
  double Y;   /* half length of Y [m] */
  double L;   /* Length of Magnet */
  double EffL; /* Effective Length*/
  double Pos[XYZ]; /* Position of Center of Magnet */
};

struct Sect
{
  double B0;  /* Magnetic Field [kGauss] */
  double Rho; /* bending radius */
  double Width; /* Full Width of X [m] */
  double Y;   /* half length of Y [m] */
  double L;   /* Length of Magnet */
  double EffL; /* Effective Length*/
  double Pos[XYZ]; /* Position of rotation Center */
  double BendAngle; /* BendingAngle */
  double alpha; /* angle between central track and entrance pole */
  double beta;  /* angle between central track and exit pole */
};

struct Magnet {
  Quad  quad;
  Bend  bend;
  Sect  sect;
  int   ID;       /* ID represents the kind of the magnet */
  char  name[50]; /* name of magnet */
  double theta;   /* slope against beamline */
};

#endif
