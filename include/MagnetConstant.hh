#ifndef MagnetConstant_h
#define MagnetConstant_h 1

#include "globals.hh"
#include "G4SystemOfUnits.hh"

//par to make area
//const double AreaX = 0;//1.05*m;
//const double AreaY = 0;//-1.5*m;
//const double AreaAngle = -15.;

const double slope = 0.2*mm;

// Q1 magnet
const double Q1x = 0.5*m;
const double Q1y = 0.5*m;
const double Q1z = 0.88*m;
const double Q1Yx = 2.4*m;
const double Q1Yy = 2.4*m;
const double Q1Yz = 0.88*m;
const double Q1zEff = 1.024*m;
const double Q1B0 = -13.95; // [kG]
//const double Q1B0 = -12.555; // [kG] 90% of kato val
const double Q1a0 = 0.155*m;    // [m]
const double Q1cy = Q1x*sqrt(2.)/2. - Q1a0;//radius

// Q2 magnet
const double Q2Yx = 2.0*m;
const double Q2Yy = 2.0*m;
const double Q2Yz = 0.54*m;
const double Q2zEff = 0.716*m;
//const double Q2B0 = -10.924; // [kG]
const double Q2B0 = 9.92444; // [kG]
//const double Q2B0 = 9.837; //[kG] target-Q1 46 cm
//const double Q2a0 = 0.2186*m;      // [m]
const double Q2a0 = 0.18*m;// Corrected, 16Mar2015
const double Q2x = 0.8*m;//2.25*Q2a0*2;
const double Q2y = Q2a0*2;//1.25*Q2a0*2;
//const double Q2z = 0.54*m;
const double Q2z = 0.5*m; // Corrected, 16Mar2015
const double Q2a = Q2x/2 - 0.35*Q2a0;
const double Q2b = Q2y/2 - 0.5*Q2a0;

// D magnet
const double rhoD = 3.000*m; // [m] 
const double bendAngleD = 70.000; // [degree] 
const double DSpaceX = 0.8*m;  // [m]
const double Dr1 = rhoD - DSpaceX/2.;
const double Dr2 = rhoD + DSpaceX/2.;
//const double DHalfGap = 0.1572*m;
const double DHalfGap = 0.16*m; // Corrected, 16Mar2015
const double DLength = 0;//3.6651914*m;
const double DB0 = -13.5;//-14.44; // [kG]
const double alphaD = 0; // [degree]
const double betaD = 0; // [degree]
// D magnet surface
const double rhoDf = 3.000*m; // [m] 
const double bendAngleDf = 70.000; // [degree] 
const double DfSpaceX = 2.4*m;//1.824*m;  // [m]
const double Dfr1 = rhoDf - DfSpaceX/2.;
const double Dfr2 = rhoDf + 0.8*m + 0.35*m;
const double DfHalfGap = 0.6*m;
const double DfLength = 3.66519*m;
const double DfB0 = 0; // [kG]
const double alphaDf = 0; // [degree]
const double betaDf = 0; // [degree]

// DRIFT0 VI-Q1
const double driftL0 = 0.6*m;
// DRIFT1 Q1-Q2
const double driftL1 = 0.4*m;
// DRIFT1 Q2-D
const double driftL2 = 0.4*m;
// DRIFT1 Q13-VO
const double driftL5 = 0.270*m;

#endif
