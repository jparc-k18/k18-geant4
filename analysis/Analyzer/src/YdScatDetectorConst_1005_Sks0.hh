#ifndef YdScatDetectorConst_H
#define YdScatDetectorConst_H 1

#include "RadDeg.hh"

const double LH2TgtX = 20.;   //mm 
const double LH2TgtY = 300.;   //mm 
const double LH2TgtZ = 40.;   //mm 

const double LD2TgtX = 20.;   //mm 
const double LD2TgtY = 300.;   //mm 
const double LD2TgtZ = 40.;   //mm 

// Slice target
const double SLH2TgtX = 40.;   //mm 
const double SLH2TgtY = 20.;   //mm 
const double SLH2TgtZ = 20.;   //mm 

const double SLD2TgtX = 40.;   //mm 
//const double SLD2TgtY = 20.;
const double SLD2TgtY = 30.;   //mm 
const double SLD2TgtZ = 20.;   //mm 

const double STgtSpace = 5.;   //mm 
const int    STgtLayNum = 5;

// Cylindrical target
const double CyLH2TgtR = 40.; //mm //tyokkei
const double CyLH2TgtZ = 300.; // mm
const double TargetVessThickness = 1.5; // mm

// Cylindrical target
const double VaccumChamWinR = 60.; //mm //tyokkei
const double VaccumChamWinZ = 400.; //mm
const double VaccumChamThickness = 1.5; //mm;

const int    FiberLayerNum = 4;
const int    FiberSegMax   = 600;

const double OffsetCylindricalFromTgt = 140.;  //mm
const double rMinCylindricalDetCont = 40.;   //mm 
const double rMaxCylindricalDetCont = 120.;   //mm 
const double zCylindricalDetCont = 400.;   //mm 

const double OffsetFiberFromCly = -90.; // mm
const double rMinFiberDetCont = 50.;    // mm
const double rMaxFiberDetCont = 90.;    // mm
const double zFiberDetCont = 400.;      // mm

const double rFiberLayer1 = 60.;   //mm 
const double startAngle1 = 0.;        // degree
const double durationAngle1 = 360.;   // degree

const double FiberSize1_XY = 1.0;   //mm 
const double FiberSize1_X  = 0.5;   //mm 
const double FiberSize1_Y  = 1.0;   //mm 
const double FiberSize1_Z  = 400.0;   //mm 
//const double FiberAngle1   = 2.*atan((FiberSize1_XY/2.)/(rFiberLayer1-FiberSize1_XY/2.))*Rad2Deg; // degree
const double FiberAngle1   = 2.*atan((FiberSize1_Y/2.)/(rFiberLayer1-FiberSize1_X/2.))*Rad2Deg; // degree
//const int    segNumFiber1 = (int)(durationAngle1/FiberAngle1);
//Calculate segNumFiber1 using the above equation
const int    segNumFiber1 = 375;

const double rFiberLayer2 = 61.;   //mm 
const double startAngle2 = 0.;        // degree
const double durationAngle2 = 360.;   // degree

const double FiberSize2_XY = 1.0;   //mm 
const double FiberSize2_X  = 0.5;   //mm 
const double FiberSize2_Y  = 1.0;   //mm 
const double FiberSize2_Z  = 400.0;   //mm 
//const int    segNumFiber2 = (int)(zCylindricalDetCont/FiberSize2_XY);
//Calculate segNumFiber1 using the above equation
const int    segNumFiber2 = 400;

const double rFiberLayer3 = 70.;   //mm 
const double startAngle3 = 0.;        // degree
const double durationAngle3 = 360.;   // degree

const double FiberSize3_XY = 1.0;   //mm 
const double FiberSize3_X  = 0.5;   //mm 
const double FiberSize3_Y  = 1.0;   //mm 
const double FiberSize3_Z  = 400.0;   //mm 
const double FiberAngle3   = 2.*atan((FiberSize3_Y/2.)/(rFiberLayer3-FiberSize3_X/2.))*Rad2Deg; // degree
//const int    segNumFiber3 = (int)(durationAngle3/FiberAngle3);
//Calculate segNumFiber3 using the above equation
const int    segNumFiber3 = 438;

const double rFiberLayer4 = 71.;   //mm 
const double startAngle4 = 0.;        // degree
const double durationAngle4 = 360.;   // degree

const double FiberSize4_XY = 1.0;   //mm 
const double FiberSize4_X  = 0.5;   //mm 
const double FiberSize4_Y  = 1.0;   //mm 
const double FiberSize4_Z  = 400.0;   //mm 
//const int    segNumFiber4 = (int)(zCylindricalDetCont/FiberSize4_XY);
//Calculate segNumFiber1 using the above equation
const int    segNumFiber4 = 400;

const double rFiberLayer5 = 90.;   //mm 
const double startAngle5 = 0.;        // degree
const double durationAngle5 = 360.;   // degree

const double FiberSize5_XY = 1.0;   //mm 
const double FiberSize5_X  = 0.5;   //mm 
const double FiberSize5_Y  = 1.0;   //mm 
const double FiberSize5_Z  = 400.0;   //mm 
const double FiberAngle5   = 2.*atan((FiberSize5_Y/2.)/(rFiberLayer5-FiberSize5_X/2.))*Rad2Deg; // degree
//const int    segNumFiber5 = (int)(durationAngle5/FiberAngle5);
//Calculate segNumFiber5 using the above equation
const int    segNumFiber5 = 563;

const double rFiberLayer6 = 91.;   //mm 
const double startAngle6 = 0.;        // degree
const double durationAngle6 = 360.;   // degree

const double FiberSize6_XY = 1.0;   //mm 
const double FiberSize6_X  = 0.5;   //mm 
const double FiberSize6_Y  = 1.0;   //mm 
const double FiberSize6_Z  = 400.0;   //mm 
//const int    segNumFiber6 = (int)(zCylindricalDetCont/FiberSize6_XY);
//Calculate segNumFiber6 using the above equation
const int    segNumFiber6 = 400;


const int    CdcLayerNum = 6;
const double CdcZ = 580.0;      //mm
const double CdcRMin = 100.0;   //mm
const double CdcRMax = 280.0;   //mm
const double rCdcLayer1 = 115.; //mm
const double rCdcLayer2 = 145.; //mm
const double rCdcLayer3 = 175.; //mm
const double rCdcLayer4 = 205.; //mm
const double rCdcLayer5 = 235.; //mm
const double rCdcLayer6 = 265.; //mm


const double rCrystal = 300.;   //mm 
const double startAngleCrystal = 0.;        // degree
const double durationAngleCrystal = 360.;   // degree

const double CrystalSize_X = 54.0;   //mm 
const double CrystalSize_Y = 60.0;   //mm 
const double CrystalSize_Z  = 580.0;   //mm 

const double CrystalAngle   = 2.*atan((CrystalSize_X/2.)/(rCrystal-CrystalSize_Y/2.))*Rad2Deg; // degree
const double CrystalSize_X2 = 2.*(rCrystal+CrystalSize_Y/2.)*tan(CrystalAngle/2.*Deg2Rad);
//const int    segNumCrystal = (int)(durationAngleCrystal/CrystalAngle)+1;
//Calculate segNumFiber5 using the above equation
const int    segNumCrystal = 32;

#endif
