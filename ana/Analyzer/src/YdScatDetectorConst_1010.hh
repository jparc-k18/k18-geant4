#ifndef YdScatDetectorConst_H
#define YdScatDetectorConst_H 1

#include "RadDeg.hh"

const double LH2TgtX = 20.;    //mm
const double LH2TgtY = 300.;   //mm
const double LH2TgtZ = 40.;    //mm

const double LD2TgtX = 20.;    //mm 
const double LD2TgtY = 300.;   //mm 
const double LD2TgtZ = 40.;    //mm 

// Slice target
const double SLH2TgtX = 40.;   //mm 
const double SLH2TgtY = 20.;	 //mm 
const double SLH2TgtZ = 20.;	 //mm 

const double SLD2TgtX = 40.;   //mm 
//const double SLD2TgtY = 2.*cm;
const double SLD2TgtY = 30.;   //mm 
const double SLD2TgtZ = 20.;   //mm 

const double STgtSpace = 5;    //mm 
const int    STgtLayNum = 5;

const int    FiberLayerNum = 6;
const int    FiberSegMax   = 500;

const double rMinCylindricalDetCont = 40.;     //mm 
const double rMaxCylindricalDetCont = 120.;    //mm 
const double zCylindricalDetCont = 400.;       //mm 

const double rFiberLayer1 = 60.;        //mm 
const double startAngle1 = 45.;        // degree
const double durationAngle1 = 270.;   // degree

const double FiberSize1_XY = 1.0;     //mm 
const double FiberSize1_Z  = 400.0;   //mm 
const double FiberAngle1   = 2.*atan((FiberSize1_XY/2.)/(rFiberLayer1-FiberSize1_XY/2.))*Rad2Deg; // degree
//const int    segNumFiber1 = (int)(durationAngle1/FiberAngle1);
//Calculate segNumFiber1 using the above equation
const int    segNumFiber1 = 280;

const double rFiberLayer2 = 61.;       //mm 
const double startAngle2 = 45.;       // degree
const double durationAngle2 = 270.;   // degree

const double FiberSize2_XY = 1.0;     //mm 
const double FiberSize2_Z  = 400.0;   //mm 
//const int    segNumFiber2 = (int)(zCylindricalDetCont/FiberSize2_XY);
//Calculate segNumFiber1 using the above equation
const int    segNumFiber2 = 400;

const double rFiberLayer3 = 70.;      //mm 
const double startAngle3 = 45.;      // degree
const double durationAngle3 = 270.;   // degree

const double FiberSize3_XY = 1.0;     //mm 
const double FiberSize3_Z  = 400.0;   //mm 
const double FiberAngle3   = 2.*atan((FiberSize3_XY/2.)/(rFiberLayer3-FiberSize3_XY/2.))*Rad2Deg; // degree
//const int    segNumFiber3 = (int)(durationAngle3/FiberAngle3);
//Calculate segNumFiber3 using the above equation
const int    segNumFiber3 = 327;

const double rFiberLayer4 = 71.;      //mm 
const double startAngle4 = 45.;      // degree
const double durationAngle4 = 270.;   // degree

const double FiberSize4_XY = 1.0;     //mm 
const double FiberSize4_Z  = 400.0;   //mm 
//const int    segNumFiber4 = (int)(zCylindricalDetCont/FiberSize4_XY);
//Calculate segNumFiber1 using the above equation
const int    segNumFiber4 = 400;

const double rFiberLayer5 = 90.;      //mm  
const double startAngle5 = 45.;      // degree
const double durationAngle5 = 270.;   // degree

const double FiberSize5_XY = 1.0;     //mm 
const double FiberSize5_Z  = 400.0;   //mm 
const double FiberAngle5   = 2.*atan((FiberSize5_XY/2.)/(rFiberLayer5-FiberSize5_XY/2.))*Rad2Deg; // degree
//const int    segNumFiber5 = (int)(durationAngle5/FiberAngle5);
//Calculate segNumFiber5 using the above equation
const int    segNumFiber5 = 421;

const double rFiberLayer6 = 91.;      //mm  
const double startAngle6 = 45.;      // degree
const double durationAngle6 = 270.;   // degree

const double FiberSize6_XY = 1.0;     //mm 
const double FiberSize6_Z  = 400.0;   //mm 
//const int    segNumFiber6 = (int)(zCylindricalDetCont/FiberSize6_XY);
//Calculate segNumFiber6 using the above equation
const int    segNumFiber6 = 400;

const double rCrystal = 110.;         //mm 
const double startAngleCrystal = 45.;        // degree
const double durationAngleCrystal = 270.;   // degree

const double CrystalSize_X = 20.0;      //mm 
const double CrystalSize_Y = 20.0;      //mm 
const double CrystalSize_Z  = 400.0;    //mm 

const double CrystalAngle   = 2.*atan((CrystalSize_X/2.)/(rCrystal-CrystalSize_Y/2.))*Rad2Deg; // degree
const double CrystalSize_X2 = 2.*(rCrystal+CrystalSize_Y/2.)*tan(CrystalAngle/2.*Deg2Rad);
//const int    segNumCrystal = (int)(durationAngleCrystal/CrystalAngle)+1;
//Calculate segNumFiber5 using the above equation
const int    segNumCrystal = 24;

#endif
