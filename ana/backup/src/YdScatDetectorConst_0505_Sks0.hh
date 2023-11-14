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
const double TargetVessThickness = 0.25; // mm

const int    FiberLayerNum = 10;
const int    FiberSegMax   = 1200;

const double OffsetCylindricalFromTgt = 100.;  //mm
const double zFiberFromTgt = 150; //mm
const double rMinCylindricalDetCont = 25.;   //mm 
const double rMaxCylindricalDetCont = 160.;   //mm 
const double zCylindricalDetCont = 520.;   //mm 

const double OffsetFiberFromCly = -50.; // mm
const double rMinFiberDetCont = 50.;    // mm
const double rMaxFiberDetCont = 90.;    // mm
const double zFiberDetCont = 400.;      // mm

const double rFiberLayer1 = 50.;   //mm 
const double startAngle1 = 0.;        // degree
const double durationAngle1 = 360.;   // degree

const double FiberSize1_X  = 0.5;   //mm 
const double FiberSize1_Y  = 0.5;   //mm 
const double FiberSize1_Z  = 400.0;   //mm 
//const double FiberAngle1   = 2.*atan((FiberSize1_XY/2.)/(rFiberLayer1-FiberSize1_XY/2.))*Rad2Deg; // degree
const double FiberAngle1   = 2.*atan((FiberSize1_X/2.)/(rFiberLayer1-FiberSize1_Y/2.))*Rad2Deg; // degree
//const int    segNumFiber1 = (int)(durationAngle1/FiberAngle1);
//Calculate segNumFiber1 using the above equation
const int    segNumFiber1 = 625;

const double rFiberLayer2 = 51.;   //mm 
const double startAngle2 = 0.;        // degree
const double durationAngle2 = 360.;   // degree

const double FiberSize2_X  = 0.5;   //mm 
const double FiberSize2_Y  = 0.5;   //mm 
const double FiberSize2_Z  = 400.0;   //mm 
//const int    segNumFiber2 = (int)(zCylindricalDetCont/FiberSize2_XY);
//Calculate segNumFiber1 using the above equation
const double tiltAngleFiber2 = atan(2.*Pi*rFiberLayer2/FiberSize2_Z)*Rad2Deg;
const int    segNumFiber2 = 500;

const double rFiberLayer3 = 60.;   //mm 
const double startAngle3 = 0.;        // degree
const double durationAngle3 = 360.;   // degree

const double FiberSize3_X  = 0.5;   //mm 
const double FiberSize3_Y  = 0.5;   //mm 
const double FiberSize3_Z  = 400.0;   //mm 
const double FiberAngle3   = 2.*atan((FiberSize3_X/2.)/(rFiberLayer3-FiberSize3_Y/2.))*Rad2Deg; // degree
//const int    segNumFiber3 = (int)(durationAngle3/FiberAngle3);
//Calculate segNumFiber3 using the above equation
const int    segNumFiber3 = 750;

const double rFiberLayer4 = 61.;   //mm 
const double startAngle4 = 0.;        // degree
const double durationAngle4 = 360.;   // degree

const double FiberSize4_X  = 0.5;   //mm 
const double FiberSize4_Y  = 0.5;   //mm 
const double FiberSize4_Z  = 400.0;   //mm 
//const int    segNumFiber4 = (int)(zCylindricalDetCont/FiberSize4_XY);
//Calculate segNumFiber1 using the above equation
const double tiltAngleFiber4 = -atan(2.*Pi*rFiberLayer4/FiberSize4_Z)*Rad2Deg;
const int    segNumFiber4 = 552;

const double rFiberLayer5 = 70.;   //mm 
const double startAngle5 = 0.;        // degree
const double durationAngle5 = 360.;   // degree

const double FiberSize5_X  = 0.5;   //mm 
const double FiberSize5_Y  = 0.5;   //mm 
const double FiberSize5_Z  = 400.0;   //mm 
const double FiberAngle5   = 2.*atan((FiberSize5_X/2.)/(rFiberLayer5-FiberSize5_Y/2.))*Rad2Deg; // degree
//const int    segNumFiber5 = (int)(durationAngle5/FiberAngle5);
//Calculate segNumFiber5 using the above equation
const int    segNumFiber5 = 876;

const double rFiberLayer6 = 71.;   //mm 
const double startAngle6 = 0.;        // degree
const double durationAngle6 = 360.;   // degree

const double FiberSize6_X  = 0.5;   //mm 
const double FiberSize6_Y  = 0.5;   //mm 
const double FiberSize6_Z  = 400.0;   //mm 
//const int    segNumFiber6 = (int)(zCylindricalDetCont/FiberSize6_XY);
//Calculate segNumFiber6 using the above equation
const double tiltAngleFiber6 = atan(2.*Pi*rFiberLayer6/FiberSize6_Z)*Rad2Deg;
const int    segNumFiber6 = 596;

const double rFiberLayer7 = 80.; //mm 
const double startAngle7 = 0.;        // degree
const double durationAngle7 = 360.;   // degree

const double FiberSize7_X  = 0.5;  //mm
const double FiberSize7_Y  = 0.5;  //mm
const double FiberSize7_Z  = 400.0; //mm
const double FiberAngle7   = 2.*atan((FiberSize7_X/2.)/(rFiberLayer7-FiberSize7_Y/2.))*Rad2Deg; // degree
//const int    segNumFiber7 = (int)(durationAngle7/FiberAngle7);
//Calculate segNumFiber7 using the above equation
const int    segNumFiber7 = 1002;

const double rFiberLayer8 = 81.; //mm
const double startAngle8 = 0.;        // degree
const double durationAngle8 = 360.;   // degree

const double FiberSize8_X  = 0.5; //mm
const double FiberSize8_Y  = 0.5; //mm
const double FiberSize8_Z  = 400.0; //mm
//const int    segNumFiber8 = (int)(zCylindricalDetCont/FiberSize8_XY);
//Calculate segNumFiber8 using the above equation
const double tiltAngleFiber8 = -atan(2.*Pi*rFiberLayer8/FiberSize8_Z)*Rad2Deg;
const int    segNumFiber8 = 628;

const double rFiberLayer9 = 90.;  //mm
const double startAngle9 = 0.;        // degree
const double durationAngle9 = 360.;   // degree

const double FiberSize9_X  = 0.5; //mm
const double FiberSize9_Y  = 0.5; //mm
const double FiberSize9_Z  = 400.0;
const double FiberAngle9   = 2.*atan((FiberSize9_X/2.)/(rFiberLayer9-FiberSize9_Y/2.))*Rad2Deg; // degree
//const int    segNumFiber9 = (int)(durationAngle9/FiberAngle9);
//Calculate segNumFiber9 using the above equation
const int    segNumFiber9 = 1127;

const double rFiberLayer10 = 91.; //mm
const double startAngle10 = 0.;        // degree
const double durationAngle10 = 360.;   // degree

const double FiberSize10_X  = 0.5; //mm
const double FiberSize10_Y  = 0.5; //mm
const double FiberSize10_Z  = 400.0; //mm
//const int    segNumFiber10 = (int)(zCylindricalDetCont/FiberSize10_XY);
//Calculate segNumFiber10 using the above equation
const double tiltAngleFiber10 = atan(2.*Pi*rFiberLayer10/FiberSize10_Z)*Rad2Deg;
const int    segNumFiber10 = 654;

const int segNumFiber[FiberLayerNum+1] = {0, segNumFiber1, segNumFiber2,
					  segNumFiber3, segNumFiber4,
					  segNumFiber5, segNumFiber6,
					  segNumFiber7, segNumFiber8,
					  segNumFiber9, segNumFiber10};

const double rFiber[FiberLayerNum+1] = {0., rFiberLayer1, rFiberLayer2,
					  rFiberLayer3, rFiberLayer4,
					  rFiberLayer5, rFiberLayer6,
					  rFiberLayer7, rFiberLayer8,
					  rFiberLayer9, rFiberLayer10};

const double tiltAngleFiber[FiberLayerNum+1] =   {0., 
						    0., tiltAngleFiber2,
						    0., tiltAngleFiber4,
						    0., tiltAngleFiber6,
						    0., tiltAngleFiber8,
						    0., tiltAngleFiber10};

const double fiberStartAngle[FiberLayerNum+1] = {0,
						 startAngle1, startAngle2,
						 startAngle3, startAngle4,
						 startAngle5, startAngle6,
						 startAngle7, startAngle8,
						 startAngle9, startAngle10};

const double fiberDurationAngle[FiberLayerNum+1] = {0,
						 durationAngle1, durationAngle2,
						 durationAngle3, durationAngle4,
						 durationAngle5, durationAngle6,
						 durationAngle7, durationAngle8,
						 durationAngle9, durationAngle10};

const double FiberAngle[FiberLayerNum+1] = {0.,
					    FiberAngle1, 0.,
					    FiberAngle3, 0.,
					    FiberAngle5, 0.,
					    FiberAngle7, 0.,
					    FiberAngle9, 0.};

const double FiberSizeX[FiberLayerNum+1] = { 0.,
					     FiberSize1_X, FiberSize2_X,
					     FiberSize3_X, FiberSize4_X,
					     FiberSize5_X, FiberSize6_X,
					     FiberSize7_X, FiberSize8_X,
					     FiberSize9_X, FiberSize10_X};

const double FiberSizeY[FiberLayerNum+1] = { 0.,
					     FiberSize1_Y, FiberSize2_Y,
					     FiberSize3_Y, FiberSize4_Y,
					     FiberSize5_Y, FiberSize6_Y,
					     FiberSize7_Y, FiberSize8_Y,
					     FiberSize9_Y, FiberSize10_Y};

const double FiberSizeZ[FiberLayerNum+1] = { 0.,
					     FiberSize1_Z, FiberSize2_Z,
					     FiberSize3_Z, FiberSize4_Z,
					     FiberSize5_Z, FiberSize6_Z,
					     FiberSize7_Z, FiberSize8_Z,
					     FiberSize9_Z, FiberSize10_Z};

const int CrystalLayerNum = 4;
// Crystal, ForwardCrystal, Rc, ForwardRc

const double OffsetCrystalFromCly = -50.; //mm
const double rCrystal = 100.;   //mm 
const double startAngleCrystal = 0.;        // degree
const double durationAngleCrystal = 360.;   // degree

const double CrystalSize_X = 31.6; //mm
const double CrystalSize_Y = 25.0; //mm
const double CrystalSize_Z  = 400.0; //mm

const double CrystalAngle   = 2.*atan((CrystalSize_X/2.)/rCrystal)*Rad2Deg; // degree
const double CrystalSize_X2 = 2.*(rCrystal+CrystalSize_Y)*tan(CrystalAngle/2.*Deg2Rad);
//const int    segNumCrystal = (int)(durationAngleCrystal/CrystalAngle)+1;
//Calculate segNumFiber5 using the above equation
const int    segNumCrystal = 20;


const double OffsetForwardCrystalFromCly = 190.; //mm
const double rForwardCrystal = 70; //mm
const double startAngleForwardCrystal = 0.;        // degree
const double durationAngleForwardCrystal = 360.;   // degree

const double ForwardCrystalSize_Y = 55.0; //mm
const double ForwardCrystalSize_Z  = 80.0; //mm
const double ForwardCrystalSize_X = 2.*(rForwardCrystal)*tan(CrystalAngle/2.*Deg2Rad);
const double ForwardCrystalSize_X2 = 2.*(rForwardCrystal+ForwardCrystalSize_Y)*tan(CrystalAngle/2.*Deg2Rad);
//const int    segNumCrystal = (int)(durationAngleCrystal/CrystalAngle)+1;
//Calculate segNumFiber5 using the above equation
const int    segNumForwardCrystal = 20;


const double OffsetRcFromCly = 0.; //mm
const double rRc = 130.; //mm
const double startAngleRc = 0.;        // degree
const double durationAngleRc = 360.;   // degree

const double RcAngle = CrystalAngle;

const double RcSize_Y = 20.0; //mm
const double RcSize_Z  = 500.0; //mm
const double RcSize_X = 2.*(rRc)*tan(RcAngle/2.*Deg2Rad);
const double RcSize_X2 = 2.*(rRc+RcSize_Y)*tan(RcAngle/2.*Deg2Rad);
const int    segNumRc = 20;

const double OffsetForwardRcFromCly = 240.; //mm
const double rForwardRc = 70.; //mm
const double startAngleForwardRc = 0.;        // degree
const double durationAngleForwardRc = 360.;   // degree

const double ForwardRcSize_Y = 60.0; //mm
const double ForwardRcSize_Z  = 20.0; //mm
const double ForwardRcSize_X = 2.*(rForwardRc)*tan(RcAngle/2.*Deg2Rad);
const double ForwardRcSize_X2 = 2.*(rForwardRc+ForwardRcSize_Y)*tan(RcAngle/2.*Deg2Rad);
const int    segNumForwardRc = 20;

const int segNumCalori[CrystalLayerNum+1] = {0,
					 segNumCrystal, segNumForwardCrystal,
					 segNumRc, segNumForwardRc};

const double rCaloriCenter[CrystalLayerNum+1] = {0,
						 rCrystal+CrystalSize_Y/2.,
						 rForwardCrystal+ForwardCrystalSize_Y/2.,
						 rRc+RcSize_Y/2.,
						 rForwardRc+ForwardRcSize_Y/2.};

const double CaloriAngle[CrystalLayerNum+1] = {0,
					       CrystalAngle, CrystalAngle,
					       RcAngle, RcAngle};

#endif
