#include "CylindricalAna.hh"
#include "YdScatDetectorConst.hh"

#include <string>
#include <cmath>
#include <iostream>

//#include "RadDeg.hh"

CylindricalAna::CylindricalAna()
  : layerNumForScatPart_(0), layerNumForScatPartXY_(0), layerNumForScatPartUV_(0),
    xyFitFlag_(-1), statusXY_(false), statusZ_(false), chisqrXY_(-1.), chisqrZ_(-1.),
    weightByEdep(true)
{
  for (int i=0; i<FiberLayerNum+1; i++) {
    fclForScatPart_[i].clusterSize = 0;
    fclForScatPart_[i].clusterCenter = 0.;
    fclForScatPart_[i].centerPos[0] = -999.;
    fclForScatPart_[i].centerPos[1] = -999.;
    fclForScatPart_[i].centerPos[2] = -999.;
    fclForScatPart_[i].edep = 0.;
    fclForScatPart_[i].residual = -999.;
    fclForScatPart_[i].crossPosIn[0] = -999.;
    fclForScatPart_[i].crossPosIn[1] = -999.;
    fclForScatPart_[i].crossPosIn[2] = -999.;
    fclForScatPart_[i].crossPosOut[0] = -999.;
    fclForScatPart_[i].crossPosOut[1] = -999.;
    fclForScatPart_[i].crossPosOut[2] = -999.;
    fclForScatPart_[i].pathLength = -999.;
    fclForScatPart_[i].edepNorm = -999.;
  }
  cclForScatPart_.clusterSize = 0;
  cclForScatPart_.clusterCenter = 0.;
  cclForScatPart_.centerPos[0] = -999.;
  cclForScatPart_.centerPos[1] = -999.;
  cclForScatPart_.centerPos[2] = -999.;
  cclForScatPart_.edep = 0.;

  crcForScatPart_.clusterSize = 0;
  crcForScatPart_.clusterCenter = 0.;
  crcForScatPart_.centerPos[0] = -999.;
  crcForScatPart_.centerPos[1] = -999.;
  crcForScatPart_.centerPos[2] = -999.;
  crcForScatPart_.edep = 0.;

}

CylindricalAna::~CylindricalAna()
{}

void CylindricalAna::SetFiberDataForScatPart(int layer, int seg, double de)
{
  if (weightByEdep) {
    fclForScatPart_[layer].clusterCenter = 
      (fclForScatPart_[layer].clusterCenter*fclForScatPart_[layer].edep +
       de*(double)seg)/(fclForScatPart_[layer].edep+de);
    fclForScatPart_[layer].clusterSize += 1;
    fclForScatPart_[layer].edep += de;

  } else {
    fclForScatPart_[layer].clusterCenter = 
      (fclForScatPart_[layer].clusterCenter*fclForScatPart_[layer].clusterSize +
       (double)seg)/(double)(fclForScatPart_[layer].clusterSize+1);
    fclForScatPart_[layer].clusterSize += 1;
    fclForScatPart_[layer].edep += de;
  }
}

void CylindricalAna::SetCrystalDataForScatPart( int seg, double de)
{
  cclForScatPart_.clusterCenter = 
    (cclForScatPart_.clusterCenter*cclForScatPart_.clusterSize + seg)/(double)(cclForScatPart_.clusterSize+1);
  cclForScatPart_.clusterSize += 1;
  cclForScatPart_.edep += de;
}

void CylindricalAna::SetRcDataForScatPart( int seg, double de)
{
  crcForScatPart_.clusterCenter = 
    (crcForScatPart_.clusterCenter*crcForScatPart_.clusterSize + seg)/(double)(crcForScatPart_.clusterSize+1);
  crcForScatPart_.clusterSize += 1;
  crcForScatPart_.edep += de;
}


int CylindricalAna::Evaluate()
{
  for (int i=1; i<=FiberLayerNum; i++) {
    if (fclForScatPart_[i].clusterSize >= 1) {
      layerNumForScatPart_++;

      double seg = fclForScatPart_[i].clusterCenter;

      if (i==1 || i==3 || i==5 || i==7 || i==9) {
	layerNumForScatPartXY_++;
	double phi = (fiberStartAngle[i]+FiberAngle[i]*seg); // Degree
	fclForScatPart_[i].centerPos[0] = rFiber[i]*cos(phi*Deg2Rad);
	fclForScatPart_[i].centerPos[1] = rFiber[i]*sin(phi*Deg2Rad);
	fclForScatPart_[i].centerPos[2] = 0.;
	hitIndexXY_.push_back(i);
      } else if (i==2 || i==4 || i==6 || i==8 || i==10) {
	hitIndexUV_.push_back(i);
	layerNumForScatPartUV_++;
      }
    }
  }

  if (cclForScatPart_.clusterSize >= 1) {
    double seg = cclForScatPart_.clusterCenter;
    double phi = startAngleCrystal+CrystalAngle*seg;

    cclForScatPart_.centerPos[0] = rCrystal*cos(phi*Deg2Rad);
    cclForScatPart_.centerPos[1] = rCrystal*sin(phi*Deg2Rad);
    cclForScatPart_.centerPos[2] = 0.;

  }

  if (crcForScatPart_.clusterSize >= 1) {
    double seg = crcForScatPart_.clusterCenter;
    double phi = startAngleRc+CrystalAngle*seg;

    crcForScatPart_.centerPos[0] = rRc*cos(phi*Deg2Rad);
    crcForScatPart_.centerPos[1] = rRc*sin(phi*Deg2Rad);
    crcForScatPart_.centerPos[2] = 0.;

  }

  return layerNumForScatPart_;

}

bool CylindricalAna::FitXYPlane()
{
  static const std::string funcname = "CylindricalAna::FitXYPlane ";


  bool fiberFlag = false;
  if (layerNumForScatPartXY_ >= 2)
    fiberFlag = true;

  if (!fiberFlag) {
    std::cerr << funcname << " Hit Number is not enough," 
	      << " L1:" << fclForScatPart_[1].clusterSize 
	      << " L3:" << fclForScatPart_[3].clusterSize 
	      << " L5:" << fclForScatPart_[5].clusterSize 
	      << " L7:" << fclForScatPart_[7].clusterSize 
	      << " L9:" << fclForScatPart_[9].clusterSize 
	      << std::endl;
    statusXY_ = false;
    return statusXY_;
  }

  double slopeFiber;
  slopeFiber = (fclForScatPart_[hitIndexXY_[0]].centerPos[1]-fclForScatPart_[hitIndexXY_[layerNumForScatPartXY_-1]].centerPos[1])/
    (fclForScatPart_[hitIndexXY_[0]].centerPos[0]-fclForScatPart_[hitIndexXY_[layerNumForScatPartXY_-1]].centerPos[0]);

  std::vector <double> x, y, s;
  for (int i=0; i<layerNumForScatPartXY_; i++) {
    x.push_back(fclForScatPart_[hitIndexXY_[i]].centerPos[0]);
    y.push_back(fclForScatPart_[hitIndexXY_[i]].centerPos[1]);
    s.push_back(FiberSizeX[hitIndexXY_[i]]/sqrt(12.));
  }

  double A=0., B=0., C=0., D=0., E=0., F=0.;
  if (fabs(slopeFiber) <= 1.) {
    // y = ax + b
    for (int i=0; i<layerNumForScatPartXY_; i++) {
      A += x[i]/(s[i]*s[i]);
      B += 1./(s[i]*s[i]);
      C += y[i]/(s[i]*s[i]);
      D += x[i]*x[i]/(s[i]*s[i]);
      E += x[i]*y[i]/(s[i]*s[i]);
      F += y[i]*y[i]/(s[i]*s[i]);
    }

    xyFitFlag_=0;

    Axy_=(E*B-C*A)/(D*B-A*A);
    Bxy_=(D*C-E*A)/(D*B-A*A);

    chisqrXY_ = 0.;
    for (int i=0; i<layerNumForScatPartXY_; i++) {
      double ycal;
      ycal = Axy_*x[i] + Bxy_;
      fclForScatPart_[hitIndexXY_[i]].residual = (y[i]-ycal);
      chisqrXY_ += (y[i]-ycal)*(y[i]-ycal)/(s[i]*s[i]);
    }
    chisqrXY_ /= layerNumForScatPartXY_;
  } else {
    // x = ay + b
    for (int i=0; i<layerNumForScatPartXY_; i++) {
      A += y[i]/(s[i]*s[i]);
      B += 1./(s[i]*s[i]);
      C += x[i]/(s[i]*s[i]);
      D += y[i]*y[i]/(s[i]*s[i]);
      E += x[i]*y[i]/(s[i]*s[i]);
      F += x[i]*x[i]/(s[i]*s[i]);
    }

    xyFitFlag_=1;

    Axy_=(E*B-C*A)/(D*B-A*A);
    Bxy_=(D*C-E*A)/(D*B-A*A);

    chisqrXY_ = 0.;
    for (int i=0; i<layerNumForScatPartXY_; i++) {
      double xcal;
      xcal = Axy_*y[i] + Bxy_;
      fclForScatPart_[hitIndexXY_[i]].residual = (x[i]-xcal);
      chisqrXY_ += (x[i]-xcal)*(x[i]-xcal)/(s[i]*s[i]);
    }
    chisqrXY_ /= layerNumForScatPartXY_;
  }

  statusXY_ = true;
  return statusXY_;
}


bool CylindricalAna::FitZPlane()
{
  static const std::string funcname = "CylindricalAna::FitZPlane ";

  if (!statusXY_) {
    std::cerr << funcname << "statusXY_ is " << statusXY_;
    statusZ_=false;
    return statusZ_;
  }

  bool fiberFlag=false;
  bool cdcFlag=false;

  if (layerNumForScatPartUV_>=2)
    fiberFlag = true;

  if (!fiberFlag) {
    std::cerr << funcname << " Hit Number is not enough," 
	      << " L2:" << fclForScatPart_[2].clusterSize 
	      << " L4:" << fclForScatPart_[4].clusterSize 
	      << " L6:" << fclForScatPart_[6].clusterSize 
	      << " L8:" << fclForScatPart_[8].clusterSize 
	      << " L10:" << fclForScatPart_[10].clusterSize 
	      << std::endl;
    statusZ_ = false;
    return statusZ_;
  }

  double A=0., B=0., C=0., D=0., E=0., F=0.;

  if (xyFitFlag_ == 0) {
    double xtmp1,xtmp2;
    double a=Axy_;
    double b=Bxy_;
    double r;
    double c;

    double x0, xcal, ycal, zcal;
    std::vector <double> x, z, s;
    for (int i=0; i<layerNumForScatPartUV_; i++) {
      r = rFiber[hitIndexUV_[i]];
      c = (a*b)*(a*b)-(a*a+1.)*(b*b-r*r);
      if (c < 0.) {
	std::cerr << funcname << "Cannot find the answer of 2d pol at Layer"
		  << hitIndexUV_[i] << std::endl;
	statusZ_=false;
	return statusZ_;
      }
      xtmp1 = (-a*b+sqrt(c))/(a*a+1.);
      xtmp2 = (-a*b-sqrt(c))/(a*a+1.);

      if (fclForScatPart_[hitIndexUV_[i]-1].clusterSize>0)
	x0 = fclForScatPart_[hitIndexUV_[i]-1].centerPos[0];
      else
	x0 = fclForScatPart_[hitIndexXY_[0]].centerPos[0];

      if (fabs(xtmp1-x0) < fabs(xtmp2-x0)) {
	xcal = xtmp1;
      } else {
	xcal = xtmp2;
      }
      fclForScatPart_[hitIndexUV_[i]].centerPos[0] = xcal;
      fclForScatPart_[hitIndexUV_[i]].centerPos[1] = Axy_*xcal+Bxy_;
      ycal = Axy_*xcal+Bxy_;
      // calculate z position;
      zcal = calcZPos(hitIndexUV_[i], xcal, ycal);
      fclForScatPart_[hitIndexUV_[i]].centerPos[2] = zcal;

      x.push_back(xcal);
      z.push_back(zcal);
      s.push_back(FiberSizeX[hitIndexUV_[i]]/sqrt(12.));
    }

    for (int i=0; i<layerNumForScatPartUV_; i++) {
      A += x[i]/(s[i]*s[i]);
      B += 1./(s[i]*s[i]);
      C += z[i]/(s[i]*s[i]);
      D += x[i]*x[i]/(s[i]*s[i]);
      E += x[i]*z[i]/(s[i]*s[i]);
      F += z[i]*z[i]/(s[i]*s[i]);
    }

    Az_=(E*B-C*A)/(D*B-A*A);
    Bz_=(D*C-E*A)/(D*B-A*A);

    ScatVec_ = ThreeVector(1., Axy_, Az_);
    ScatPos_ = ThreeVector(0., Bxy_, Bz_);

    if (z[0] < z[layerNumForScatPartUV_-1]) {
      if (ScatVec_.z() < 0.)
	ScatVec_*=-1.;
    } else {
      if (ScatVec_.z() > 0.)
	ScatVec_*=-1.;
    }

    chisqrZ_ = 0.;
    for (int i=0; i<layerNumForScatPartUV_; i++) {
      double zcal2;
      zcal2 = Az_*x[i] + Bz_;
      fclForScatPart_[hitIndexUV_[i]].residual = (z[i]-zcal2);
      chisqrZ_ += (z[i]-zcal2)*(z[i]-zcal2)/(s[i]*s[i]);
    }
    chisqrZ_ /= layerNumForScatPartUV_;

    statusZ_=true;

  } else if (xyFitFlag_ == 1) {
    double ytmp1,ytmp2;
    double a=Axy_;
    double b=Bxy_;
    double r;
    double c;

    double y0, xcal, ycal, zcal;
    std::vector <double> y, z, s;
    for (int i=0; i<layerNumForScatPartUV_; i++) {
      r = rFiber[hitIndexUV_[i]];
      c = (a*b)*(a*b)-(a*a+1.)*(b*b-r*r);
      if (c < 0.) {
	std::cerr << funcname << "Cannot find the answer of 2d pol at Layer" 
		  << hitIndexUV_[i] << std::endl;
	statusZ_=false;
	return statusZ_;
      }
      ytmp1 = (-a*b+sqrt(c))/(a*a+1.);
      ytmp2 = (-a*b-sqrt(c))/(a*a+1.);

      if (fclForScatPart_[hitIndexUV_[i]-1].clusterSize>0)
	y0 = fclForScatPart_[hitIndexUV_[i]-1].centerPos[1];
      else
	y0 = fclForScatPart_[hitIndexXY_[0]].centerPos[1];

      if (fabs(ytmp1-y0) < fabs(ytmp2-y0)) {
	ycal = ytmp1;
      } else {
	ycal = ytmp2;
      }
      fclForScatPart_[hitIndexUV_[i]].centerPos[1] = ycal;
      fclForScatPart_[hitIndexUV_[i]].centerPos[0] = Axy_*ycal+Bxy_;
      xcal = Axy_*ycal+Bxy_;

      // calculate z position;
      zcal = calcZPos(hitIndexUV_[i], xcal, ycal);
      fclForScatPart_[hitIndexUV_[i]].centerPos[2] = zcal;

      y.push_back(ycal);
      z.push_back(zcal);
      s.push_back(FiberSizeX[hitIndexUV_[i]]/sqrt(12.));
    }

    for (int i=0; i<layerNumForScatPartUV_; i++) {
      A += y[i]/(s[i]*s[i]);
      B += 1./(s[i]*s[i]);
      C += z[i]/(s[i]*s[i]);
      D += y[i]*y[i]/(s[i]*s[i]);
      E += y[i]*z[i]/(s[i]*s[i]);
      F += z[i]*z[i]/(s[i]*s[i]);
    }

    Az_=(E*B-C*A)/(D*B-A*A);
    Bz_=(D*C-E*A)/(D*B-A*A);
  
    ScatVec_ = ThreeVector(Axy_, 1., Az_);
    ScatPos_ = ThreeVector(Bxy_, 0., Bz_);

    if (z[0] < z[layerNumForScatPartUV_-1]) {
      if (ScatVec_.z() < 0.)
	ScatVec_*=-1.;
    } else {
      if (ScatVec_.z() > 0.)
	ScatVec_*=-1.;
    }

    chisqrZ_ = 0.;
    for (int i=0; i<layerNumForScatPartUV_; i++) {
      double zcal2;
      zcal2 = Az_*y[i] + Bz_;
      fclForScatPart_[hitIndexUV_[i]].residual = (z[i]-zcal2);
      chisqrZ_ += (z[i]-zcal2)*(z[i]-zcal2)/(s[i]*s[i]);
    }
    chisqrZ_ /= layerNumForScatPartUV_;

    statusZ_=true;
  }


  CalcCrossPoint();

  return statusZ_;
  
}

bool CylindricalAna::CalcCrossPoint()
{
  static const std::string funcname = "CylindricalAna::FitZPlane ";

  for (int i=1; i<=FiberLayerNum; i++) {
    if (xyFitFlag_ == 0 && fclForScatPart_[i].clusterSize>0) {
      double x0 = fclForScatPart_[i].centerPos[0];
      double xtmp1,xtmp2;
      double a=Axy_;
      double b=Bxy_;
      double r;
      double c;

      double xIn, yIn, zIn;
      double xOut, yOut, zOut;

      // Inner surface
      r = rFiber[i]-FiberSizeX[i]/2.;
      c = (a*b)*(a*b)-(a*a+1.)*(b*b-r*r);
      if (c < 0.) {
	std::cerr << funcname << "Cannot find the answer of 2d pol at Layer"
		  << i  << std::endl;
	return false;
      }
      xtmp1 = (-a*b+sqrt(c))/(a*a+1.);
      xtmp2 = (-a*b-sqrt(c))/(a*a+1.);
      if (fabs(xtmp1-x0) < fabs(xtmp2-x0)) {
	xIn = xtmp1;
	yIn = Axy_*xIn+Bxy_;
	zIn = Az_*xIn + Bz_;
      } else {
	xIn = xtmp2;
	yIn = Axy_*xIn+Bxy_;
	zIn = Az_*xIn + Bz_;
      }

      // Outer surface
      r = rFiber[i]+FiberSizeX[i]/2.;
      c = (a*b)*(a*b)-(a*a+1.)*(b*b-r*r);
      if (c < 0.) {
	std::cerr << funcname << "Cannot find the answer of 2d pol at Layer"
		  << i << std::endl;
	return false;
      }
      xtmp1 = (-a*b+sqrt(c))/(a*a+1.);
      xtmp2 = (-a*b-sqrt(c))/(a*a+1.);
      if (fabs(xtmp1-x0) < fabs(xtmp2-x0)) {
	xOut = xtmp1;
	yOut = Axy_*xOut+Bxy_;
	zOut = Az_*xOut + Bz_;
      } else {
	xOut = xtmp2;
	yOut = Axy_*xOut+Bxy_;
	zOut = Az_*xOut + Bz_;
      }

      fclForScatPart_[i].crossPosIn[0] = xIn;
      fclForScatPart_[i].crossPosIn[1] = yIn;
      fclForScatPart_[i].crossPosIn[2] = zIn;
      fclForScatPart_[i].crossPosOut[0] = xOut;
      fclForScatPart_[i].crossPosOut[1] = yOut;
      fclForScatPart_[i].crossPosOut[2] = zOut;
      fclForScatPart_[i].pathLength = 
	sqrt((xIn - xOut)*(xIn - xOut) + (yIn - yOut)*(yIn - yOut) + 
	     (zIn - zOut)*(zIn - zOut));
      fclForScatPart_[i].edepNorm = 
	fclForScatPart_[i].edep/fclForScatPart_[i].pathLength;
	
    } else if (xyFitFlag_ == 1 && fclForScatPart_[i].clusterSize>0) {
      double y0 = fclForScatPart_[i].centerPos[1];
      double ytmp1,ytmp2;
      double a=Axy_;
      double b=Bxy_;
      double r;
      double c;

      double xIn, yIn, zIn;      
      double xOut, yOut, zOut;

      // Inner surface
      r = rFiber[i]-FiberSizeX[i]/2.;
      c = (a*b)*(a*b)-(a*a+1.)*(b*b-r*r);
      if (c < 0.) {
	std::cerr << funcname << "Cannot find the answer of 2d pol at Layer"
		  << i << std::endl;
	return false;
      }
      ytmp1 = (-a*b+sqrt(c))/(a*a+1.);
      ytmp2 = (-a*b-sqrt(c))/(a*a+1.);
      if (fabs(ytmp1-y0) < fabs(ytmp2-y0)) {
	yIn = ytmp1;
	xIn = Axy_*yIn+Bxy_;
	zIn = Az_*yIn+Bz_;
      } else {
	yIn = ytmp2;
	xIn = Axy_*yIn+Bxy_;
	zIn = Az_*yIn+Bz_;
      }
      // Outer surface
      r = rFiber[i]+FiberSizeX[i]/2.;
      c = (a*b)*(a*b)-(a*a+1.)*(b*b-r*r);
      if (c < 0.) {
	std::cerr << funcname << "Cannot find the answer of 2d pol at Layer"
		  << i << std::endl;
	return false;
      }
      ytmp1 = (-a*b+sqrt(c))/(a*a+1.);
      ytmp2 = (-a*b-sqrt(c))/(a*a+1.);
      if (fabs(ytmp1-y0) < fabs(ytmp2-y0)) {
	yOut = ytmp1;
	xOut = Axy_*yOut+Bxy_;
	zOut = Az_*yOut+Bz_;
      } else {
	yOut = ytmp2;
	xOut = Axy_*yOut+Bxy_;
	zOut = Az_*yOut+Bz_;
      }

      fclForScatPart_[i].crossPosIn[0] = xIn;
      fclForScatPart_[i].crossPosIn[1] = yIn;
      fclForScatPart_[i].crossPosIn[2] = zIn;
      fclForScatPart_[i].crossPosOut[0] = xOut;
      fclForScatPart_[i].crossPosOut[1] = yOut;
      fclForScatPart_[i].crossPosOut[2] = zOut;
      fclForScatPart_[i].pathLength = 
	sqrt((xIn - xOut)*(xIn - xOut) + (yIn - yOut)*(yIn - yOut) + 
	     (zIn - zOut)*(zIn - zOut));
      fclForScatPart_[i].edepNorm = 
	fclForScatPart_[i].edep/fclForScatPart_[i].pathLength;

    }
    
  }

  return true;
}

double CylindricalAna::calcPhi(double x, double y)
{
  if (x>=0. && y>=0)
    return atan(y/x);
  else if (x<0. && y>=0. )
    return Pi + atan(y/x);
  else if (x<0. && y<0. )
    return Pi + atan(y/x);
  else if (x>=0. && y<0. )
    return 2.*Pi + atan(y/x);
}

double CylindricalAna::calcZPos(int layer, double x, double y)
{
  double phi = calcPhi(x, y); // radian
  double r = rFiber[layer]*phi;
  // local position
  double lx = fclForScatPart_[layer].clusterCenter*FiberSizeX[layer];

  double theta = tiltAngleFiber[layer];
  double x0 = lx*cos(theta*Deg2Rad);
  double z0 = lx*sin(theta*Deg2Rad);

  double z =  -1./tan(theta*Deg2Rad)*(r-x0)+z0;
  if (z<0)
    z = z + FiberSizeZ[layer];


  /*
  std::cout << "layer:" << layer << " seg:" 
	    << fclForScatPart_[layer].clusterCenter 
	    << " phi:" << phi << " r:" << r
	    << " " << x << " " << y << " " << z
	    << std::endl;
  */


  return z - zFiberFromTgt;

}
