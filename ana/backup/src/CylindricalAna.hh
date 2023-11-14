#ifndef CylindricalAna_h
#define CylindricalAna_h 1

#include "ThreeVector.hh"
#include "YdScatDetectorConst.hh"

#include <vector>

struct FiberCluster
{
  int     clusterSize;
  double  clusterCenter;
  double  centerPos[3];
  double  edep;
  double  residual;
  double crossPosIn[3];
  double crossPosOut[3];
  double pathLength;
  double  edepNorm;
};


struct CrystalCluster
{
  int     clusterSize;
  double  clusterCenter;
  double  centerPos[3];
  double  edep;
};

class CylindricalAna
{

private:
  int          layerNumForScatPart_;
  int          layerNumForScatPartXY_;
  int          layerNumForScatPartUV_;
  FiberCluster fclForScatPart_[FiberLayerNum+1];
  CrystalCluster cclForScatPart_;
  CrystalCluster crcForScatPart_;
  int          xyFitFlag_;
  double       chisqrXY_;
  double       chisqrZ_;
  double       Axy_, Bxy_;
  double       Az_, Bz_;
  bool         statusXY_;
  bool         statusZ_;
  ThreeVector  ScatVec_, ScatPos_;
  bool         weightByEdep;
  std::vector <int> hitIndexXY_, hitIndexUV_;

public:
  CylindricalAna();
  ~CylindricalAna();
  void SetFiberDataForScatPart(int layer, int seg, double de);
  void SetCrystalDataForScatPart( int seg, double de);
  void SetRcDataForScatPart( int seg, double de);
  int  Evaluate();
  bool FitXYPlane();
  bool FitZPlane();
  ThreeVector GetScatVector() { return ScatVec_;};
  ThreeVector GetScatPosition() { return ScatPos_;};
  double GetChisqrXY() { return chisqrXY_;};
  double GetChisqrZ() { return chisqrZ_;};
  inline double GetResidual(int i);
  inline double GetFiberClusterEdep(int i);
  inline double GetFiberClusterEdepNorm(int i);
  inline double GetFiberPathLength(int i);
  double GetCrystalClusterEdep() { return cclForScatPart_.edep;};
  double GetRcClusterEdep() { return crcForScatPart_.edep;};
  bool CalcCrossPoint();
  double calcPhi(double x, double y);
  double calcZPos(int layer, double x, double y);
  bool   GetHitPosition(int layer, double *pos);
};

inline double CylindricalAna::GetResidual(int i)
{
  if (i>=1 && i<=FiberLayerNum)
    return fclForScatPart_[i].residual;
  else {
    std::cout << "CylindricalAna::GetResidual  Invalid layerId " << i << std::endl;
    return fclForScatPart_[0].residual;
  }
}

inline double CylindricalAna::GetFiberClusterEdep(int i)
{
  if (i>=1 && i<=FiberLayerNum)
    return fclForScatPart_[i].edep;
  else {
    std::cout << "CylindricalAna::GetFiberClusterEdep  Invalid layerId " << i << std::endl;
    return fclForScatPart_[0].edep;
  }
}

inline double CylindricalAna::GetFiberClusterEdepNorm(int i)
{
  if (i>=1 && i<=FiberLayerNum)
    return fclForScatPart_[i].edepNorm;
  else {
    std::cout << "CylindricalAna::GetFiberClusterEdepNorm  Invalid layerId " << i << std::endl;
    return fclForScatPart_[0].edepNorm;
  }
}

inline double CylindricalAna::GetFiberPathLength(int i)
{
  if (i>=1 && i<=FiberLayerNum)
    return fclForScatPart_[i].pathLength;
  else {
    std::cout << "CylindricalAna::GetFiberPathLength  Invalid layerId " << i << std::endl;
    return fclForScatPart_[0].pathLength;
  }
}

inline bool CylindricalAna::GetHitPosition(int layer, double *pos)
{
  if (layer>=1 && layer<=FiberLayerNum) {
    if (fclForScatPart_[layer].centerPos[2]>-990.) {
      pos[0] = fclForScatPart_[layer].centerPos[0];
      pos[1] = fclForScatPart_[layer].centerPos[1];
      pos[2] = fclForScatPart_[layer].centerPos[2];
      return true;
    } else {
      return false;
    }
  } else {
    return false;
  }
}

#endif
