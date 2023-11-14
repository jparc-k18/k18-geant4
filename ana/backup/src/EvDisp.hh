#ifndef EvDisp_h
#define EvDisp_h 1

#include "TROOT.h"
#include "TApplication.h"
#include "TRint.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TGeometry.h"
#include "TMixture.h"
#include "TBRIK.h"
#include "TTRD1.h"
#include "TTRD2.h"
#include "TTUBS.h"
#include "TTUBS.h"
#include "TRotMatrix.h"
#include "TNode.h"
#include "TPolyMarker3D.h"
#include "TPolyLine3D.h"
#include "TView.h"
#include "TPad.h"
#include "TButton.h"
#include "TMarker3DBox.h"
#include "TPave.h"
#include "TPaveLabel.h"
#include "TPaveText.h"
#include "TGraph.h"
#include "TH2F.h"
#include "TLine.h"
#include "TArc.h"

#include "YdScatDetectorConst.hh"

#ifndef ROOT_TObject
#include <TObject.h>
#endif

#include "ThreeVector.hh"

extern  TApplication *theApp;

const int MaxParticleType = 15;
const int MaxHitEv = 100;

const int LayerNumUVFiber = 5;

class EvDisp : public TObject {
private:
  TPad      *tp_[3];
  TCanvas   *tc_;
  TH2F      *hbase[3];
  TArc      *tgtXY;

  TLine     *fl[FiberLayerNum+1][FiberSegMax][4];
  TLine     *cl[segNumCrystal][4];
  TLine     *fcl[segNumCrystal][4];
  TLine     *rl[segNumRc][4];
  int       hitFiberNum1_; 
  TArc      *farc1[MaxHitEv];
  int       hitCaloriNum1_; 
  TArc      *carc1[MaxHitEv];

  TLine     *tgtlXZ[4];
  TLine     *flXZ[FiberLayerNum+1][2]; // 2 is up and down
  TLine     *clXZ1[4]; 
  TLine     *clXZ2[4]; 
  TLine     *rclXZ1[4]; 
  TLine     *rclXZ2[4]; 

  TLine     *tgtlYZ[4];
  TLine     *flYZ[FiberLayerNum+1][2]; // 2 is up and down
  TLine     *clYZ1[4]; 
  TLine     *clYZ2[4]; 
  TLine     *fclYZ1[4]; 
  TLine     *fclYZ2[4]; 
  TLine     *rclYZ1[4]; 
  TLine     *rclYZ2[4]; 
  TLine     *frclYZ1[4]; 
  TLine     *frclYZ2[4]; 

  TPad      *tp2_[LayerNumUVFiber];
  TCanvas   *tc2_;
  TH2F      *hbase2[LayerNumUVFiber];
  TLine     *uvfl[LayerNumUVFiber][4];
  TLine     *huvfl[MaxHitEv][2];
  TLine     *hrfl[MaxHitEv];
  int       hitUVFiberNum_; 
  int       hitRFiberNum_; 

  TArc      *hitposXZ[MaxHitEv];
  TArc      *hitposYZ[MaxHitEv];
  int       hitPosXZNum_; 
  int       hitPosYZNum_; 

  TLine     *trlineYZ[MaxHitEv];
  TLine     *trlineXZ[MaxHitEv];
  TLine     *trlineXY[MaxHitEv];
  int       trLineNum_; 

public:
  EvDisp(void);
  ~EvDisp(void);

  void Initialize(void);
  void InitializeXZPlane(void);
  void InitializeYZPlane(void);
  void EndOfEvent();
  void DrawFiber(int layer, int seg, double x0, double y0, double theta);
  void DrawCalori(int seg, double x0, double y0, double thetaCal);
  void DrawCaloriForward(int seg, double x0, double y0, double thetaCal);
  void DrawRc(int seg, double x0, double y0, double thetaCal);
  void DrawUVFiber(void);
  void DrawBox(TLine *l[4], double p[4][2], int color);
  void calcRotate(double theta, double x, double y, double *rotpos);
  void DrawHitFiber(int layer, int seg, int pid, double de);
  void DrawHitCalori(int layer, int seg, int pid, double de);
  void DrawHitUVFiber(int layer, int seg, int pid, double de);
  void DrawHitRFiber(int layer, int seg, int pid, double de);
  void DrawHitXZPos(double *pos, int pid);
  void DrawHitYZPos(double *pos, int pid);
  void DrawTrack(double *pos1, double *pos2, int pid);
  int  GetPidColor(int pid);
  void get_command(void) const;

private:
  static TApplication *theApp;
};


#endif
