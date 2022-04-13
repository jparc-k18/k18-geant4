/* 
   Unit is mm.
*/

#include "EvDisp.hh"

#include <iostream>
#include <iomanip>
#include <cstdio>
#include <sstream>
#include <string>

#include "RadDeg.hh"

const int MaxChar = 200;


TApplication *EvDisp::theApp =0;

EvDisp::EvDisp(void)
{

  for (int i=0; i<=FiberLayerNum; i++) {
    for (int j=0; j<FiberSegMax; j++) {
      for (int k=0; k<4; k++) {
	fl[i][j][k] = NULL;
      }
    }
  }

  for (int i=0; i<segNumCrystal; i++)
    for (int j=0; j<4; j++) {
      cl[i][j] = NULL;
      fcl[i][j] = NULL;
    }

  for (int i=0; i<segNumRc; i++)
    for (int j=0; j<4; j++) {
      rl[i][j] = NULL;
    }

  for (int i=0; i<MaxHitEv; i++)
    farc1[i] = NULL;

  for (int i=0; i<MaxHitEv; i++)
    carc1[i] = NULL;

  hitFiberNum1_ = 0;
  hitCaloriNum1_ = 0;

  for (int i=0; i<LayerNumUVFiber; i++) 
    for (int j=0; j<4; j++)
      uvfl[i][j] = NULL;

  hitUVFiberNum_ = 0;
  hitRFiberNum_ = 0;
  for (int i=0; i<MaxHitEv; i++) {
    for (int j=0; j<2; j++) {
      huvfl[i][j] = NULL;
    }
    hrfl[i] = NULL;
  }

  for (int i=0; i<MaxHitEv; i++) {
    hitposXZ[i] = NULL;
    hitposYZ[i] = NULL;
  }
  hitPosXZNum_ = 0; 
  hitPosYZNum_ = 0; 

  for (int i=0; i<MaxHitEv; i++) {
    trlineXZ[i] = NULL;
    trlineYZ[i] = NULL;
    trlineXY[i] = NULL;
  }
  trLineNum_ = 0;

}

EvDisp::~EvDisp(void)
{
  for (int i=0; i<=FiberLayerNum; i++) {
    for (int j=0; j<FiberSegMax; j++) {
      for (int k=0; k<4; k++) {
	if (fl[i][j][k])
	  delete fl[i][j][k];
      }
    }
  }

  for (int i=0; i<segNumCrystal; i++)
    for (int j=0; j<4; j++) {
      if (cl[i][j])
	delete cl[i][j];
      if (fcl[i][j])
	delete fcl[i][j];
    }

  for (int i=0; i<segNumRc; i++)
    for (int j=0; j<4; j++) {
      if (rl[i][j])
	delete rl[i][j];
    }

  if (tgtXY)
    delete tgtXY;

  delete hbase[0];
  delete hbase[1];
  delete hbase[2];
  delete tp_[0];
  delete tp_[1];
  delete tp_[2];
  delete tc_;

  for (int i=0; i<LayerNumUVFiber; i++) 
    for (int j=0; j<4; j++)
      if (uvfl[i][j])
	delete uvfl[i][j];

  for (int i=0; i<MaxHitEv; i++) {
    for (int j=0; j<2; j++) {
      if (huvfl[i][j])
	delete huvfl[i][j];
    }
    if (hrfl[i])
      delete hrfl[i];

  }

  for (int i=0; i<LayerNumUVFiber; i++) 
    delete tp2_[i];
  delete tc2_;
}

void EvDisp::Initialize(void)
{
  if( !theApp ){
    theApp=new TApplication( "App", 0, 0 );
  }
  
  tc_ = new TCanvas("canvas","Event Display",800,400);
  tp_[0] = new TPad("pad0","Event",0.0, 0.0, 0.5, 1.0, 10);
  tp_[0]->Draw();

  tp_[1] = new TPad("pad1","Event",0.5, 0.0, 1.0, 0.5, 10);
  tp_[1]->Draw();

  tp_[2] = new TPad("pad2","Event",0.5, 0.5, 1.0, 1.0, 10);
  tp_[2]->Draw();

  tp_[0]->cd();
  std::cout << "tp_->Draw();" << std::endl;

  hbase[0] = new TH2F("hbase1","",600, -180, 180, 600, -180,180);
  hbase[0]->GetXaxis()->SetTitle("X (mm)");
  hbase[0]->GetYaxis()->SetTitle("Y (mm)");
  hbase[0]->Draw();

  tgtXY = new TArc(0., 0., CyLH2TgtR/2.);
  tgtXY->Draw("same");

  char buf[100];
  for (int l=1; l<=FiberLayerNum; l += 2) {
    double theta = FiberAngle[l];
    int SegNum = segNumFiber[l];
    for (int i=0; i<SegNum; i++) {
      double x0, y0;
      x0 = rFiber[l]*cos((double)i*theta*Deg2Rad);
      y0 = rFiber[l]*sin((double)i*theta*Deg2Rad);
      DrawFiber(l, i, x0, y0, theta);
    }
  }

  for (int i=0; i<segNumCrystal; i++) {
    double x0, y0;
    x0 = (rCrystal+CrystalSize_Y/2.)*cos((double)i*CrystalAngle*Deg2Rad);
    y0 = (rCrystal+CrystalSize_Y/2.)*sin((double)i*CrystalAngle*Deg2Rad);
    DrawCalori(i, x0, y0, CrystalAngle);
  }
  
  for (int i=0; i<segNumForwardCrystal; i++) {
    if ((i>=2 && i<=8) || i>=12&&i<=18) {
      double x0, y0;
      x0 = (rForwardCrystal+ForwardCrystalSize_Y/2.)*cos((double)i*CrystalAngle*Deg2Rad);
      y0 = (rForwardCrystal+ForwardCrystalSize_Y/2.)*sin((double)i*CrystalAngle*Deg2Rad);
      double Offset;
      Offset = 30.;
      DrawCaloriForward(i, x0, y0, CrystalAngle);
    }
  }

  for (int i=0; i<segNumRc; i++) {
    double x0, y0;
    x0 = (rRc+RcSize_Y/2.)*cos((double)i*RcAngle*Deg2Rad);
    y0 = (rRc+RcSize_Y/2.)*sin((double)i*RcAngle*Deg2Rad);
    DrawRc(i, x0, y0, RcAngle);
  }

  tc_->cd();
  tp_[2]->cd();
  hbase[2] = new TH2F("hbase3","", 100, -200, 400, 100, -160, 160);
  hbase[2]->GetXaxis()->SetTitle("Z (mm)");
  hbase[2]->GetYaxis()->SetTitle("Y (mm)");
  hbase[2]->Draw();

  InitializeYZPlane();

  tc_->cd();
  tp_[1]->cd();
  hbase[1] = new TH2F("hbase2","",100, -200, 400, 100, -160, 160);
  hbase[1]->GetXaxis()->SetTitle("Z (mm)");
  hbase[1]->GetYaxis()->SetTitle("X (mm)");
  hbase[1]->Draw();

  InitializeXZPlane();

  tc_->cd();
  tc_->Update();

  std::cout << "tc_->Update();" << std::endl;

  tc2_ = new TCanvas("canvas2","Event Display",300,700);
  for (int i=0; i<LayerNumUVFiber; i++) {
    sprintf(buf, "tp2_%d", i);
    tp2_[i] = new TPad(buf,"Event",0.0, i*0.2, 1.0, (i+1)*0.2, 10);
    tp2_[i]->Draw();
  }

  DrawUVFiber();


}

void EvDisp::InitializeXZPlane(void)
{

  double p[4][2];
  
  p[0][0] = -150;   p[0][1] = 20; 
  p[1][0] = -150;   p[1][1] = -20; 
  p[2][0] =  150;   p[2][1] = -20; 
  p[3][0] =  150;   p[3][1] = 20; 
  
  DrawBox(tgtlXZ, p, kBlack);

  for (int i=1; i<=FiberLayerNum; i++) {
    p[0][0] = -150.;               p[0][1] = rFiber[i];
    p[1][0] = -150.+FiberSizeZ[i]; p[1][1] = rFiber[i];

    flXZ[i][0] = new TLine(p[0][0], p[0][1], p[1][0], p[1][1]);
    flXZ[i][0]->Draw("same");

    p[0][0] = -150.;               p[0][1] = -rFiber[i];
    p[1][0] = -150.+FiberSizeZ[i]; p[1][1] = -rFiber[i];

    flXZ[i][1] = new TLine(p[0][0], p[0][1], p[1][0], p[1][1]);
    flXZ[i][1]->Draw("same");
  }

  p[0][0] = -150.;               p[0][1] = rCrystal; 
  p[1][0] = -150.+CrystalSize_Z; p[1][1] = rCrystal; 
  p[2][0] = -150.+CrystalSize_Z; p[2][1] = rCrystal+CrystalSize_Y; 
  p[3][0] = -150.;               p[3][1] = rCrystal+CrystalSize_Y; 

  DrawBox(clXZ1, p, 9);

  p[0][0] = -150.;               p[0][1] = -rCrystal; 
  p[1][0] = -150.+CrystalSize_Z; p[1][1] = -rCrystal; 
  p[2][0] = -150.+CrystalSize_Z; p[2][1] = -rCrystal-CrystalSize_Y; 
  p[3][0] = -150.;               p[3][1] = -rCrystal-CrystalSize_Y; 

  DrawBox(clXZ2, p, 9);

  p[0][0] = -150.;          p[0][1] = rRc; 
  p[1][0] = -150.+RcSize_Z; p[1][1] = rRc; 
  p[2][0] = -150.+RcSize_Z; p[2][1] = rRc+RcSize_Y; 
  p[3][0] = -150.;          p[3][1] = rRc+RcSize_Y; 
  DrawBox(rclXZ1, p, 7);

  p[0][0] = -150.;          p[0][1] = -rRc; 
  p[1][0] = -150.+RcSize_Z; p[1][1] = -rRc; 
  p[2][0] = -150.+RcSize_Z; p[2][1] = -rRc-RcSize_Y; 
  p[3][0] = -150.;          p[3][1] = -rRc-RcSize_Y; 
  DrawBox(rclXZ1, p, 7);

}



void EvDisp::InitializeYZPlane(void)
{
  double p[4][2];

  p[0][0] = -150;   p[0][1] = 20; 
  p[1][0] = -150;   p[1][1] = -20; 
  p[2][0] =  150;   p[2][1] = -20; 
  p[3][0] =  150;   p[3][1] = 20; 
  DrawBox(tgtlYZ, p, kBlack);

  for (int i=1; i<=FiberLayerNum; i++) {
    p[0][0] = -150.;               p[0][1] = rFiber[i];
    p[1][0] = -150.+FiberSizeZ[i]; p[1][1] = rFiber[i];

    flYZ[i][0] = new TLine(p[0][0], p[0][1], p[1][0], p[1][1]);
    flYZ[i][0]->Draw("same");

    p[0][0] = -150.;               p[0][1] = -rFiber[i];
    p[1][0] = -150.+FiberSizeZ[i]; p[1][1] = -rFiber[i];

    flYZ[i][1] = new TLine(p[0][0], p[0][1], p[1][0], p[1][1]);
    flYZ[i][1]->Draw("same");
  }

  p[0][0] = -150.;               p[0][1] = rCrystal; 
  p[1][0] = -150.+CrystalSize_Z; p[1][1] = rCrystal; 
  p[2][0] = -150.+CrystalSize_Z; p[2][1] = rCrystal+CrystalSize_Y; 
  p[3][0] = -150.;               p[3][1] = rCrystal+CrystalSize_Y; 
  DrawBox(clYZ1, p, 9);

  p[0][0] = -150.;               p[0][1] = -rCrystal; 
  p[1][0] = -150.+CrystalSize_Z; p[1][1] = -rCrystal; 
  p[2][0] = -150.+CrystalSize_Z; p[2][1] = -rCrystal-CrystalSize_Y; 
  p[3][0] = -150.;               p[3][1] = -rCrystal-CrystalSize_Y; 
  DrawBox(clYZ2, p, 9);

  p[0][0] = -150.+CrystalSize_Z;                      p[0][1] = rForwardCrystal; 
  p[1][0] = -150.+CrystalSize_Z+ForwardCrystalSize_Z; p[1][1] = rForwardCrystal; 
  p[2][0] = -150.+CrystalSize_Z+ForwardCrystalSize_Z; p[2][1] = rCrystal+CrystalSize_Y; 
  p[3][0] = -150.+CrystalSize_Z;                      p[3][1] = rCrystal+CrystalSize_Y; 

  DrawBox(fclYZ1, p, 6);

  p[0][0] = -150.+CrystalSize_Z;                      p[0][1] = -(rForwardCrystal); 
  p[1][0] = -150.+CrystalSize_Z+ForwardCrystalSize_Z; p[1][1] = -(rForwardCrystal); 
  p[2][0] = -150.+CrystalSize_Z+ForwardCrystalSize_Z; p[2][1] = -(rCrystal+CrystalSize_Y); 
  p[3][0] = -150.+CrystalSize_Z;                      p[3][1] = -(rCrystal+CrystalSize_Y); 

  DrawBox(fclYZ2, p, 6);

  p[0][0] = -150.;          p[0][1] = rRc; 
  p[1][0] = -150.+RcSize_Z; p[1][1] = rRc; 
  p[2][0] = -150.+RcSize_Z; p[2][1] = rRc+RcSize_Y; 
  p[3][0] = -150.;          p[3][1] = rRc+RcSize_Y; 
  DrawBox(rclYZ1, p, 7);

  p[0][0] = -150.;          p[0][1] = -rRc; 
  p[1][0] = -150.+RcSize_Z; p[1][1] = -rRc; 
  p[2][0] = -150.+RcSize_Z; p[2][1] = -rRc-RcSize_Y; 
  p[3][0] = -150.;          p[3][1] = -rRc-RcSize_Y; 
  DrawBox(rclYZ1, p, 7);

  p[0][0] = -150.+RcSize_Z-ForwardRcSize_Z; p[0][1] = rForwardRc; 
  p[1][0] = -150.+RcSize_Z;                 p[1][1] = rForwardRc; 
  p[2][0] = -150.+RcSize_Z;                 p[2][1] = rRc; 
  p[3][0] = -150.+RcSize_Z-ForwardRcSize_Z; p[3][1] = rRc;

  DrawBox(frclYZ1, p, 7);

  p[0][0] = -150.+RcSize_Z-ForwardRcSize_Z; p[0][1] = -(rForwardRc); 
  p[1][0] = -150.+RcSize_Z;                 p[1][1] = -(rForwardRc); 
  p[2][0] = -150.+RcSize_Z;                 p[2][1] = -rRc; 
  p[3][0] = -150.+RcSize_Z-ForwardRcSize_Z; p[3][1] = -rRc;

  DrawBox(frclYZ1, p, 7);


}

void EvDisp::EndOfEvent()
{
  tc_->Update();
  tc2_->Update();
  get_command();

  for (int i=0; i<hitFiberNum1_; i++)
    if (farc1[i])
      delete farc1[i];
  hitFiberNum1_ = 0;

  for (int i=0; i<hitCaloriNum1_; i++)
    if (carc1[i])
      delete carc1[i];
  hitCaloriNum1_ = 0;

  for (int i=0; i<hitUVFiberNum_; i++)
    for (int j=0; j<2; j++) 
      if (huvfl[i][j])
	delete huvfl[i][j];
  hitUVFiberNum_ = 0;

  for (int i=0; i<hitRFiberNum_; i++)
    if (hrfl[i])
      delete hrfl[i];
  hitRFiberNum_ = 0;


  for (int i=0; i<hitPosXZNum_; i++)
    if (hitposXZ[i])
      delete hitposXZ[i];
  hitPosXZNum_ = 0;

  for (int i=0; i<hitPosYZNum_; i++)
    if (hitposYZ[i])
      delete hitposYZ[i];
  hitPosYZNum_ = 0;

  for (int i=0; i<trLineNum_; i++) {
    if (trlineXZ[i])
      delete trlineXZ[i];
    if (trlineYZ[i])
      delete trlineYZ[i];
    if (trlineXY[i])
      delete trlineXY[i];
  }
  trLineNum_ = 0;
}


void EvDisp::DrawFiber(int layer, int seg, double x0, double y0, double theta)
{
  double p[4][2];
  double x=FiberSizeX[layer];
  double y=FiberSizeY[layer];

  calcRotate((double)(seg)*theta, y/2., x/2., p[0]);
  calcRotate((double)(seg)*theta, -y/2., x/2., p[1]);
  calcRotate((double)(seg)*theta, -y/2., -x/2., p[2]);
  calcRotate((double)(seg)*theta, y/2., -x/2., p[3]);

  for (int i=0; i<4; i++) {
    p[i][0] += x0;
    p[i][1] += y0;
  }

  for (int i=0; i<4; i++) {
    if (i>=0&&i<=2)
      fl[layer][seg][i] = new TLine(p[i][0], p[i][1], p[i+1][0], p[i+1][1]);
    else if (i==3)
      fl[layer][seg][i] = new TLine(p[3][0], p[3][1], p[0][0], p[0][1]);
    fl[layer][seg][i]->Draw();
  }
}

void EvDisp::DrawCalori(int seg, double x0, double y0, double thetaCal)
{
  double p[4][2];
  double x1=CrystalSize_X;
  double x2=CrystalSize_X2;

  calcRotate((double)(seg)*thetaCal, CrystalSize_Y/2., x2/2., p[0]);
  calcRotate((double)(seg)*thetaCal, -CrystalSize_Y/2., x1/2., p[1]);
  calcRotate((double)(seg)*thetaCal, -CrystalSize_Y/2., -x1/2., p[2]);
  calcRotate((double)(seg)*thetaCal, CrystalSize_Y/2., -x2/2., p[3]);

  for (int i=0; i<4; i++) {
    p[i][0] += x0;
    p[i][1] += y0;
  }

  for (int i=0; i<4; i++) {
    if (i>=0&&i<=2)
      cl[seg][i] = new TLine(p[i][0], p[i][1], p[i+1][0], p[i+1][1]);
    else if (i==3)
      cl[seg][i] = new TLine(p[3][0], p[3][1], p[0][0], p[0][1]);
    cl[seg][i]->SetLineColor(9);
    cl[seg][i]->SetLineWidth(2);
    cl[seg][i]->Draw();
  }
}


void EvDisp::DrawCaloriForward(int seg, double x0, double y0, double thetaCal)
{
  double p[4][2];
  double x1= ForwardCrystalSize_X;
  double x2= ForwardCrystalSize_X2;

  calcRotate((double)(seg)*thetaCal, ForwardCrystalSize_Y/2., x2/2., p[0]);
  calcRotate((double)(seg)*thetaCal, -ForwardCrystalSize_Y/2., x1/2., p[1]);
  calcRotate((double)(seg)*thetaCal, -ForwardCrystalSize_Y/2., -x1/2., p[2]);
  calcRotate((double)(seg)*thetaCal, ForwardCrystalSize_Y/2., -x2/2., p[3]);

  for (int i=0; i<4; i++) {
    p[i][0] += x0;
    p[i][1] += y0;
  }

  for (int i=0; i<4; i++) { 
   if (i>=0&&i<=2)
     fcl[seg][i] = new TLine(p[i][0], p[i][1], p[i+1][0], p[i+1][1]);
   else if (i==3)
     fcl[seg][i] = new TLine(p[3][0], p[3][1], p[0][0], p[0][1]);
   fcl[seg][i]->SetLineColor(6);
   fcl[seg][i]->SetLineWidth(2);
   fcl[seg][i]->Draw();
  }
}

void EvDisp::DrawRc(int seg, double x0, double y0, double thetaCal)
{
  double p[4][2];
  double x1=RcSize_X;
  double x2=RcSize_X2;

  calcRotate((double)(seg)*thetaCal, RcSize_Y/2., x2/2., p[0]);
  calcRotate((double)(seg)*thetaCal, -RcSize_Y/2., x1/2., p[1]);
  calcRotate((double)(seg)*thetaCal, -RcSize_Y/2., -x1/2., p[2]);
  calcRotate((double)(seg)*thetaCal, RcSize_Y/2., -x2/2., p[3]);

  for (int i=0; i<4; i++) {
    p[i][0] += x0;
    p[i][1] += y0;
  }

  for (int i=0; i<4; i++) {
    if (i>=0&&i<=2)
      rl[seg][i] = new TLine(p[i][0], p[i][1], p[i+1][0], p[i+1][1]);
    else if (i==3)
      rl[seg][i] = new TLine(p[3][0], p[3][1], p[0][0], p[0][1]);

    rl[seg][i]->SetLineColor(7);
    rl[seg][i]->SetLineWidth(2);
    rl[seg][i]->Draw();
  }
}

void EvDisp::DrawUVFiber(void)
{
  char buf[100];

  tc2_->cd();
  for (int i=0; i<LayerNumUVFiber; i++) {
    tp2_[i]->cd();
    sprintf(buf, "hbase2_%d", i);
    hbase2[i] = new TH2F(buf, "", 100, -50, 600, 100, -50, 450);
    hbase2[i]->GetXaxis()->SetTitle("R (mm)");
    hbase2[i]->GetYaxis()->SetTitle("Z (mm)");
    hbase2[i]->Draw();

    double pos[4][2];
    
    double R = 2.*Pi*rFiber[(i+1)*2];
    double Z = FiberSizeZ[(i+1)*2];
    pos[0][0] = 0.;      pos[0][1] = 0.;
    pos[1][0] = R;       pos[1][1] = 0.;
    pos[2][0] = R;       pos[2][1] = Z;
    pos[3][0] = 0.;      pos[3][1] = Z;

    for (int j=0; j<4; j++) {
      if (j==3)
	uvfl[i][j] = new TLine(pos[3][0], pos[3][1], pos[0][0], pos[0][1]);
      else
	uvfl[i][j] = new TLine(pos[j][0], pos[j][1], pos[j+1][0], pos[j+1][1]);
      uvfl[i][j]->Draw("same");
    }
  }
  tc2_->Update();

}

void EvDisp::DrawBox(TLine *l[4], double p[4][2], int color)
{
  for (int i=0; i<4; i++) {
    if (i>=0&&i<=2)
      l[i] = new TLine(p[i][0], p[i][1], p[i+1][0], p[i+1][1]);
    else if (i==3)
      l[i] = new TLine(p[3][0], p[3][1], p[0][0], p[0][1]);
    l[i]->SetLineColor(color);
    l[i]->Draw();
  }
}

void EvDisp::calcRotate(double theta, double x, double y, double *rotpos)
{
  rotpos[0] = cos(theta*Deg2Rad)*x - sin(theta*Deg2Rad)*y; 
  rotpos[1] = sin(theta*Deg2Rad)*x + cos(theta*Deg2Rad)*y; 
}

void EvDisp::DrawHitFiber(int layer, int seg, int pid, double de)
{
  if (!(layer == 1 || layer == 3 || layer == 5 || layer == 7 || layer == 9))
    return;
  
  double x0 = rFiber[layer]*cos((double)seg*FiberAngle[layer]*Deg2Rad);
  double y0 = rFiber[layer]*sin((double)seg*FiberAngle[layer]*Deg2Rad);
  farc1[hitFiberNum1_] = new TArc(x0, y0, 1.*de/0.1);
  if (1.*de/0.1>10.) {
    farc1[hitFiberNum1_]->SetLineColor(GetPidColor(pid));
    farc1[hitFiberNum1_]->SetFillStyle(0);
  } else 
    farc1[hitFiberNum1_]->SetFillColor(GetPidColor(pid));

  tc_->cd();
  tp_[0]->cd();
  farc1[hitFiberNum1_]->Draw("same");
  hitFiberNum1_++;
}

void EvDisp::DrawHitCalori(int layer, int seg, int pid, double de)
{

  double x0 = rCaloriCenter[layer]*cos((double)seg*CaloriAngle[layer]*Deg2Rad);
  double y0 = rCaloriCenter[layer]*sin((double)seg*CaloriAngle[layer]*Deg2Rad);

  carc1[hitCaloriNum1_] = new TArc(x0, y0, 1.*de/10.);
  if (1.*de/10.>10.) {
    carc1[hitCaloriNum1_]->SetLineColor(GetPidColor(pid));
    carc1[hitCaloriNum1_]->SetFillStyle(0);
  } else 
    carc1[hitCaloriNum1_]->SetFillColor(GetPidColor(pid));

  tc_->cd();
  tp_[0]->cd();
  carc1[hitCaloriNum1_]->Draw("same");
  hitCaloriNum1_++;
}

void EvDisp::DrawHitUVFiber(int layer, int seg, int pid, double de)
{
  if (!(layer == 2 || layer == 4 || layer == 6 || layer == 8 || layer == 10))
    return;

  tc2_->cd();
  tp2_[layer/2-1]->cd();

  double theta = tiltAngleFiber[layer];
  double r0 = (double)seg*FiberSizeX[layer]*cos(theta*Deg2Rad);
  double z0 = (double)seg*FiberSizeX[layer]*sin(theta*Deg2Rad);
  double R = 2.*Pi*rFiber[layer];
  double Z = FiberSizeZ[layer];
  double p[4][2];

  if (tiltAngleFiber[layer]>=0) {
    p[0][0] = z0*tan(theta*Deg2Rad)+r0; p[0][1] = 0.;
    p[1][0] = 0.;                       p[1][1] = r0/tan(theta*Deg2Rad)+z0;
    p[2][0] = R;                        p[2][1] = r0/tan(theta*Deg2Rad)+z0;
    p[3][0] = z0*tan(theta*Deg2Rad)+r0; p[3][1] = Z;
  } else {
    p[0][0] = z0*tan(theta*Deg2Rad)+r0; p[0][1] = 0.;
    p[1][0] = R;                        p[1][1] = -(R-r0)/tan(theta*Deg2Rad)+z0;
    p[2][0] = 0.;                       p[2][1] = -(R-r0)/tan(theta*Deg2Rad)+z0;
    p[3][0] = z0*tan(theta*Deg2Rad)+r0; p[3][1] = Z;
  }

  huvfl[hitUVFiberNum_][0] = new TLine(p[0][0], p[0][1], p[1][0], p[1][1]) ;
  huvfl[hitUVFiberNum_][1] = new TLine(p[2][0], p[2][1], p[3][0], p[3][1]) ;

  for (int i=0; i<2; i++) {
    huvfl[hitUVFiberNum_][i]->SetLineColor(GetPidColor(pid));
    huvfl[hitUVFiberNum_][i]->Draw("same");
  }
  hitUVFiberNum_++;
}

void EvDisp::DrawHitRFiber(int layer, int seg, int pid, double de)
{
  if (!(layer == 1 || layer == 3 || layer == 5 || layer == 7 || layer == 9))
    return;

  tc2_->cd();
  tp2_[(int)(layer/2)]->cd();

  double r0 = (double)seg*FiberSizeX[layer];
  double Z = FiberSizeZ[layer];
  double p[2][2];

  p[0][0] = r0; p[0][1] = 0.;
  p[1][0] = r0; p[1][1] = Z;

  hrfl[hitRFiberNum_] = new TLine(p[0][0], p[0][1], p[1][0], p[1][1]) ;
  hrfl[hitRFiberNum_]->SetLineColor(GetPidColor(pid));
  hrfl[hitRFiberNum_]->Draw("same");

  hitRFiberNum_++;
}

void EvDisp::DrawHitXZPos(double *pos, int pid)
{
  hitposXZ[hitPosXZNum_] = new TArc(pos[2], pos[0], 3.);
  hitposXZ[hitPosXZNum_]->SetFillColor(GetPidColor(pid));

  tc_->cd();
  tp_[1]->cd();
  hitposXZ[hitPosXZNum_]->Draw("same");
  hitPosXZNum_++;

}

void EvDisp::DrawHitYZPos(double *pos, int pid)
{
  hitposYZ[hitPosYZNum_] = new TArc(pos[2], pos[1], 3.);
  hitposYZ[hitPosYZNum_]->SetFillColor(GetPidColor(pid));

  tc_->cd();
  tp_[2]->cd();
  hitposYZ[hitPosYZNum_]->Draw("same");
  hitPosYZNum_++;

}

void EvDisp::DrawTrack(double *pos1, double *pos2, int pid)
{
  tc_->cd();

  trlineXY[trLineNum_] = new TLine(pos1[0], pos1[1], pos2[0], pos2[1]);
  trlineXY[trLineNum_]->SetLineColor(GetPidColor(pid));
  tp_[0]->cd();
  trlineXY[trLineNum_]->Draw("same");

  trlineXZ[trLineNum_] = new TLine(pos1[2], pos1[0], pos2[2], pos2[0]);
  trlineXZ[trLineNum_]->SetLineColor(GetPidColor(pid));
  tp_[1]->cd();
  trlineXZ[trLineNum_]->Draw("same");

  trlineYZ[trLineNum_] = new TLine(pos1[2], pos1[1], pos2[2], pos2[1]);
  trlineYZ[trLineNum_]->SetLineColor(GetPidColor(pid));
  tp_[2]->cd();
  trlineYZ[trLineNum_]->Draw("same");

  trLineNum_++;

}

int EvDisp::GetPidColor(int pid)
{
  if (pid == 1) // pi+
    return kBlue;
  else if (pid == 2) // pi-
    return kGreen;
  else if (pid == 4 || pid == 5 || pid == 9) //  e-, e+, gamma
    return kCyan;
  else if (pid == 10) //  proton
    return kRed;
  else if (pid == 14) //  proton_d
    return kMagenta;
  else {
    if (pid ==1)
      std::cout << "pi+" << std::endl;
    else if (pid ==2)
      std::cout << "pi-" << std::endl;
    else if (pid ==3)
      std::cout << "kaon+" << std::endl;
    else if (pid ==4)
      std::cout << "kaon-" << std::endl;
    else if (pid ==5)
      std::cout << "e-" << std::endl;
    else if (pid ==6)
      std::cout << "e+" << std::endl;
    else if (pid ==7)
      std::cout << "mu-" << std::endl;
    else if (pid ==8)
      std::cout << "mu+" << std::endl;
    else if (pid ==9)
      std::cout << "gamma" << std::endl;
    else if (pid ==10)
      std::cout << "proton" << std::endl;
    else if (pid ==11)
      std::cout << "deuteron" << std::endl;
    else if (pid ==12)
      std::cout << "sigma+" << std::endl;
    else if (pid ==13)
      std::cout << "sigma-" << std::endl;
    else if (pid ==14)
      std::cout << "proton_d" << std::endl;
    else 
      std::cout << pid << std::endl;
    return kBlack;
  }
}

void EvDisp::get_command(void) const
{
  char ch;
  char data[100];
  static int stat=0;
  static int Nevent=0;
  static int ev=0;
  
  if (stat == 1 && Nevent > 0 && ev<Nevent) {
    ev++;
    return;
  } 
  if (ev==Nevent) {
    stat=0;
    ev=0;
  }

  if (stat == 0) {
    printf("q|n|p>");

    /* get command */
    scanf("%c",&ch);
    if (ch!='\n')
      while(getchar() != '\n');

    switch (ch) {
    case 'q': exit(0);
    case 'n':
      stat = 1;
      do {
	printf("event#>");
	scanf("%s",data);
      } while ((Nevent=atoi(data))<=0);
      std::cout << "Continue " << Nevent << "event" << std::endl;
      break;
    case 'p':
      theApp->Run(kTRUE);
      break;
    }
  }
}
