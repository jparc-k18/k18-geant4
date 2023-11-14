/*
  tof_check.cc
  Geometrical check for TOF detector 
  by using BC5 information.
  
  Toshi Gogami , 19Nov2014
*/

// ~~~ BC5 information from DCgeom.RealSts ~~~
const int n1=6;
double x1[n1] = {1260.56, 1262.61, 1264.66, 1266.71, 1268.77, 1270.82};
double val1[n1] = {3463.35, 3468.99, 3474.63, 3480.27, 3485.90, 3491.54};
// ~~~ TOF information from DCgeom.RealSts ~~~
const int n2=1;
double x2[n2] = {1351.19};
double val2[n2] = {3712.37};

void tof_check(){
  // ====== General conditions =======
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  
  // ===== Create Graph ===========
  TGraph* gr1 = new TGraph(n1,x1,val1);  
  gr1->SetName("gr1");
  TGraph* gr2 = new TGraph(n2,x2,val2);  
  gr2->SetName("gr2");
  
  // === Create empty histogram ==============
  TH2F* h_emp = new TH2F("h_emp","",100,1260.0,1360.0,100,3400.0,3750.0);
  
  // ===== Draw Histogram =====
  TCanvas* c1 = new TCanvas("c1","c1");
  h_emp->Draw();
  gr1->Draw("*lsame");
  gr2->Draw("*lsame");
  
  // ===== Options =======================
  gr1->SetMarkerStyle(25);
  gr2->SetMarkerStyle(24);
  cout << " " << gr1->Eval(x2[0]) << " mm" << endl;
  cout << " " << gr2->Eval(x2[0]) << " mm" << endl;
  
  // ==== if... =======
  double newpos = 1320.0;
  cout << "Now position=(" << newpos << "," << gr1->Eval(newpos) << ") mm" << endl;; 
  
  
}
