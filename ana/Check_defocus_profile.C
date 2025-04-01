void Check_defocus_profile(){

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1111);
  TFile *fin = new TFile("Gene6_Q1_100_Q2_100.root");
  TTree * t= new TTree;
  t=(TTree*)fin->Get("g4s2s");
  
  TCut cut0 = "TOF[0].fE>0.&&SDC1[0].fE>0.&&SDC2[0].fE>0.&&SDC3[0].fE>0.&&SDC4[0].fE>0.&&SDC5[0].fE>0.&&WC[0].fE>0.";

  TCanvas *c0 = new TCanvas("c0","c0");
  c0->Divide(2,1);
  c0->cd(1)->SetGrid(1);
  TH2D * hVP2_0 = new TH2D("hVP2_0","hVP2_0"
			   ,200,-800,800,300,-300,300);
  t->Draw("VP[1].fVy:VP[1].fVx>>hVP2_0",cut0,"",1);
  hVP2_0->SetTitle("VP2 Defocus run");
  hVP2_0->GetXaxis()->SetTitle("x [mm]");
  hVP2_0->GetYaxis()->SetTitle("y [mm]");
  hVP2_0->GetXaxis()->SetRangeUser(-400,400);
  hVP2_0->GetYaxis()->SetRangeUser(-100,100);
 

  c0->cd(2)->SetGrid(1);
  TH2D * hVP4_0 = new TH2D("hVP4_0","hVP4_0"
			   ,200,-800,800,300,-300,300);
  t->Draw("VP[3].fVy:VP[3].fVx>>hVP4_0",cut0,"",1);
  hVP4_0->SetTitle("VP4 Defocus run");
  hVP4_0->GetXaxis()->SetTitle("x [mm]");
  hVP4_0->GetYaxis()->SetTitle("y [mm]");
  hVP4_0->GetXaxis()->SetRangeUser(-400,400);
  hVP4_0->GetYaxis()->SetRangeUser(-100,100);

  

}
