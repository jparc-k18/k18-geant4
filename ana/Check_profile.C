void Check_profile(){

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1111);
  //  TFile *fin = new TFile("Gene5_test.root");
    TFile *fin = new TFile("Gene7_test.root");
  TTree * t= new TTree;
  t=(TTree*)fin->Get("g4s2s");
  
  TCut cut0 = "SDC1[0].fE>0.&&SDC2[0].fE>0.";

  TCanvas *c0 = new TCanvas("c0","c0");
  c0->Divide(3,2);
  c0->cd(1);
  TH1D *hmom0 = new TH1D("hmom0","hmom0",800,800,1600);
  TH1D *hmom1 = new TH1D("hmom1","hmom1",800,800,1600);
  hmom0->GetXaxis()->SetTitle("Momentum [MeV/c]");
  hmom0->SetTitle("Momentum");
  t->Draw("p0>>hmom0");
  t->Draw("p0>>hmom1",cut0);
  hmom0->SetLineColor(1);
  hmom1->SetLineColor(2);
  hmom0->Draw();
  hmom1->Draw("same");
  
  c0->cd(2);
  TH1D *hang0 = new TH1D("hang0","hang0",500,0,16);
  TH1D *hang1 = new TH1D("hang1","hang1",500,0,16);
  hang0->GetXaxis()->SetTitle("#theta [degree]");
  hang0->SetTitle("#theta");
  t->Draw("theta0*180./acos(-1)>>hang0");
  t->Draw("theta0*180./acos(-1)>>hang1",cut0);
  hang0->SetLineColor(1);
  hang1->SetLineColor(2);
  hang0->Draw();
  hang1->Draw("same");

  c0->cd(3);
  TH2D * hxy0 = new TH2D("hxy0","hxy0"
			 ,200,-50,50,200,-30,30);
  t->Draw("y0:x0>>hxy0","","");
  hxy0->SetTitle("x vs y at Target");
  hxy0->GetXaxis()->SetTitle("x [mm]");
  hxy0->GetYaxis()->SetTitle("y [mm]");

  
  
  c0->cd(4);
  TH2D * hVP2_0 = new TH2D("hVP2_0","hVP2_0"
			   ,200,-800,800,300,-300,300);
  t->Draw("VP[1].fVy:VP[1].fVx>>hVP2_0",cut0,"");
  hVP2_0->SetTitle("VP2 No cut");
  hVP2_0->GetXaxis()->SetTitle("x [mm]");
  hVP2_0->GetYaxis()->SetTitle("y [mm]");


  c0->cd(5);
  TH2D * hVP4_0 = new TH2D("hVP4_0","hVP4_0"
			   ,200,-800,800,300,-300,300);
  t->Draw("VP[3].fVy:VP[3].fVx>>hVP4_0",cut0);
  hVP4_0->SetTitle("VP4 No cut");
  hVP4_0->GetXaxis()->SetTitle("x [mm]");
  hVP4_0->GetYaxis()->SetTitle("y [mm]");

  TCanvas *c1 = new TCanvas("c1","c1");
  c1->Divide(2,1);
  c1->cd(1);
  TH2D *hVP2_1[5];
  for(int i=0; i<5; ++i){
    std::cout<<"i"<<i<<std::endl;
    hVP2_1[i] = new TH2D(Form("hVP2_1_%d",i),
			 Form("hVP2_1_%d",i),
			 200,-800,800,
			 300,-300,300);
    TCut cut1 = cut0;
    cut1 += "theta0*180./acos(-1)<1.";
    cut1 += Form("p0>%f",1000.+(double)(100*i)-10.);
    cut1 += Form("p0<%f",1000.+(double)(100*i)+10.);
    t->Draw(Form("VP[1].fVy:VP[1].fVx>>hVP2_1_%d",i),cut1,"");
    hVP2_1[i]->SetTitle("VP2 #theta < 1deg");
    hVP2_1[i]->SetLineColor(i+1);
    hVP2_1[i]->SetMarkerColor(i+1);
    hVP2_1[i]->GetXaxis()->SetTitle("x [mm]");
    hVP2_1[i]->GetYaxis()->SetTitle("y [mm]");
  }
  
  hVP2_1[0]->Draw("scat");
  hVP2_1[0]->Draw("scat");
  hVP2_1[1]->Draw("scat&&same");
  hVP2_1[2]->Draw("scat&&same");
  hVP2_1[3]->Draw("scat&&same");
  hVP2_1[4]->Draw("scat&&same");


  c1->cd(2);
  TH2D *hVP4_1[5];
  for(int i=0; i<5; ++i){
    hVP4_1[i] = new TH2D(Form("hVP4_1_%d",i),
			 Form("hVP4_1_%d",i),
			 200,-800,800,
			 300,-300,300);
    TCut cut1 = cut0;
    cut1 += "theta0*180./acos(-1)<1.";
    cut1 += Form("p0>%f",1000.+(double)(100*i)-10.);
    cut1 += Form("p0<%f",1000.+(double)(100*i)+10.);
    t->Draw(Form("VP[3].fVy:VP[3].fVx>>hVP4_1_%d",i),cut1,"");
    hVP4_1[i]->SetTitle("VP4 #theta < 1deg");
    hVP4_1[i]->SetMarkerColor(i+1);
    hVP4_1[i]->SetLineColor(i+1);
    hVP4_1[i]->GetXaxis()->SetTitle("x [mm]");
    hVP4_1[i]->GetYaxis()->SetTitle("y [mm]");
  }
  hVP4_1[0]->Draw("scat");
  hVP4_1[1]->Draw("scat&&same");
  hVP4_1[2]->Draw("scat&&same");
  hVP4_1[3]->Draw("scat&&same");
  hVP4_1[4]->Draw("scat&&same");


  TCanvas *c2 = new TCanvas("c2","c2");
  c2->Divide(2,1);
  c2->cd(1);
  TH2D *hVP2_2[5];
  for(int i=0; i<5; ++i){
    std::cout<<"i"<<i<<std::endl;
    hVP2_2[i] = new TH2D(Form("hVP2_2_%d",i),
			 Form("hVP2_2_%d",i),
			 200,-800,800,
			 300,-300,300);
    TCut cut1 = cut0;
    cut1 += Form("theta0*180./acos(-1)>%f",(double)(i*2));
    cut1 += Form("theta0*180./acos(-1)<%f",(double)(i*2+2));
    cut1 += Form("p0>%f",1350.-10.);
    cut1 += Form("p0<%f",1350.+10.);
    t->Draw(Form("VP[1].fVy:VP[1].fVx>>hVP2_2_%d",i),cut1,"");
    hVP2_2[i]->SetTitle("VP2 #theta  dep");
    hVP2_2[i]->SetLineColor(i+1);
    hVP2_2[i]->SetMarkerColor(i+1);
    hVP2_2[i]->GetXaxis()->SetTitle("x [mm]");
    hVP2_2[i]->GetYaxis()->SetTitle("y [mm]");
  }
  
  hVP2_2[0]->Draw("scat");
  hVP2_2[0]->Draw("scat");
  hVP2_2[1]->Draw("scat&&same");
  hVP2_2[2]->Draw("scat&&same");
  hVP2_2[3]->Draw("scat&&same");
  hVP2_2[4]->Draw("scat&&same");
  
  c2->cd(2);
  TH2D *hVP4_2[5];
  for(int i=0; i<5; ++i){
    std::cout<<"i"<<i<<std::endl;
    hVP4_2[i] = new TH2D(Form("hVP4_2_%d",i),
			 Form("hVP4_2_%d",i),
			 200,-800,800,
			 300,-300,300);
    TCut cut1 = cut0;
    cut1 += Form("theta0*180./acos(-1)>%f",(double)(i*2));
    cut1 += Form("theta0*180./acos(-1)<%f",(double)(i*2+2));
    cut1 += Form("p0>%f",1350.-10.);
    cut1 += Form("p0<%f",1350.+10.);
    t->Draw(Form("VP[3].fVy:VP[3].fVx>>hVP4_2_%d",i),cut1,"");
    hVP4_2[i]->SetTitle("VP4 #theta  dep");
    hVP4_2[i]->SetLineColor(i+1);
    hVP4_2[i]->SetMarkerColor(i+1);
    hVP4_2[i]->GetXaxis()->SetTitle("x [mm]");
    hVP4_2[i]->GetYaxis()->SetTitle("y [mm]");
  }
  
  hVP4_2[0]->Draw("scat");
  hVP4_2[0]->Draw("scat");
  hVP4_2[1]->Draw("scat&&same");
  hVP4_2[2]->Draw("scat&&same");
  hVP4_2[3]->Draw("scat&&same");
  hVP4_2[4]->Draw("scat&&same");

  

}
