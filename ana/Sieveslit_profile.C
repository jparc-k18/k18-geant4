double curcle (double x,double y){
  return x*x + y*y;
}

void Sieveslit_profile(){
 
  gStyle->SetOptStat(0);
  // gStyle->SetOptFit(1111);
  gStyle->SetStatX(0.9);
  gStyle->SetStatY(0.9);

  
  //  TFile *f = new TFile("Gene6_Q1_100_Q2_100_v2.root");
  TFile *f = new TFile("Gene8_test_v2.root");
  TTree *t = new TTree;
  t=(TTree*)f->Get("g4s2s");

  TTreeReader reader("g4s2s",f);
  TTreeReaderValue<std::vector<TParticle>> SDC1(reader,"SDC1");
  TTreeReaderValue<std::vector<TParticle>> SDC2(reader,"SDC2");
  TTreeReaderValue<std::vector<TParticle>> VP(reader,"VP");
  TTreeReaderValue<Double_t> x0(reader,"x0");
  TTreeReaderValue<Double_t> y0(reader,"y0");
  TTreeReaderValue<Double_t> xp0(reader,"xp0");
  TTreeReaderValue<Double_t> yp0(reader,"yp0");
  
  TCut cut0 = "SDC1[0].fE>0.&&SDC2[0].fE>0.";
  TCut cut1 = "";

  TH2D *hist = new TH2D("hist","hist",1600,-800,800,300,-300,300);
  TH2D *hist2 = new TH2D("his2","hist2",1600,-800,800,300,-300,300);
  TH2D *VP2 = new TH2D("VP2","VP2",1600,-800,800,600,-300,300);
  TH2D *VP4 = new TH2D("VP4","VP4",1600,-800,800,600,-300,300);
  TH1D *prjx2 = new TH1D("prjx2","prjx2",600,-300,300);
  TH1D *prjy2 = new TH1D("prjy2","prjy2",300,-150,150);
  TH1D *prjx4 = new TH1D("prjx4","prjx4",600,-300,300);
  TH1D *prjy4 = new TH1D("prjy4","prjy4",300,-150,150);

  // TF1 *gaus = new TF1("gaus","[0]*TMath::Gaus(-(x-[1])*(x-[1])/(2*[2]*[2]))",-300.,300.);
  
  double z0 = -4926.9;
  double z2 = -3246.2;
  double z4 = -2303.57;
  double sigma = 0.002;
  TRandom rand;
    
  while (reader.Next()){
    int ievent = reader.GetCurrentEntry();
    int SDC1_multi = 0;
    int SDC2_multi = 0;
    for (const auto &p:(*SDC1)){
      if (true && p.GetPDG()->Charge() != 0){
	SDC1_multi ++ ;
      }
    }
    for (const auto &p:(*SDC2)){
      if (true && p.GetPDG()->Charge() != 0){
	SDC2_multi ++ ;
      }
    }

    if(true){
      reader.SetEntry(ievent);
      double x2 = (*x0)+(*xp0)*(z2-z0);
      double y2 = (*y0)+(*yp0)*(z2-z0);
      double x4 = (*x0)+(*xp0)*(z4-z0);
      double y4 = (*y0)+(*yp0)*(z4-z0);
      double xr2 = (*x0)+(rand.Gaus(*xp0,sigma))*(z2-z0);
      double yr2 = (*y0)+(rand.Gaus(*yp0,sigma))*(z2-z0);
      double xr4 = (*x0)+(rand.Gaus(*xp0,sigma))*(z4-z0);
      double yr4 = (*y0)+(rand.Gaus(*yp0,sigma))*(z4-z0);

      //VP2->Fill(x2,y2);
      //VP4->Fill(x4,y4);
      
      for (const auto &p : (*VP)){
	if (true && p.GetPDG()->Charge() != 0){
	  int layer = p.GetMother(1);

	  if (layer == 2){
	    VP2->Fill(xr2,yr2);
	    for (int i=0; i<24; i++){
	      for (int j=0; j<11; j++){
		if (curcle(x2-(-230.+20.*i), y2-(-100.+20.*j))<=1.*1.){
		  hist->Fill(xr2, yr2);
		  if (y2 >= -10. && y2 <= 10.){
		    prjx2->Fill(xr2);
		  }
		  if (x2 >= 0. && x2 <= 20.){
		    prjy2->Fill(yr2);
		  }
		}		
	      }
	    }
	  }
	  if (layer == 4){
	    VP4->Fill(xr4,yr4);
	    for (int i=0; i<39; i++){
	      for (int j=0; j<11; j++){
		if (curcle(x4-(-380.+20.*i),y4-(-100.+20.*j))<=1.*1.){
		  if (y4 >= -10. && y4 <= 10.){
		    prjx4->Fill(xr4);
		  }
		  if (x4 >= -10. && x4 <= 10.){
		    prjy4->Fill(yr4);
		  }
		  hist2->Fill(xr4,yr4);
		}
	      }
	    }
	  }
	}	  
      }	
    }
  }

  
  TCanvas *c0 = new TCanvas("c0","c0",1500,900);
  TCanvas *c1 = new TCanvas("c1","c1",1500,900);
  TCanvas *c2 = new TCanvas("c2","c2",1500,900);
  // TCanvas *c3 = new TCanvas("c3","c3",1600,1200);
  // TCanvas *c4 = new TCanvas("c4","c4");
  // TCanvas *c5 = new TCanvas("c5","c5");
  // TCanvas *c6 = new TCanvas("c6","c6");
  // TCanvas *c7 = new TCanvas("c7","c7");
  c0->Divide(2,2);
  c1->Divide(2,2);
  c2->Divide(2,2);
  // c3->Divide(2,2);
  
  c0->cd(3);
  hist->SetTitle("SS1");
  hist->GetXaxis()->SetRangeUser(-300,300);
  hist->GetYaxis()->SetRangeUser(-150,150);
  hist->GetXaxis()->SetTitle("x[mm]");
  hist->GetYaxis()->SetTitle("y[mm]");
  hist->Draw("colz");
  
  c0->cd(4);
  hist2->SetTitle("SS2");
  hist2->GetXaxis()->SetRangeUser(-300,300);
  hist2->GetYaxis()->SetRangeUser(-150,150);
  hist2->GetXaxis()->SetTitle("x[mm]");
  hist2->GetYaxis()->SetTitle("y[mm]");
  hist2->Draw("colz");
  
  c0->cd(1);
  VP2->Draw("colz");
  VP2->GetXaxis()->SetRangeUser(-300,300);
  VP2->GetYaxis()->SetRangeUser(-150,150);
  VP2->GetXaxis()->SetTitle("x[mm]");
  VP2->GetYaxis()->SetTitle("y[mm]");
  VP2->SetTitle("SS1 No Cut");
  
  c0->cd(2);
  VP4->Draw("colz");
  VP4->GetXaxis()->SetRangeUser(-300,300);
  VP4->GetYaxis()->SetRangeUser(-150,150);
  VP4->GetXaxis()->SetTitle("x[mm]");
  VP4->GetYaxis()->SetTitle("y[mm]");
  VP4->SetTitle("SS2 No Cut");
    
  TGraphErrors *grx2 = new TGraphErrors();  
  c1->cd(1);
  prjx2->SetTitle("SS1x");
  prjx2->GetXaxis()->SetTitle("x[mm]");
  prjx2->GetXaxis()->SetRangeUser(-250.,250.);
  prjx2->Draw();
  grx2->SetPoint(0,-400.,0.);
  grx2->SetPoint(1,400.,0.);
  std::cout << std::endl << "*prjx2 fitting result*" << std::endl;
  for (int i=-9; i<11; i++){
    std::string func_name = "gaus" + std::to_string(i);
    TF1 *gaus = new TF1(func_name.c_str(),"[0]*TMath::Gaus(-(x-[1])*(x-[1])/(2*[2]*[2]))",20.*i-10.-7.,20.*i-10.+7.);
    gaus->SetParameters(120.,20.*i-10.,3.);
    prjx2->Fit(gaus,"Q","",20.*i-10.-7.,20.*i-10.+7.);
    gaus->Draw("same");
    double position = gaus->GetParameter(1);
    double error = gaus->GetParError(1);
    double sigma = gaus->GetParameter(2);
    std::cout << position << " " << error << " " << sigma << std::endl;

    double x = 20.*i-10.;
    double y = x-position;
    double ey = error;
    grx2->SetPoint(i+11,x,y);
    grx2->SetPointError(i+11,0.,ey);
  }
  
  c1->cd(3);
  grx2->SetTitle("position accurancy of SS1x");
  grx2->GetXaxis()->SetRangeUser(-250.,250.);
  grx2->GetYaxis()->SetRangeUser(-1.,1.);
  grx2->GetXaxis()->SetTitle("x[mm]");
  grx2->GetYaxis()->SetTitle("#Delta x[mm]");
  grx2->SetMarkerStyle(20); 
  grx2->SetMarkerSize(1);
  grx2->SetMarkerColor(kBlue);
  grx2->Draw("AP");

  TGraphErrors *gry2 = new TGraphErrors();
  c1->cd(2);
  prjy2->SetTitle("SS1y");
  prjy2->GetXaxis()->SetTitle("y[mm]");
  prjy2->GetXaxis()->SetRangeUser(-120.,120.);
  prjy2->Draw();
  gry2->SetPoint(0,-400.,0.);
  gry2->SetPoint(1,400.,0.);
  std::cout << std::endl << "*prjy2 fitting result*" << std::endl;
  for(int j=-5; j<6; j++){
    std::string func_name = "gaus" + std::to_string(j);
    TF1 *gaus = new TF1(func_name.c_str(),"[0]*TMath::Gaus(-(x-[1])*(x-[1])/(2*[2]*[2]))",20.*j-7.,20.*j+7.);
    gaus->SetParameters(120.,20.*j,4.);
    prjy2->Fit(gaus,"Q","",20.*j-7.,20.*j+7.);
    gaus->Draw("same");
    double position = gaus->GetParameter(1);
    double error = gaus->GetParError(1);
    double sigma = gaus->GetParameter(2);
    std::cout << position << " " << error << " " << sigma << std::endl;
    
    double x = 20.*j;
    double y = x-position;
    double ey = error;
    gry2->SetPoint(j+7,x,y);
    gry2->SetPointError(j+7,0.,ey);
  }
  
  c1->cd(4);
  gry2->SetTitle("position accurancy of SS1y");
  gry2->GetXaxis()->SetRangeUser(-120.,120.);
  gry2->GetYaxis()->SetRangeUser(-1.,1.);
  gry2->GetXaxis()->SetTitle("y[mm]");
  gry2->GetYaxis()->SetTitle("#Delta y[mm]");
  gry2->SetMarkerStyle(20); 
  gry2->SetMarkerSize(1);
  gry2->SetMarkerColor(kBlue);
  gry2->Draw("AP");
  
  TGraphErrors *grx4 = new TGraphErrors();
  c2->cd(1);
  prjx4->SetTitle("VP4x");
  prjx4->GetXaxis()->SetTitle("x[mm]");
  prjx4->GetXaxis()->SetRangeUser(-250.,250.);
  prjx4->Draw();
  grx4->SetPoint(0,-400.,0.);
  grx4->SetPoint(1,400.,0.);
  std::cout << std::endl << "*prjx4 fitting result*" << std::endl;
  for (int i=-10; i<11; i++){
    std::string func_name = "gaus" + std::to_string(i);
    TF1 *gaus = new TF1(func_name.c_str(),"[0]*TMath::Gaus(-(x-[1])*(x-[1])/(2*[2]*[2]))",20.*i-7.,20.*i+7.);
    gaus->SetParameters(120.,20.*i,5.);
    prjx4->Fit(gaus,"Q","",20.*i-7.,20.*i+7.);
    gaus->Draw("same");
    double position = gaus->GetParameter(1);
    double error = gaus->GetParError(1);
    double sigma = gaus->GetParameter(2);
    std::cout << position << " " << error << " " << sigma << std::endl;

    double x = 20.*i;
    double y = x-position;
    double ey = error;
    grx4->SetPoint(i+12,x,y);
    grx4->SetPointError(i+12,0.,ey);
  }
  
  c2->cd(3);
  grx4->SetTitle("position accurancy of SS2x");
  grx4->GetXaxis()->SetRangeUser(-250.,250.);
  grx4->GetYaxis()->SetRangeUser(-1.,1.);
  grx4->GetXaxis()->SetTitle("x[mm]");
  grx4->GetYaxis()->SetTitle("#Delta x[mm]");
  grx4->SetMarkerStyle(20); 
  grx4->SetMarkerSize(1);
  grx4->SetMarkerColor(kBlue);
  grx4->Draw("AP");

  TGraphErrors *gry4 = new TGraphErrors();
  c2->cd(2);
  prjy4->SetTitle("SS2y");
  prjy4->GetXaxis()->SetTitle("y[mm]");
  prjy4->GetXaxis()->SetRangeUser(-120.,120.);
  prjy4->Draw();
  gry4->SetPoint(0,-400.,0.);
  gry4->SetPoint(1,400.,0.);
  std::cout << std::endl << "*prjy4 fitting result*" << std::endl;
  for(int j=-5; j<6; j++){
    std::string func_name = "gaus" + std::to_string(j);
    TF1 *gaus = new TF1(func_name.c_str(),"[0]*TMath::Gaus(-(x-[1])*(x-[1])/(2*[2]*[2]))",20.*j-7.,20.*j+7.);
    gaus->SetParameters(120.,20.*j,4.);
    prjy4->Fit(gaus,"Q","",20.*j-7.,20.*j+7.);
    gaus->Draw("same");
    double position = gaus->GetParameter(1);
    double error = gaus->GetParError(1);
    double sigma = gaus->GetParameter(2);
    std::cout << position << " " << error << " " << sigma << std::endl;

    double x = 20.*j;
    double y = x-position;
    double ey = error;
    gry4->SetPoint(j+7,x,y);
    gry4->SetPointError(j+7,0.,ey);
  }

  c2->cd(4);
  gry4->SetTitle("position accurancy of SS2y");
  gry4->GetYaxis()->SetRangeUser(-1.,1.);
  gry4->GetXaxis()->SetRangeUser(-120.,120.);
  gry4->GetXaxis()->SetTitle("y[mm]");
  gry4->GetYaxis()->SetTitle("#Delta y[mm]");
  gry4->SetMarkerStyle(20); 
  gry4->SetMarkerSize(1);
  gry4->SetMarkerColor(kBlue);
  gry4->Draw("AP");

  c0->SaveAs("Sieveslit_Profile_Gene8.png");
  c1->SaveAs("SS1_PositionAccuracy_Gene8.png");
  c2->SaveAs("SS2_PositionAccuracy_Gene8.png");
}
