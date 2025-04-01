#include <TFile.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>

//_____________________________________________________________________________
void
Check_Q1effect2( void )
{
  gStyle->SetOptStat(0);

  //  TFile::Open( "Gene6_Q1_100_Q2_100.root" );
  TFile *f1 = new TFile( "Gene6_Q1_100_Q2_100.root" );
  //TTreeReader reader( "g4s2s", gFile );
  TTreeReader reader( "g4s2s", f1 );
  TTreeReaderValue<Double_t> theta0( reader, "theta0" );
  TTreeReaderValue<Double_t> p0( reader, "p0" );
  TTreeReaderValue<std::vector<TParticle>> TOF( reader, "TOF" );
  TTreeReaderValue<std::vector<TParticle>> SDC1( reader, "SDC1" );
  TTreeReaderValue<std::vector<TParticle>> VP1( reader, "VP" );

  TFile *f2 = new TFile( "Gene6_Q1_101_Q2_100.root" );
  TFile *f3 = new TFile( "Gene6_Q1_102_Q2_100.root" );
  TFile *f4 = new TFile( "Gene6_Q1_104_Q2_100.root" );
  TFile *f5 = new TFile( "Gene6_Q1_108_Q2_100.root" );
  //TTreeReader reader( "g4s2s", gFile );
  TTreeReader reader2( "g4s2s", f2 );
  TTreeReader reader3( "g4s2s", f3 );
  TTreeReader reader4( "g4s2s", f4 );
  TTreeReader reader5( "g4s2s", f5 );    
  TTreeReaderValue<std::vector<TParticle>> VP2( reader2, "VP" );
  TTreeReaderValue<std::vector<TParticle>> VP3( reader3, "VP" );
  TTreeReaderValue<std::vector<TParticle>> VP4( reader4, "VP" );
  TTreeReaderValue<std::vector<TParticle>> VP5( reader5, "VP" );

  TH2D * hVP2x[4];
  TH2D * hVP2y[4];
  TH2D * hVP4x[4];
  TH2D * hVP4y[4];
  for(int i=0; i<4; ++i){
    hVP2x[i] = new TH2D(Form("hVP2x_%d",i),
			Form("hVP2x_%d",i),
			200,-400,400,100,-20,20);
    hVP2y[i] = new TH2D(Form("hVP2y_%d",i),
			Form("hVP2y_%d",i),
			150,-150,150,100,-20,20);
    hVP4x[i] = new TH2D(Form("hVP4x_%d",i),
			Form("hVP4x_%d",i),
			200,-600,600,100,-40,40);
    hVP4y[i] = new TH2D(Form("hVP4y_%d",i),
			Form("hVP4y_%d",i),
			150,-100,100,100,-20,20);
    hVP2x[i]->SetLineColor(i+1);
    hVP2x[i]->SetMarkerColor(i+1);
    hVP2y[i]->SetLineColor(i+1);
    hVP2y[i]->SetMarkerColor(i+1);

    hVP4x[i]->SetLineColor(i+1);
    hVP4x[i]->SetMarkerColor(i+1);
    hVP4y[i]->SetLineColor(i+1);
    hVP4y[i]->SetMarkerColor(i+1);
  }

  
  while( reader.Next() ){
    Int_t ievent = reader.GetCurrentEntry();
    
    int tof_multi=0;
    for(const auto& p : (*TOF)){
      if(true
	 // && p.GetWeight() > 0.01 /* energy deposit */
	 && p.GetPDG()->Charge() != 0
	 ){
	tof_multi++;
      }
    }
    int SDC1_multi=0;
    for(const auto& p : (*SDC1)){
      if(true
	 // && p.GetWeight() > 0.01 /* energy deposit */
	 && p.GetPDG()->Charge() != 0
	 ){
	SDC1_multi++;
      }
    }
    // std::cout<<"tof_multi:"<<tof_multi
    // 	     <<", SDC1_multi"<<SDC1_multi<<std::endl;

    if(tof_multi>0&&SDC1_multi>0){
      reader2.SetEntry(ievent);
      reader3.SetEntry(ievent);
      reader4.SetEntry(ievent);
      reader5.SetEntry(ievent);

      double VP2x[5]={0}, VP2y[5]={0};
      double VP4x[5]={0}, VP4y[5]={0};
      for(const auto& p : (*VP1)){
	if(true
	   && p.GetPDG()->Charge() != 0
	   ){
	  int layer = p.GetMother(1);
	  if(layer==2){
	    VP2x[0]=p.Vx();
	    VP2y[0]=p.Vy();
	    // std::cout<<"VP2x:"<<VP2x[0]<<", y:"<<VP2y[0]<<std::endl;
	  }
	  if(layer==4){
	    VP4x[0]=p.Vx();
	    VP4y[0]=p.Vy();
	    //std::cout<<"VP4x:"<<VP4x[0]<<", y:"<<VP4y[0]<<std::endl;
	  }
	}
      }
      for(const auto& p : (*VP2)){
	if(true
	   && p.GetPDG()->Charge() != 0
	   ){
	  int layer = p.GetMother(1);
	  if(layer==2){
	    VP2x[1]=p.Vx();
	    VP2y[1]=p.Vy();
	    // std::cout<<"VP2x:"<<VP2x[0]<<", y:"<<VP2y[0]<<std::endl;
	  }
	  if(layer==4){
	    VP4x[1]=p.Vx();
	    VP4y[1]=p.Vy();
	    //std::cout<<"VP4x:"<<VP4x[0]<<", y:"<<VP4y[0]<<std::endl;
	  }
	}
      }
      for(const auto& p : (*VP3)){
	if(true
	   && p.GetPDG()->Charge() != 0
	   ){
	  int layer = p.GetMother(1);
	  if(layer==2){
	    VP2x[2]=p.Vx();
	    VP2y[2]=p.Vy();
	    // std::cout<<"VP2x:"<<VP2x[0]<<", y:"<<VP2y[0]<<std::endl;
	  }
	  if(layer==4){
	    VP4x[2]=p.Vx();
	    VP4y[2]=p.Vy();
	    //std::cout<<"VP4x:"<<VP4x[0]<<", y:"<<VP4y[0]<<std::endl;
	  }
	}
      }
      for(const auto& p : (*VP4)){
	if(true
	   && p.GetPDG()->Charge() != 0
	   ){
	  int layer = p.GetMother(1);
	  if(layer==2){
	    VP2x[3]=p.Vx();
	    VP2y[3]=p.Vy();
	    // std::cout<<"VP2x:"<<VP2x[0]<<", y:"<<VP2y[0]<<std::endl;
	  }
	  if(layer==4){
	    VP4x[3]=p.Vx();
	    VP4y[3]=p.Vy();
	    //std::cout<<"VP4x:"<<VP4x[0]<<", y:"<<VP4y[0]<<std::endl;
	  }
	}
      }
      for(const auto& p : (*VP5)){
	if(true
	   && p.GetPDG()->Charge() != 0
	   ){
	  int layer = p.GetMother(1);
	  if(layer==2){
	    VP2x[4]=p.Vx();
	    VP2y[4]=p.Vy();
	    // std::cout<<"VP2x:"<<VP2x[0]<<", y:"<<VP2y[0]<<std::endl;
	  }
	  if(layer==4){
	    VP4x[4]=p.Vx();
	    VP4y[4]=p.Vy();
	    //std::cout<<"VP4x:"<<VP4x[0]<<", y:"<<VP4y[0]<<std::endl;
	  }
	}
      }

      if(VP2x[0]!=0.&&VP2x[1]!=0.&&VP2x[2]!=0.&&VP2x[3]!=0.&&VP2x[4]!=0.){
	std::cout<<"VP2x:"<<VP2x[0]
		 <<", "<<VP2x[1]
		 <<", "<<VP2x[2]
		 <<", "<<VP2x[3]
		 <<", "<<VP2x[4]<<std::endl;
	for(int i=0; i<4; ++i){
	  hVP2x[i]->Fill(VP2x[0],VP2x[0]-VP2x[i+1]);
	  hVP2y[i]->Fill(VP2y[0],VP2y[0]-VP2y[i+1]);
	}
      }

      if(VP4x[0]!=0.&&VP4x[1]!=0.&&VP4x[2]!=0.&&VP4x[3]!=0.&&VP4x[4]!=0.){
	std::cout<<"VP4x:"<<VP4x[0]
		 <<", "<<VP4x[1]
		 <<", "<<VP4x[2]
		 <<", "<<VP4x[3]
		 <<", "<<VP4x[4]<<std::endl;
	for(int i=0; i<4; ++i){
	  hVP4x[i]->Fill(VP4y[0],VP4y[0]-VP4y[i+1]);
	  hVP4y[i]->Fill(VP4x[0],VP4x[0]-VP4x[i+1]);
	}
      }
    }
  }
  TCanvas *c0 = new TCanvas("c0","c0");
  c0->Divide(2,2);
  c0->cd(1)->SetGrid(1);
  for(int i=0; i<4; ++i){
    hVP2x[3-i]->SetTitle("VP2 X#DeltaY");
    hVP2x[3-i]->GetXaxis()->SetTitle("X [mm]");
    hVP2x[3-i]->GetYaxis()->SetTitle("#DeltaY [mm]");
    if(i==0)
      hVP2x[3-i]->Draw("scat");
    else
      hVP2x[3-i]->Draw("scat&&same");
  }

  c0->cd(2)->SetGrid(1);
  for(int i=0; i<4; ++i){
    hVP2y[3-i]->SetTitle("VP2 Y#DeltaX");
    hVP2y[3-i]->GetXaxis()->SetTitle("Y [mm]");
    hVP2y[3-i]->GetYaxis()->SetTitle("#DeltaX [mm]");
    if(i==0)
      hVP2y[3-i]->Draw("scat");
    else
      hVP2y[3-i]->Draw("scat&&same");
  }

  c0->cd(3)->SetGrid(1);
  for(int i=0; i<4; ++i){
    hVP4x[3-i]->SetTitle("VP4 X#DeltaY");
    hVP4x[3-i]->GetXaxis()->SetTitle("X [mm]");
    hVP4x[3-i]->GetYaxis()->SetTitle("#DeltaY [mm]");
    if(i==0)
      hVP4x[3-i]->Draw("scat");
    else
      hVP4x[3-i]->Draw("scat&&same");
  }
  c0->cd(4)->SetGrid(1);
  for(int i=0; i<4; ++i){
    hVP4y[3-i]->SetTitle("VP4 Y#DeltaX");
    hVP4y[3-i]->GetXaxis()->SetTitle("Y [mm]");
    hVP4y[3-i]->GetYaxis()->SetTitle("#DeltaX [mm]");
    if(i==0)
      hVP4y[3-i]->Draw("scat");
    else
      hVP4y[3-i]->Draw("scat&&same");
  }
}
