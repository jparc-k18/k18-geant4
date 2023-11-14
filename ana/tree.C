// -*- C++ -*-

#include "DetectorID.hh"

void tree()
{
  if(!gFile || !gFile->IsOpen()){
    std::cout << "Usage: " << std::endl
              << "   $ root foo.root" << std::endl
              << "   root [0] .x " << __func__ << ".C" << std::endl;
    gSystem->Exit(EXIT_FAILURE);
  }

  TTreeReader reader("g4s2s", gFile);
  TTreeReaderValue<Int_t> evnum(reader, "evnum");
  TTreeReaderValue<std::vector<TParticle>> PRM(reader, "PRM");
  TTreeReaderValue<std::vector<TParticle>> SDC(reader, "SDC");
  TTreeReaderValue<std::vector<TParticle>> TOF(reader, "TOF");
  TTreeReaderValue<std::vector<TParticle>> AC1(reader, "AC1");
  TTreeReaderValue<std::vector<TParticle>> WC(reader, "WC");
  TTreeReaderValue<std::vector<TParticle>> VP(reader, "VP");

  auto hTOFHitPat = new TH1D("hTOFHitPat", "TOF HitPat",
                             NumOfSegTOF, 0.5, NumOfSegTOF+0.5);
  auto hTOFMulti = new TH1D("hTOFMulti", "TOF Multiplicity",
                            NumOfSegTOF+1, -0.5, NumOfSegTOF+0.5);

  while(reader.Next()){
    Int_t tof_multiplicity = 0;
    for(const auto& p : (*TOF)){
      hTOFHitPat->Fill(p.GetMother(1)); /* segment */
      if(true
         // && p.GetWeight() > 0.01 /* energy deposit */
         && p.GetPDG()->Charge() != 0
         ){
        tof_multiplicity++;
      }
    }
    hTOFMulti->Fill(tof_multiplicity);
  }

  auto c1 = new TCanvas("c1", "c1", 800, 400);
  c1->Divide(2, 1);
  c1->cd(1); hTOFMulti->Draw();
  c1->cd(2); hTOFHitPat->Draw();
}
