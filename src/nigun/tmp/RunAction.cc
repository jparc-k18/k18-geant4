// ====================================================================
//   RunAction.cc
//
// ====================================================================
#include "G4Run.hh"
#include "RunAction.hh"

#include "TROOT.h"
#include "TH1.h"
#include "TFile.h"

//////////////////////
RunAction::RunAction()
//////////////////////
{
  // initialize root
  gROOT-> Reset();
  
  // define histograms
  hist_shower = new TH1D("shower", "Shower Development", 50, 0., 50.);
  hist_shower-> GetXaxis()-> SetTitle("Depth");
  hist_shower-> GetYaxis()-> SetTitle("Deposit Energy (MeV)");
}

///////////////////////
RunAction::~RunAction()
///////////////////////
{
  delete hist_shower;
}

///////////////////////////////////////////////////
void RunAction::BeginOfRunAction(const G4Run* /*aRun*/)
///////////////////////////////////////////////////
{
  hist_shower-> Reset();
}

/////////////////////////////////////////////////
void RunAction::EndOfRunAction(const G4Run* /*aRun*/)
/////////////////////////////////////////////////
{
  TFile* file = new TFile("outputs2s.root",
                          "RECREATE", "Geant4 ROOT analysis");
  hist_shower-> Write();
  
  file-> Close();
  delete file;
}

