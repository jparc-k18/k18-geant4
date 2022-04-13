/*
  "S2SRunAction.cc"
  
  T.Gogami, 10 Nov 2014
*/

#include "S2SRunAction.hh"
//#include "S2SParamManager.hh"
#include "S2SAnalysis.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"

#include <TFile.h>
#include <TTree.h>
#include <TObjArray.h>
#include <TH1F.h>
#include <stdlib.h>

S2SRunAction::S2SRunAction()
//:RfileName("../analysis/root/DRAGON.root")
{}

S2SRunAction::S2SRunAction(S2SAnalysis* ana)
  ://RfileName("../analysis/root/DRAGON.root"),
  analysis(ana)
{}

S2SRunAction::~S2SRunAction()
{}

//S2SRunAction::S2SRunAction(S2SParamManager* PManDragon)
//  :paramMan(PManDragon),RfileName("DRAGON.root")
//{
//  RfileName = paramMan->GetROOTFileName();
//  if( paramMan->GetRFnum()!=0 ){
//    char tmp[80];
//    sprintf(tmp,"%d",paramMan->GetRFnum());
//    RfileName = RfileName+"_"+tmp;
//
//  }
//}

void S2SRunAction::BeginOfRunAction(const G4Run* aRun)
{
  //G4cout << " " << G4endl;
  //G4cout << "// Simulation Start //" << G4endl;
  //G4cout << " " << G4endl;
  //file = new TFile(RfileName,"recreate"); //Open ROOT file (Recreate)
  //TTree *tree = new TTree("tree","S-2S simulation");//Tree
  //Tree = (TTree*)file->Get("tree");
  analysis->BeginOfRun( aRun );
}

void S2SRunAction::EndOfRunAction(const G4Run* aRun)
{
  //Tree->Write();  //Write tree into ROOT file
  //file -> Close();//Close ROOT file
  analysis->EndOfRun( aRun );
  G4cout << "  END... bye ! " << G4endl;
}
