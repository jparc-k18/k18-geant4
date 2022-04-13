// ====================================================================
//   EventAction.cc
//
// ====================================================================
#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4SDManager.hh"

#include "EventAction.hh"
#include "CalorimeterSD.hh"

#include "TROOT.h"
#include "TH1.h"

//////////////////////////
EventAction::EventAction()
//////////////////////////
{
}

///////////////////////////
EventAction::~EventAction()
///////////////////////////
{
}

////////////////////////////////////////////////////////////
void EventAction::BeginOfEventAction(const G4Event* /*anEvent*/)
////////////////////////////////////////////////////////////
{
}

//////////////////////////////////////////////////////////
void EventAction::EndOfEventAction(const G4Event* anEvent)
//////////////////////////////////////////////////////////
{
  G4SDManager* SDManager= G4SDManager::GetSDMpointer();

  // get SD
  //CalorimeterSD* calSD= (CalorimeterSD*) 
  //SDManager-> FindSensitiveDetector("/calorimeter");
  //calSD-> PrintAll();

  // get "Hit Collection of This Event"
  G4HCofThisEvent* HCTE= anEvent-> GetHCofThisEvent();
  if(! HCTE) return;

  // get a hit collection
  static G4int idcal= -1;
  if(idcal<0)  idcal= SDManager-> GetCollectionID("calorimeter");
  CalHitsCollection* hccal= (CalHitsCollection*)HCTE-> GetHC(idcal);
  if (!hccal) return; // no hit collection

  // get hits
  G4int nhits= hccal-> entries();
  for(G4int idx=0; idx< nhits; idx++) {
    G4int ich= (*hccal)[idx]-> GetID();
    G4double edep= (*hccal)[idx]-> GetEdep();

    // fill a histogram
    TH1D* hist_shower= (TH1D*)gROOT-> FindObject("shower");
    hist_shower-> Fill(ich, edep/MeV);
  }
}

