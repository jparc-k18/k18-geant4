/*
  VPSD.cc
  2007/4  K.Shirotori
*/

#include "VPSD.hh"

#include "G4HCofThisEvent.hh"
#include "G4VPhysicalVolume.hh"
#include "G4TouchableHistory.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4VVisManager.hh"
#include "G4TouchableHandle.hh"
#include "G4SystemOfUnits.hh"

#include "ConfMan.hh"

namespace
{
const auto& confMan = ConfMan::GetInstance();
//const G4double DEThreshold = 0.001*keV;
const G4double DEThreshold = 0.00001*keV;
const double PositionSeparationThreshold = 2.0*cm;
const double TimeSeparationThreshold     = 5.0*ns;
}

//_____________________________________________________________________________
VPSD::VPSD( G4String name )
  : G4VSensitiveDetector(name),
    EMFlag(0)
{
  collectionName.insert(name);
  EMFlag = confMan.GetEMFlag();
}

//_____________________________________________________________________________
VPSD::~VPSD()
{
}

//_____________________________________________________________________________
void
VPSD::Initialize( G4HCofThisEvent *HCE )
{
  static int HCID = -1;
  VPCollection = new VPHitsCollection(SensitiveDetectorName,
                                      collectionName[0]);
  if( HCID<0 )
    HCID = GetCollectionID(0);
  HCE->AddHitsCollection( HCID, VPCollection );
}

//_____________________________________________________________________________
G4bool
VPSD::ProcessHits( G4Step *aStep, G4TouchableHistory *ROhist )
{
  //ConfMan *confMan = ConfMan::GetConfManager();
  //  G4double edep = aStep->GetTotalEnergyDeposit();

  //if( edep<=DEThreshold && EMFlag == 1) return true;

  G4Track *aTrack = aStep->GetTrack();
  G4StepPoint *preStepPoint = aStep->GetPreStepPoint();
  G4TouchableHandle theTouchable = preStepPoint->GetTouchableHandle();
  G4int hitLayer = theTouchable->GetReplicaNumber();
  //  G4int trackNo = aTrack->GetTrackID();
  //All Perticle
  G4double hittime = aTrack->GetGlobalTime();
  G4ThreeVector hitmom = aTrack->GetMomentum();
  //  G4int nHits = VPCollection->entries();
  G4ThreeVector hitpos = aStep->GetPreStepPoint()->GetPosition();
  G4ThreeVector hitposl = theTouchable->GetHistory()->
    GetTopTransform().TransformPoint( hitpos );

  //Decay Particle Tracking
  G4String decayName= aTrack->GetDefinition()->GetParticleName();
  G4double path=aTrack->GetTrackLength();
  //G4cout<< "pathlength=" << path<<G4endl;
  if(1){
  //  if(decayName=="kaon-"){
    //    G4cout<<"decayName=\""<<decayName<<"\""<<G4endl;

//     for( G4int i=0; i<nHits; ++i ){
//       VPHit *aHit = (*VPCollection)[i];
//       if( hitLayer==aHit->GetLayerID() ){
// 	G4double time = aHit->GetTime();
// 	G4double lposx = aHit->GetXLocal();
// 	G4double lposy = aHit->GetYLocal();
// 	//G4cout << hitposl.x() << " <==> " << lposx << ", " << hitposl.z() << "<==>" << lposy << G4endl;
// 	if( fabs(hittime-time)<=TimeSeparationThreshold  &&
// 	    fabs(hitposl.x()-lposx) <= PositionSeparationThreshold &&
// 	    fabs(hitposl.z()-lposy) <= PositionSeparationThreshold ){
// 	  //aHit->AddEdep( edep );
// 	  if( hittime<time ){
// 	    aHit->SetTime( hittime );
// 	    //          aHit->SetTrackNo( trackNo );
// 	    aHit->SetPos( hitpos );
// 	    aHit->SetMom( hitmom );
// 	  }
// 	  return true;
// 	}
//       }
//     }

    VPHit *aHit = new VPHit();
    aHit->SetPath( path );
    //    if(hitLayer==0) G4cout<<"path[0]="<<path<<G4endl;
    aHit->SetLayerID( hitLayer );//Need
    aHit->SetTime( hittime );//Need
    //  aHit->SetPos( hitpos );
    aHit->SetMom( hitmom );
    //  aHit->SetTrackNo( trackNo );
    aHit->SetLocalPos( hitposl.x(), hitposl.y() );
    aHit->SetDecayParticleName( decayName );
    VPCollection->insert( aHit );
  }
#if 0
  G4cout << "[VPSD] " << "Layer=" << hitLayer
    //	 << " edep=" << edep/keV << "keV"
	 << " G: " << hitpos << "  L: " << hitposl
	 <<" P: "<< hitmom << G4endl;
#endif

  return true;
}

//_____________________________________________________________________________
void VPSD::EndOfEvent( G4HCofThisEvent *HCE )
{
}

//_____________________________________________________________________________
void VPSD::clear()
{
  G4int nHits = VPCollection->entries();
  for( G4int i=nHits-1; i>=0; --i )
    delete (*VPCollection)[i];
}

// void VPSD::DrawAll() const
// {
//   G4VVisManager *pVisManager = G4VVisManager::GetConcreteInstance();

//   if( pVisManager ){
//     G4int nHits = VPCollection->entries();
//     for( G4int i=0; i<nHits; ++i )
//       (*VPCollection)[i]->Draw();
//   }
// }

// void VPSD::PrintAll() const
// {
//   G4int nHits = VPCollection->entries();
//   for( G4int i=0; i<nHits; ++i)
//     (*VPCollection)[i]->Print();
// }
