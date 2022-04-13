/*
  SlSD.cc
  2007/4  K.Shirotori
*/

#include "SlSD.hh"

#include "G4HCofThisEvent.hh"
#include "G4VPhysicalVolume.hh"
#include "G4TouchableHistory.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4VVisManager.hh"
#include "G4TouchableHandle.hh"
#include "G4SystemOfUnits.hh"

#include "ConfMan.hh"

//const G4double DEThreshold = 0.001*keV;
const G4double DEThreshold = 0.00001*keV;
const double PositionSeparationThreshold = 2.0*cm;
const double TimeSeparationThreshold     = 5.0*ns;

SlSD::SlSD( G4String name )
  : G4VSensitiveDetector(name), EMFlag(0)
{
  collectionName.insert( name/*G4String( "SlCollection" )*/ );

  ConfMan *confMan = ConfMan::GetConfManager();
  EMFlag = confMan->GetEMFlag();
}

SlSD::~SlSD()
{
}

void SlSD::Initialize( G4HCofThisEvent *HCE )
{
  static int HCID = -1;
  SlCollection =
    new SlHitsCollection( SensitiveDetectorName,
			     collectionName[0] );
  if( HCID<0 )
    HCID = GetCollectionID(0);

  HCE->AddHitsCollection( HCID, SlCollection );

}

G4bool SlSD::ProcessHits( G4Step *aStep,
			     G4TouchableHistory *ROhist )
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
  //  G4int nHits = SlCollection->entries();
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
//       SlHit *aHit = (*SlCollection)[i];
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

    SlHit *aHit = new SlHit();
    aHit->SetPath( path );
    //    if(hitLayer==0) G4cout<<"path[0]="<<path<<G4endl;
    aHit->SetLayerID( hitLayer );//Need
    aHit->SetTime( hittime );//Need
    //  aHit->SetPos( hitpos );
    aHit->SetMom( hitmom );
    //  aHit->SetTrackNo( trackNo );
    aHit->SetLocalPos( hitposl.x(), hitposl.y() );
    aHit->SetDecayParticleName( decayName );
    SlCollection->insert( aHit );
  }
#if 0
  G4cout << "[SlSD] " << "Layer=" << hitLayer 
    //	 << " edep=" << edep/keV << "keV"  
	 << " G: " << hitpos << "  L: " << hitposl 
	 <<" P: "<< hitmom << G4endl;
#endif

  return true;
}

void SlSD::EndOfEvent( G4HCofThisEvent *HCE )
{
}

void SlSD::clear()
{
  G4int nHits = SlCollection->entries();
  for( G4int i=nHits-1; i>=0; --i )
    delete (*SlCollection)[i];
}

// void SlSD::DrawAll() const
// {
//   G4VVisManager *pVisManager = G4VVisManager::GetConcreteInstance();

//   if( pVisManager ){
//     G4int nHits = SlCollection->entries();
//     for( G4int i=0; i<nHits; ++i )
//       (*SlCollection)[i]->Draw();
//   }
// }

// void SlSD::PrintAll() const
// {
//   G4int nHits = SlCollection->entries();
//   for( G4int i=0; i<nHits; ++i)
//     (*SlCollection)[i]->Print();
// }

