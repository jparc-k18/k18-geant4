/*
  TOFSD.cc
  2007/4  K.Shirotori
*/

#include "ACSD.hh"

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

ACSD::ACSD( G4String name )
  : G4VSensitiveDetector(name)
{
  collectionName.insert( name/*G4String( "ACCollection" )*/ );

  const auto& confMan = ConfMan::GetInstance();
}

ACSD::~ACSD()
{
}

void ACSD::Initialize( G4HCofThisEvent *HCE )
{
  static int HCID = -1;
  ACCollection =
    new ACHitsCollection( SensitiveDetectorName,
			     collectionName[0] );
  if( HCID<0 )
    HCID = GetCollectionID(0);

  HCE->AddHitsCollection( HCID, ACCollection );
}

G4bool ACSD::ProcessHits( G4Step *aStep,
			     G4TouchableHistory *ROhist )
{
  //  ConfMan *confMan = ConfMan::GetConfManager();
  G4double edep = aStep->GetTotalEnergyDeposit();

  //if( edep<=DEThreshold && EMFlag == 1) return true;

  G4Track *aTrack = aStep->GetTrack();
  G4StepPoint *preStepPoint = aStep->GetPreStepPoint();
  G4TouchableHandle theTouchable = preStepPoint->GetTouchableHandle();
  G4int hitLayer = theTouchable->GetReplicaNumber();
  //  G4int trackNo = aTrack->GetTrackID();
  //All Perticle
  G4double hittime = aTrack->GetGlobalTime();
  G4ThreeVector hitmom = aTrack->GetMomentum();
  G4int nHits = ACCollection->entries();
  G4ThreeVector hitpos = aStep->GetPreStepPoint()->GetPosition();
  G4ThreeVector hitposl = theTouchable->GetHistory()->
    GetTopTransform().TransformPoint( hitpos );

  //Decay Particle Tracking
  G4String decayName= aTrack->GetDefinition()->GetParticleName();
  //  G4double path=aTrack->GetTrackLength();
  //G4cout<< "pathlength=" << path<<G4endl;

  for( G4int i=0; i<nHits; ++i ){
    ACHit *aHit = (*ACCollection)[i];
    if( hitLayer==aHit->GetLayerID() ){
      G4double time = aHit->GetTime();
      G4double lposx = aHit->GetXLocal();
      G4double lposy = aHit->GetYLocal();
      //G4cout << hitposl.x() << " <==> " << lposx << ", " << hitposl.z() << "<==>" << lposy << G4endl;
      if( fabs(hittime-time)<=TimeSeparationThreshold  &&
	  fabs(hitposl.x()-lposx) <= PositionSeparationThreshold &&
	  fabs(hitposl.z()-lposy) <= PositionSeparationThreshold ){
	aHit->AddEdep( edep );
	if( hittime<time ){
          aHit->SetTime( hittime );
	  //          aHit->SetTrackNo( trackNo );
          aHit->SetPos( hitpos );
          aHit->SetMom( hitmom );
        }
        return true;
      }
    }
  }

  ACHit *aHit = new ACHit();
  aHit->SetPass();
  aHit->SetLayerID( hitLayer );
  aHit->SetTime( hittime );
  //  aHit->SetPos( hitpos );
  aHit->SetMom( hitmom );
  //  aHit->SetTrackNo( trackNo );
  aHit->SetLocalPos( hitposl.y(), hitposl.z() );
  aHit->SetDecayParticleName( decayName );
  ACCollection->insert( aHit );

#if 0
  G4cout << "[ACSD] " << "Layer=" << hitLayer
	 << " edep=" << edep/keV << "keV"
	 << " G: " << hitpos << "  L: " << hitposl
	 <<" P: "<< hitmom << G4endl;
#endif

  return true;
}

void ACSD::EndOfEvent( G4HCofThisEvent *HCE )
{
}

void ACSD::clear()
{
  G4int nHits = ACCollection->entries();
  for( G4int i=nHits-1; i>=0; --i )
    delete (*ACCollection)[i];
}

// void ACSD::DrawAll() const
// {
//   G4VVisManager *pVisManager = G4VVisManager::GetConcreteInstance();

//   if( pVisManager ){
//     G4int nHits = ACCollection->entries();
//     for( G4int i=0; i<nHits; ++i )
//       (*ACCollection)[i]->Draw();
//   }
// }

// void ACSD::PrintAll() const
// {
//   G4int nHits = ACCollection->entries();
//   for( G4int i=0; i<nHits; ++i)
//     (*ACCollection)[i]->Print();
// }
