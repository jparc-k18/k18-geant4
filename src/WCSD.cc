// -*- C++ -*-

#include "WCSD.hh"

#include <G4HCofThisEvent.hh>
#include <G4VPhysicalVolume.hh>
#include <G4TouchableHistory.hh>
#include <G4Step.hh>
#include <G4Track.hh>
#include <G4VVisManager.hh>
#include <G4TouchableHandle.hh>
#include <G4SystemOfUnits.hh>
#include <Randomize.hh>

#include <TF1.h>
#include <TGraph2D.h>
#include <TMath.h>

#include "ConfMan.hh"

namespace
{
const auto& confMan = ConfMan::GetInstance();
//const G4double DEThreshold = 0.001*keV;
// const G4double DEThreshold = 0.00001*keV;
// const double PositionSeparationThreshold = 2.0*cm;
// const double TimeSeparationThreshold     = 5.0*ns;
const G4double Win = 1.33;
const G4double ThresholdVelocity = 1./Win;
//const G4double lmin = 300./1000000000.; // wave length at minimum
//const G4double lmax = 600./1000000000.; // wave length at maximum
const G4double lmin = 300.0*nm; // wave length at minimum
const G4double lmax = 600.0*nm; // wave length at maximum
//const G4double CE = 0.30; // light collection efficiency of a box;
const G4double QE = 0.30; // quantum efficiency of a PMT
//const G4double npefactor = 1./22.0/2.0; // to be ~30 for each PMT (Toshi, 26Nov2014)
//const G4double npefactor = 1./22.0/2.0 * 0.89; // to be ~40 for pi+ in each PMT (Toshi, 16Dec2014)
//const G4double npefactor = 1./22.0/2.0 * 1.18; // to be ~40 for pi+ in each PMT (Toshi, 22Dec2014)
//const G4double npefactor = 1./22.0/2.0 * 1.18 * 2.5; // to be ~100 for pi+ in each PMT (Toshi, 25Dec2014)
//const G4double npefactor = 1./22.0/2.0 * 1.18 * 2.5 * 1.05; // changed for test (Toshi, 27Jan2015)
//const G4double npefactor = 0.033 * 2.0; // to be ~100 NPE for each PMT (T.Gogami, 9July2015) Sum = 200
const G4double npefactor = 0.0349 * 2.0; // to be ~100 NPE for each PMT (T.Gogami, 9July2015) Sum = 200
//const G4double npefactor = 0.0349 * 2.0 / 2.0; // to be ~100/2.0 NPE for each PMT (T.Gogami, 24July2015) Sum = 200/2.0
//const G4double spl = 2.99792458; // m/s
}

//_____________________________________________________________________________
WCSD::WCSD(const G4String& name)
  : G4VSensitiveDetector(name)
{
  collectionName.insert(name);
  // ==== Define X vs. Y vs. NPE table which was ===== //
  // ==== obtained in ELPH experiment in 2014    ===== //
  DefineXYTable();                  // =============== //
  // ================================================= //
  // ==== Define Y vs. NPE table which was   ========= //
  // ==== obtained in cosmic-ray test in 2015    ===== //
  DefineYTable();                   // =============== //
  // ================================================= //
}

//_____________________________________________________________________________
WCSD::~WCSD()
{
}

//_____________________________________________________________________________
void
WCSD::DefineXYTable()
{
  //const int n_wcdat=21;
  const int n_wcdat=20;
  double x_wcdat[n_wcdat]={0, 4, 0, 0, 0,
			   0, 0, 0, 0, 0,
			   9, 9, 9, 9, 9,
			   //9, 9, 9, 0, 4,
			   9, 9, 9, 4,
			   4 };
  double y_wcdat[n_wcdat]={0, 0, 10, 20, 23,
			   25, 27, 29, 31, 33,
			   0, 10, 20, 23, 25,
			   //27, 29, 31, 0, 27,
			   27, 29, 31, 27,
			   30 };
  double z_wcdat[n_wcdat]={1, 1.02637, 0.993974, 1.00959, 1.03243,
			   1.01257, 1.01404, 1.0233, 1.03977, 1.07227,
			   0.98562, 0.978638, 0.973557, 0.986008, 0.99615,
			   //1.00098, 0.991164, 0.953529, 1, 1.02017,
			   1.00098, 0.991164, 0.953529, 1.02017,
			   1.08848 };
  xyTable = new TGraph2D(n_wcdat,x_wcdat,y_wcdat,z_wcdat);
}

//_____________________________________________________________________________
void
WCSD::DefineYTable()
{
  yTable = new TF1("func1","[0]+[1]*x+[2]*x*x",-40.0,40.0);
  //yTable = new TF1("func1","[0]+[1]*x+[2]*x*x");

  // ====== These parameters were obtained by cosmic-ray test =====
  yTable->SetParameters(1.0,-1.22502e-4,1.42488e-4);
}

//_____________________________________________________________________________
void
WCSD::Initialize(G4HCofThisEvent* HCE)
{
  WCCollection = new WCHitsCollection(SensitiveDetectorName,
                                      collectionName[0]);
  HCE->AddHitsCollection(GetCollectionID(0), WCCollection);
}

//_____________________________________________________________________________
G4bool
WCSD::ProcessHits(G4Step* aStep, G4TouchableHistory* /* ROhist */)
{
#if 0
  //  ConfMan *confMan = ConfMan::GetConfManager();
  G4double edep = aStep->GetTotalEnergyDeposit();

  //if( edep<=DEThreshold && EMFlag == 1) return true; // Comment out, 25Nov2014

  G4Track* aTrack = aStep->GetTrack();
  G4StepPoint* preStepPoint = aStep->GetPreStepPoint();
  G4TouchableHandle theTouchable = preStepPoint->GetTouchableHandle();
  G4int hitLayer = theTouchable->GetReplicaNumber();
  //  G4int trackNo = aTrack->GetTrackID();

  // ~~~~~~ All Perticle ~~~~~~~~
  G4double hittime = aTrack->GetGlobalTime();
  G4ThreeVector hitmom = aTrack->GetMomentum();
  G4int nHits = WCCollection->entries();
  //G4int nHits = WCCollection->Getentries();
  G4ThreeVector hitpos = aStep->GetPreStepPoint()->GetPosition();
  G4ThreeVector hitposl = theTouchable->GetHistory()->
    GetTopTransform().TransformPoint( hitpos );
  G4double velocity = aTrack-> GetVelocity()/CLHEP::c_light;
  //G4double velocity = aTrack-> GetVelocity()*m/ns;
  //velocity = velocity/spl;
  G4double npe = 0.0;
  G4double StepLength = 0.;
  StepLength = aStep->GetStepLength();
  G4double Charge=aStep->GetTrack()->GetDefinition()->GetPDGCharge();
  G4double thBeta = ThresholdVelocity;
  G4double tempx, tempy, tempz;
  G4double xyFactor;
  tempx = hitposl.x() /cm; // Beam direction
  tempy = hitposl.y() /cm; // Horizontal direction
  tempz = hitposl.z() /cm; // Vertical direction
  if(tempx<0) tempx = -1.0 * tempx;
  if(tempy<0) tempy = -1.0 * tempy;
  if(tempz<0) tempz = -1.0 * tempz;
  //G4cout << nHits << G4endl;
  //G4cout << nHits << " " << hitposl.x()/cm << G4endl;
  //  G4cout << " " << tempx << "," << tempy << " --> "
  //	 << GetXYdepFactor(tempx,tempy) << G4endl;

  // -------------------------------------------------------------
  //xyFactor = GetXYdepFactor(tempy,tempz); // XY table from ELPH beam test (old prototype box)
  xyFactor = GetYdepFactor(hitposl.z()/cm); // Y table from cosmic-ray test (new prototype box)
  // -------------------------------------------------------------
  //G4cout << " " << tempx << " " << xyFactor << G4endl;
//  G4cout << nHits  << ": " << hitposl.x()/cm << " " << hitposl.y()/cm << " "
//	 << hitposl.z()/cm << " --> " << xyFactor << G4endl;

  if( velocity >= thBeta && Charge!=0 ) {
  //if( velocity >= thBeta ) {
    //	  npe += 2.*pi/137.*(1.-1./1.33/(thBeta)/1.33/(thBeta))*
    //	    (StepLength/cm)*(1/0.0000003-1/0.0000006)*CE*QE;
    npe = 2.0*CLHEP::pi/137.*(1.-1./(Win*Win*velocity*velocity))*
      //(StepLength)*(1./lmin-1./lmax)*CE*QE;
      (StepLength)*(1./lmin-1./lmax)*npefactor*QE;
      //(StepLength/m)*(1./lmin-1./lmax)*CE*QE;
    //npe = npe + (npe * (G4UniformRand()-0.5)*2.0 * 0.2); // npe +/- 20%
    //npe = npe + (npe * (G4UniformRand()-0.5)*2.0 * 1.0); // npe +/- 100%
      //npe = npe + (random+/-20%);
  }
  else{
    npe = 0.0;
  }

  //~~~~~~~~ Decay Particle Tracking ~~~~~~~~
  G4String decayName= aTrack->GetDefinition()->GetParticleName();
  WCHit *aHit = new WCHit();
  aHit->SetPass();
  aHit->SetLayerID( hitLayer );
  aHit->SetTime( hittime );
  aHit->SetPos( hitpos );
  aHit->SetMom( hitmom );
  //aHit->SetNPE( npe ); //
  aHit->SetNPE( npe*xyFactor ); //
  //aHit->AddNPE( npe );
  //aHit->SetTrackNo( trackNo );
  aHit->SetLocalPos( hitposl.x(), hitposl.z() );
  aHit->SetDecayParticleName( decayName );
  aHit->SetEdep( edep );
  WCCollection->insert( aHit );
  //G4cout << aHit->GetNPE() << G4endl;

#if 0
  G4cout << "[WCSD] " << "Layer=" << hitLayer
	 << " edep=" << edep/keV << "keV"
	 << " G: " << hitpos << "  L: " << hitposl
	 <<" P: "<< hitmom << G4endl;
#endif

#endif

  const auto preStepPoint = aStep->GetPreStepPoint();
  const auto aTrack = aStep->GetTrack();
  const auto Definition = aTrack->GetDefinition();
  const G4String particleName = Definition->GetParticleName();
  const G4String particleType = Definition->GetParticleType();
  if(preStepPoint->GetStepStatus() != fGeomBoundary)
    return false;
  if(Definition->GetPDGCharge() == 0.)
    return false;
  WCCollection->insert(new WCHit(SensitiveDetectorName, aStep));
  return true;
}

//_____________________________________________________________________________
void WCSD::EndOfEvent(G4HCofThisEvent* /* HCE */)
{
}
