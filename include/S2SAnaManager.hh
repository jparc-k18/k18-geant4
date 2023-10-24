// -*- C++ -*-

#ifndef S2SAnaManager_h
#define S2SAnaManager_h 1

#include <fstream>
#include <map>
#include <vector>

#include <globals.hh>
#include <G4ThreeVector.hh>
#include <G4String.hh>

#include <TParticle.h>
#include <TRandom3.h>
#include <TVector3.h>

#include "RootHelper.hh"

class G4Run;
class G4Event;
class PrimaryInfo;

class TFile;
class TTree;

class VHitInfo;

#ifndef NumDC
#define NumDC 30
#endif
#ifndef NumTOFSeg
#define NumTOFSeg 18
#endif
#ifndef NumWCSeg
#define NumWCSeg 16
#endif
#ifndef NumDC
#define NumDC 4
#endif
#ifndef MaxHits1
#define MaxHits1 10
#endif
#ifndef MaxHits2
#define MaxHits2 60
#endif
#ifndef MaxHits3
#define MaxHits3 80
#endif

const int TOFSTOREMAX = 10;

namespace
{
using namespace root;
}

struct Event
{
  std::map<TString, std::vector<TParticle>> hits;

  double x0In;
  double y0In;
  double z0In;
  double u0In;
  double v0In;
  double uDeg;
  double vDeg;
  double phi0;
  double theta0;
  double p0;
  double t0;
  int Id;

  G4int nP, nN, nL, nSm, nSz, nSp, nXm, nXz, nXsm, nXsz;
  G4int nPim, nPiz, nPip, nKm, nKp;

  double DCXObs[NumDC];
  double DCYObs[NumDC];
  double DCX[NumDC];
  double DCY[NumDC];
  double DCt[NumDC];
  double DCp[NumDC];
  Int_t DCNh[NumDC];
  Double_t DCgPosx[NumDC][MaxHits];
  Double_t DCgPosy[NumDC][MaxHits];
  Double_t DCgPosz[NumDC][MaxHits];
  Double_t DCde[NumDC][MaxHits];

  int DCNhits;
  int DC1Hit;
  int DC2Hit;
  int DC3Hit;
  int DC4Hit;
  int DC5Hit;
  double SlituDeg[11];
  double SlitvDeg[11];
  double SlitX[11][MaxHits];
  double SlitY[11];
  double Slitt[11];
  double SlitMom[11];
  double Slitp[11];
  int SlitNh[11];
  int SlitNP[11];
  int SlitNK[11];
  int SlitNPi[11];
  int SlitF[11];

  // ---- TOF detectors ----
  //int tofco[18];
  int TOFNhits;
  double toftime[18];
  double toftime_reso[18];
  double tofdE[18];
  double tofn[18];

  // ---- Water Cherenkov detector ----
  int WCNhits;
  double wctime[12];
  double wcdE[12];
  double wcn[12];
  double wcnpe[12];
  double wctime1[6];
  double wcdE1[6];
  double wcn1[6];
  double wcnpe1[6];
  double wctime2[6];
  double wcdE2[6];
  double wcn2[6];
  double wcnpe2[6];
  //int wclayer;
  //int wcseg;

  // ----- Trigger -----
  G4bool TOFTrig;
  G4bool VDTrig;
  G4bool WCTrig;
  G4bool Q1Trig;
  G4bool Q2Trig;

  //double toftime[TOFSTOREMAX];
  //double tofdE[TOFSTOREMAX];
//   double TOFAll;
//   double TOFt[NumTOFSeg];
//   double TOFtObs[NumTOFSeg];


//   int TOFHit;
//   double ACX;
//   double ACY;
//   double ACt;
//   double ACp;
//  int ACNhits;
  ///int ACHit;
//   double WCt[NumWCSeg];
//   double WCp[NumWCSeg];
  //int WCNhits;
  //int WCHit;
};

class S2SAnaManager
{
public:
  static G4String       ClassName();
  static S2SAnaManager& GetInstance();
  virtual ~S2SAnaManager();

private:
  S2SAnaManager();
  S2SAnaManager(const S2SAnaManager&);
  S2SAnaManager& operator =(const S2SAnaManager&);

public:
  void BeginOfRun( const G4Run *aRun );
  void EndOfRun( const G4Run *aRun );
  void BeginOfPrimaryAction();
  void MakeBranch(const G4String& sd_name);
  void MakeHistogram(const G4String& sd_name);
  void SetNhits(const G4String& sd_name, G4int nhits);
  void SetHitData(const VHitInfo* hit);
  void SetPrimaryData(G4double x0, G4double y0, G4double z0,
		      G4double u0, G4double v0, G4double phi,
		      G4double theta, G4double p0, G4double t0,
		      G4int ParIdNb);
  void SetProcessData(G4int nP, G4int nN, G4int nL,
		      G4int nSm, G4int nSz, G4int nSp,
		      G4int nXm, G4int nXz, G4int nXsm,
		      G4int nXsz,G4int nPim,G4int nPiz,
		      G4int nPip,G4int nKm,G4int nKp);
  void BeginOfEvent( const G4Event *anEvent );
  void EndOfEvent( const G4Event *anEvent );
  void SetFileName(const G4String& name){ m_file_name = name; }
  //void DefineHistograms();
  void DefineTree();
  G4bool GetTriggerStatus() const { return fTriggered; }
  void SetActive() { fActive_=true; }
  void SetInActive() { fActive_=false; }

private:
  G4String m_file_name;
  G4bool fActive_;
  G4bool fTriggered;
  TRandom3* nperand;

  G4int trigNum;
  TFile *m_file;
  TTree *m_tree;

private:
  void PrintHitsInformation( const G4Event *anEvent,
			    std::ostream &ost ) const;
public:
  void InitializeEvent(void);
};

//_____________________________________________________________________________
inline G4String
S2SAnaManager::ClassName()
{
  static G4String s_name("S2SAnaManager");
  return s_name;
}


//_____________________________________________________________________________
inline S2SAnaManager&
S2SAnaManager::GetInstance()
{
  static S2SAnaManager s_instance;
  return s_instance;
}

#endif
