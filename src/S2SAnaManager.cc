// -*- C++ -*-

#include "S2SAnaManager.hh"

#include <bitset>

#include "G4Run.hh"
#include "G4Event.hh"
#include "G4Step.hh"
#include "G4ios.hh"
#include "G4SDManager.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
//#include "G4Poisson.h"

#include "Randomize.hh"

#include <TMath.h>
#include <TRandom3.h>
#include <TSystem.h>

#include "RootHelper.hh"

#include <iomanip>
#include <cmath>
#include <sstream>
#include <vector>

#include "ConfMan.hh"
#include "DCGeomMan.hh"
#include "DetectorID.hh"
#include "FuncName.hh"
#include "DCHit.hh"
#include "TOFHit.hh"
#include "HTOFHit.hh"
#include "TPCHit.hh"
#include "ACHit.hh"
#include "WCHit.hh"
#include "VPHit.hh"
#include "HistMan.hh"

#include "RCHit.hh"
#include "PDHit.hh"


namespace
{
using CLHEP::mm;
using CLHEP::MeV;
using CLHEP::ns;
const auto& confMan = ConfMan::GetInstance();
const auto& histMan = HistMan::GetInstance();
const auto qnan = TMath::QuietNaN();
Event event;
std::map<TString, TH1*> hmap;
std::vector<G4int> n_acc(kTriggerFlagSize, 0);
}

//_____________________________________________________________________________
S2SAnaManager::S2SAnaManager()
  : m_file_name("tmp.root"),
    fActive_(true),
    fTriggered(false),
    m_file(),
    m_tree()
{
}

//_____________________________________________________________________________
S2SAnaManager::~S2SAnaManager()
{
}

//_____________________________________________________________________________
void
S2SAnaManager::BeginOfRun( const G4Run* /* aRun */)
{
  fActive_=true;
  m_file = new TFile(m_file_name, "recreate");
  //static auto obj = new TNamed("conf", confMan.ConfPath()+confMan.ConfBuf());
  static auto obj = new TNamed("conf", confMan.ConfBuf()); // [seong]
  obj->Write();
  static auto git = new TNamed
    ("git", ("\n"+gSystem->GetFromPipe("git log -1")).Data());
  git->Write();
  m_tree = new TTree("g4s2s", "S-2S simulation");
  event.hits.clear();
  event.evnum = -1;
  event.trig.assign(kTriggerFlagSize, false);
  event.rctrig.assign(kRCTriggerFlagSize, false); //[seong]
  DefineTree();
  static const G4int experiment = confMan.Get<G4int>("Experiment");
  for(const auto& sd_name : std::vector<G4String>{
      "PRM", "SDC1","SDC2","SDC3","SDC4","SDC5", "TOF", "AC1", "WC", "VP"}
        // S2SDetectorConstruction::GetSDList()
    ){
    G4cout << "   make branch : " << sd_name << G4endl;
    MakeBranch(sd_name);
    MakeHistogram(sd_name);
  }
  if(experiment == 63){
    for(const auto& sd_name : std::vector<G4String>{
	"RC1", "RC2", "PD1Y", "PD1Z", "PD2Y", "PD2Z"}
	  // S2SDetectorConstruction::GetSDList()
      ){
      G4cout << "   make branch : " << sd_name << G4endl;
      MakeBranch(sd_name);
      MakeHistogram(sd_name);
    }
  }
  if(experiment == 90){
    for(const auto& sd_name : std::vector<G4String>{"TPC", "HTOF", "SAC"})
    {
      G4cout << "   make branch : " << sd_name << G4endl;
      MakeBranch(sd_name);
      MakeHistogram(sd_name);
    }
  }
  for(auto& h: hmap){
    h.second->Reset();
  }
  n_acc.clear();
  n_acc.resize(kTriggerFlagSize);
}

//_____________________________________________________________________________
void
S2SAnaManager::EndOfRun(const G4Run* /* aRun */)
{
  m_file->cd();
  if(confMan.Get<G4bool>("TREE"))
    m_tree->Write();
  for(auto& h: hmap){
    h.second->Write();
  }
  m_file->Close();
}

//_____________________________________________________________________________
void
S2SAnaManager::BeginOfPrimaryAction()
{
  event.x0In = qnan;
  event.y0In = qnan;
  event.u0In = qnan;
  event.v0In = qnan;
  event.uDeg = qnan;
  event.vDeg = qnan;
  event.p0   = qnan;
  event.pB   = qnan;
  event.t0   = qnan;
  event.phi0 = qnan;
  event.theta0 = qnan;
  // event.Id = qnan;
  //  G4cout<<"BeginOfPrimaryAction"<<G4endl;

#if 0
  // E63 for weak pi-
  event.x1 = qnan;
  event.y1 = qnan;
  event.z1 = qnan;
  event.px1 = qnan;
  event.py1 = qnan;
  event.pz1 = qnan;
  event.p1  = qnan;
  event.t1  = qnan;
#endif
}

//_____________________________________________________________________________
void
S2SAnaManager::MakeBranch(const G4String& sd_name)
{
  static const Int_t bufsize = 32000;
  m_tree->Branch(sd_name.data(),
                 "std::vector<TParticle>",
                 &event.hits[sd_name], bufsize, -1);
}

//_____________________________________________________________________________
void
S2SAnaManager::MakeHistogram(const G4String& sd_name)
{
  if(sd_name == "PRM"){
    const auto& params = histMan.Get("PTheta");
    TString key = sd_name + "PThetaGen";
    TString title = sd_name + " P%Theta (Generate); [deg]; [GeV/c]";
    hmap[key] = new TH2D(key, title,
                         params.at(0), params.at(1), params.at(2),
                         params.at(3), params.at(4), params.at(5));
    key = sd_name + "PThetaAcc";
    title = sd_name + " P%Theta (Accept); [deg]; [GeV/c]";
    hmap[key] = new TH2D(key, title,
                         params.at(0), params.at(1), params.at(2),
                         params.at(3), params.at(4), params.at(5));
    for(G4int i=0, n=TriggerFlag.size(); i<n; ++i){
      key = sd_name + "PThetaAcc" + TriggerFlag.at(i);
      title = sd_name + " P%Theta (Accept at " + TriggerFlag.at(i) + ");"
        + " [deg]; [GeV/c]";
      hmap[key] = new TH2D(key, title,
                           params.at(0), params.at(1), params.at(2),
                           params.at(3), params.at(4), params.at(5));
    }
  }else{
    for(const auto& suffix: std::vector<G4String>
          { "Nhits", "HitPat", "X", "Y", "Z", "U", "V",
            "Y%X", "V%U", "U%X", "V%Y" }){
      TString key = sd_name + suffix;
      const auto& params = histMan.Get(key);
      TString title = sd_name + " " + suffix;
      key.ReplaceAll("%", "");
      if(G4StrUtil::contains(suffix, "%")){
        hmap[key] = new TH2D(key, title,
                             params.at(0), params.at(1), params.at(2),
                             params.at(3), params.at(4), params.at(5));
      }else{
        hmap[key] = new TH1D(key, title,
                             params.at(0), params.at(1), params.at(2));
      }
    }
  }
}

//_____________________________________________________________________________
void
S2SAnaManager::SetPrimaryData(double x0, double y0, double z0,
                              double u0, double v0, double phi, double theta,
                              double p0,double pB, int /* ParIdNb */)
{
  event.x0In = x0; // generated position (x)
  event.y0In = y0; // generated position (y)
  event.z0In = z0; // generated position (z)
  //event.u0In = u0; // x' in rad
  event.u0In = -1.0 * u0; // x' in rad
  event.v0In = v0; // y' in rad
  //event.uDeg = atan(u0)*TMath::RadToDeg(); // x' in deg
  //event.vDeg = atan(v0)*TMath::RadToDeg(); // y' in deg
  //event.uDeg = -1.0 * u0 * TMath::RadToDeg(); // x' in deg
  //event.vDeg = v0*TMath::RadToDeg(); // y' in deg
  event.uDeg = u0;
  event.vDeg = v0;
  event.phi0 = phi; // Phi in rad
  event.theta0 = theta; // Theta in rad
  event.p0   = p0; // Momentum
  event.pB   = pB; // Kinetic energy
  //  event.Id = ParIdNb;
  //  G4cout<<"setPrimaryData"<<G4endl;
}

// E63 for weak pi
//_____________________________________________________________________________
void
S2SAnaManager::SetSecondaryData(double x1, double y1, double z1,
				double px1, double py1, double pz1,
				double p1,double t1)
{
  event.x1 = x1; // generated position (x)
  event.y1 = y1; // generated position (y)
  event.z1 = z1; // generated position (z)
  event.px1 = px1;
  event.py1 = py1;
  event.pz1 = pz1;
  event.p1   = p1; // Momentum
  event.t1   = t1; // Kinetic energy
  //  event.Id = ParIdNb;
  //  G4cout<<"setPrimaryData"<<G4endl;
}

//_____________________________________________________________________________
void S2SAnaManager::SetProcessData(G4int nP, G4int nN, G4int nL,
				 G4int nSm, G4int nSz, G4int nSp,
				 G4int nXm, G4int nXz, G4int nXsm,
				 G4int nXsz,G4int nPim,G4int nPiz,
				 G4int nPip,G4int nKm,G4int nKp)
{
  event.nP = nP;
  event.nN = nN;
  event.nL = nL;
  event.nSm = nSm;
  event.nSz = nSz;
  event.nSp = nSp;
  event.nXm = nXm;
  event.nXz = nXz;
  event.nXsm = nXsm;
  event.nXsz = nXsz;
  event.nPim = nPim;
  event.nPiz = nPiz;
  event.nPip = nPip;
  event.nKm = nKm;
  event.nKp = nKp;
}

void S2SAnaManager::BeginOfEvent(const G4Event *anEvent)
{
  event.evnum = anEvent->GetEventID();
}

void S2SAnaManager::EndOfEvent(const G4Event *anEvent)
{
  auto HCE = anEvent->GetHCofThisEvent();
  auto SDMan = G4SDManager::GetSDMpointer();
  std::bitset<kTriggerFlagSize> trigger_flag;
  std::bitset<kRCTriggerFlagSize> rc_trigger_flag;
  static const G4int experiment = confMan.Get<G4int>("Experiment");
  G4String particle_name = "kaon+";
  if(experiment == 63) particle_name = "pi-"; //for E63
  //G4String particle_name = "kaon-"; //for E63

  {
    //G4String name = "SDC"+std::to_string(k);
    G4String name = "SDC1";
    static const auto id = SDMan->GetCollectionID(name);
    if(id >= 0){
      //auto HC = dynamic_cast<SDCHitsCollection*>(HCE->GetHC(id));
      auto HC = dynamic_cast<DCHitsCollection*>(HCE->GetHC(id));
      for(G4int i=0, n=HC->entries(); i<n; ++i){
        SetHitData((*HC)[i]);
      }
      SetNhits(name, HC->entries());
    }
  }
  {
    G4String name = "SDC2";
    static const auto id = SDMan->GetCollectionID(name);
    if(id >= 0){
      //auto HC = dynamic_cast<SDCHitsCollection*>(HCE->GetHC(id));
      auto HC = dynamic_cast<DCHitsCollection*>(HCE->GetHC(id));
      for(G4int i=0, n=HC->entries(); i<n; ++i){
        SetHitData((*HC)[i]);
      }
      SetNhits(name, HC->entries());
    }
  }
  {
    G4String name = "SDC3";
    static const auto id = SDMan->GetCollectionID(name);
    if(id >= 0){
      //auto HC = dynamic_cast<SDCHitsCollection*>(HCE->GetHC(id));
      auto HC = dynamic_cast<DCHitsCollection*>(HCE->GetHC(id));
      for(G4int i=0, n=HC->entries(); i<n; ++i){
        SetHitData((*HC)[i]);
      }
      SetNhits(name, HC->entries());
    }
  }
  {
    G4String name = "SDC4";
    static const auto id = SDMan->GetCollectionID(name);
    if(id >= 0){
      //auto HC = dynamic_cast<SDCHitsCollection*>(HCE->GetHC(id));
      auto HC = dynamic_cast<DCHitsCollection*>(HCE->GetHC(id));
      for(G4int i=0, n=HC->entries(); i<n; ++i){
        SetHitData((*HC)[i]);
      }
      SetNhits(name, HC->entries());
    }
  }
  {
    G4String name = "SDC5";
    static const auto id = SDMan->GetCollectionID(name);
    if(id >= 0){
      //auto HC = dynamic_cast<SDCHitsCollection*>(HCE->GetHC(id));
      auto HC = dynamic_cast<DCHitsCollection*>(HCE->GetHC(id));
      for(G4int i=0, n=HC->entries(); i<n; ++i){
        SetHitData((*HC)[i]);
      }
      SetNhits(name, HC->entries());
    }
  }
  /*static const auto id = SDMan->GetCollectionID("SDC");
  if(id >= 0){
    //auto HC = dynamic_cast<SDCHitsCollection*>(HCE->GetHC(id));
    auto HC = dynamic_cast<DCHitsCollection*>(HCE->GetHC(id));
    for(G4int i=0, n=HC->entries(); i<n; ++i){
      SetHitData((*HC)[i]);
    }
    SetNhits("SDC", HC->entries());
   }*/
  {
    static const auto id = SDMan->GetCollectionID("TOF");
    if(id >= 0){
      auto HC = dynamic_cast<TOFHitsCollection*>(HCE->GetHC(id));
      for(G4int i=0, n=HC->entries(); i<n; ++i){
        auto hit = (*HC)[i];
        if(hit->Is(particle_name) && hit->IsPrimary()) trigger_flag[kTOF] = true;
        SetHitData(hit);
      }
      SetNhits("TOF", HC->entries());
    }
  }
  {
    static const auto id = SDMan->GetCollectionID("AC1");
    if(id >= 0){
      auto HC = dynamic_cast<ACHitsCollection*>(HCE->GetHC(id));
      for(G4int i=0, n=HC->entries(); i<n; ++i){
        auto hit = (*HC)[i];
        if(hit->Is(particle_name) && hit->IsPrimary()) trigger_flag[kAC1] = true;
        SetHitData(hit);
      }
      SetNhits("AC1", HC->entries());
    }
  }
  {
    static const auto id = SDMan->GetCollectionID("WC");
    if(id >= 0){
      auto HC = dynamic_cast<WCHitsCollection*>(HCE->GetHC(id));
      for(G4int i=0, n=HC->entries(); i<n; ++i){
        auto hit = (*HC)[i];
        if(hit->Is(particle_name) && hit->IsPrimary()) trigger_flag[kWC] = true;
        SetHitData(hit);
      }
      SetNhits("WC", HC->entries());
    }
  }
  {
    static const auto id = SDMan->GetCollectionID("VP");
    if(id >= 0){
      auto HC = dynamic_cast<VPHitsCollection*>(HCE->GetHC(id));
      for(G4int i=0, n=HC->entries(); i<n; ++i){
        auto hit = (*HC)[i];
        if(hit->Is(particle_name) && hit->IsPrimary()){
          trigger_flag[kVP1-1+hit->GetCopyNumber()] = true;
        }
        SetHitData(hit);
      }
      SetNhits("VP", HC->entries());
    }
  }
  if(experiment == 63)
  {
    {
      G4String name = "RC1";
      static const auto id = SDMan->GetCollectionID(name);
      if(id >= 0){
	auto HC = dynamic_cast<RCHitsCollection*>(HCE->GetHC(id));
	for(G4int i=0, n=HC->entries(); i<n; ++i){
	  auto hit = (*HC)[i];
	  if(hit->Is(particle_name) && hit->IsWeakPi()) rc_trigger_flag[kRC1] = true;
	  SetHitData((*HC)[i]);
	}
	SetNhits(name, HC->entries());
      }
    }
    {
      G4String name = "RC2";
      static const auto id = SDMan->GetCollectionID(name);
      if(id >= 0){
	auto HC = dynamic_cast<RCHitsCollection*>(HCE->GetHC(id));
	for(G4int i=0, n=HC->entries(); i<n; ++i){
	  auto hit = (*HC)[i];
	  if(hit->Is(particle_name) && hit->IsWeakPi()) rc_trigger_flag[kRC2] = true;
	  SetHitData((*HC)[i]);
	}
	SetNhits(name, HC->entries());
      }
    }
    {
      G4String name = "PD1Y";
      static const auto id = SDMan->GetCollectionID(name);
      if(id >= 0){
	auto HC = dynamic_cast<PDHitsCollection*>(HCE->GetHC(id));
	for(G4int i=0, n=HC->entries(); i<n; ++i){
	  auto hit = (*HC)[i];
	  if(hit->Is(particle_name) && hit->IsWeakPi()) rc_trigger_flag[kPD1Y] = true;
	  SetHitData((*HC)[i]);
	}
	SetNhits(name, HC->entries());
      }
    }
    {
      G4String name = "PD1Z";
      static const auto id = SDMan->GetCollectionID(name);
      if(id >= 0){
	auto HC = dynamic_cast<PDHitsCollection*>(HCE->GetHC(id));
	for(G4int i=0, n=HC->entries(); i<n; ++i){
	  auto hit = (*HC)[i];
	  if(hit->Is(particle_name) && hit->IsWeakPi()) rc_trigger_flag[kPD1Z] = true;
	  SetHitData((*HC)[i]);
	}
	SetNhits(name, HC->entries());
      }
    }
    {
      G4String name = "PD2Y";
      static const auto id = SDMan->GetCollectionID(name);
      if(id >= 0){
	auto HC = dynamic_cast<PDHitsCollection*>(HCE->GetHC(id));
	for(G4int i=0, n=HC->entries(); i<n; ++i){
	  auto hit = (*HC)[i];
	  if(hit->Is(particle_name) && hit->IsWeakPi()) rc_trigger_flag[kPD2Y] = true;
	  SetHitData((*HC)[i]);
	}
	SetNhits(name, HC->entries());
      }
    }
    {
      G4String name = "PD2Z";
      static const auto id = SDMan->GetCollectionID(name);
      if(id >= 0){
	auto HC = dynamic_cast<PDHitsCollection*>(HCE->GetHC(id));
	for(G4int i=0, n=HC->entries(); i<n; ++i){
	  auto hit = (*HC)[i];
	  if(hit->Is(particle_name) && hit->IsWeakPi()) rc_trigger_flag[kPD2Z] = true;
	  SetHitData((*HC)[i]);
	}
	SetNhits(name, HC->entries());
      }
    }
  }

  if (experiment == 90)
  {
    {
      static const auto id = SDMan->GetCollectionID("HTOF");
      if(id >= 0){
        auto HC = dynamic_cast<HTOFHitsCollection*>(HCE->GetHC(id));
        for(G4int i=0, n=HC->entries(); i<n; ++i){
          auto hit = (*HC)[i];
          SetHitData(hit);
        }
        SetNhits("HTOF", HC->entries());
      }
    }
    {
      static const auto id = SDMan->GetCollectionID("TPC");
      if(id >= 0){
        auto HC = dynamic_cast<TPCHitsCollection*>(HCE->GetHC(id));
        for(G4int i=0, n=HC->entries(); i<n; ++i){
          auto hit = (*HC)[i];
          SetHitData(hit);
        }
        SetNhits("TPC", HC->entries());
      }
    }
    {
      static const auto id = SDMan->GetCollectionID("SAC");
      if(id >= 0){
        auto HC = dynamic_cast<ACHitsCollection*>(HCE->GetHC(id));
        for(G4int i=0, n=HC->entries(); i<n; ++i){
          auto hit = (*HC)[i];
          // if(hit->Is("pi-") && hit->IsPrimary()) trigger_flag[kSAC] = true;
          SetHitData(hit);
        }
        SetNhits("SAC", HC->entries());
      }
    }
  }

  {
    // G4cout << "Acc eff." << G4endl;
    auto particle = event.hits.at("PRM").at(0);
    for(G4int i=0, n=TriggerFlag.size(); i<n; ++i){
      if(trigger_flag[i]){
        n_acc[i]++;
        hmap.at("PRMPThetaAcc"+TriggerFlag.at(i))->
          Fill(particle.Theta()/CLHEP::deg, particle.P()/CLHEP::GeV);
      }
      // G4cout << "   " << TriggerFlag.at(i) << "\t"
      //        << (G4double)n_acc[i]/(event.evnum+1) << G4endl;
    }

    if(true
       && trigger_flag[kVP1]
       && trigger_flag[kVP2]
       && trigger_flag[kVP3]
       && trigger_flag[kVP4]
       && trigger_flag[kVP5]
       && trigger_flag[kVP6]
       && trigger_flag[kVP7]
       && trigger_flag[kVP8]
       && trigger_flag[kVP9]
       && trigger_flag[kVP10]
       && trigger_flag[kTOF]
       && trigger_flag[kWC]
       ){
      hmap.at("PRMPThetaAcc")->
        Fill(particle.Theta()/CLHEP::deg, particle.P()/CLHEP::GeV);
    }
  }

  {
    for(G4int i=0; i<kTriggerFlagSize; ++i){
      event.trig[i] = trigger_flag[i];
    }
    for(G4int i=0; i<kRCTriggerFlagSize; ++i){ //[seong]
      event.rctrig[i] = rc_trigger_flag[i];
    }
  }

  if(confMan.Get<G4bool>("TREE"))
    m_tree->Fill();

  InitializeEvent();
}

//_____________________________________________________________________________
void
S2SAnaManager::SetNhits(const G4String& sd_name, G4int nhits)
{
  hmap[sd_name + "Nhits"]->Fill(nhits);
}

//_____________________________________________________________________________
void
S2SAnaManager::SetHitData(const VHitInfo* hit)
{
  if(hit && hit->GetParticle()){
    const auto& name = hit->GetDetectorName();
    const auto& p = hit->GetParticle();
    event.hits.at(name).push_back(*p);
    hmap[name + "HitPat"]->Fill(p->GetMother(1));
    hmap[name + "X"]->Fill(p->Vx());
    hmap[name + "Y"]->Fill(p->Vy());
    hmap[name + "Z"]->Fill(p->Vz());
    hmap[name + "U"]->Fill(p->Px()/p->Pz());
    hmap[name + "V"]->Fill(p->Py()/p->Pz());
    hmap[name + "YX"]->Fill(p->Vx(), p->Vy());
    hmap[name + "VU"]->Fill(p->Px()/p->Pz(), p->Py()/p->Pz());
    hmap[name + "UX"]->Fill(p->Vx(), p->Px()/p->Pz());
    hmap[name + "VY"]->Fill(p->Vy(), p->Py()/p->Pz());
  }
}

void S2SAnaManager::InitializeEvent()
{
  for(auto& pair: event.hits){
    pair.second.clear();
  }
}

void S2SAnaManager::DefineTree()
{
  m_tree->Branch("evnum", &event.evnum, "evnum/I");
  m_tree->Branch("trig", &event.trig);
  m_tree->Branch("x0",&event.x0In, "x0/D");
  m_tree->Branch("y0",&event.y0In, "y0/D");
  m_tree->Branch("z0",&event.z0In, "z0/D");
  //   m_tree->Branch("u0In",&event.u0In, "u0In/D");
  //   m_tree->Branch("v0In",&event.v0In, "v0In/D");
  m_tree->Branch("xp0",&event.uDeg, "xp0/D");
  m_tree->Branch("yp0",&event.vDeg, "yp0/D");
  m_tree->Branch("phi0",&event.phi0, "phi0/D");
  m_tree->Branch("theta0",&event.theta0, "theta0/D");
  m_tree->Branch("p0",&event.p0,   "p0/D");
  m_tree->Branch("pB",&event.pB,   "pB/D");

  // E63 for weak pi
#if 0
  m_tree->Branch("rctrig", &event.rctrig);
  m_tree->Branch("x1",  &event.x1,  "x1/D"); // [mm]
  m_tree->Branch("y1",  &event.y1,  "y1/D");
  m_tree->Branch("z1",  &event.z1,  "z1/D");
  m_tree->Branch("px1", &event.px1, "px1/D"); // [MeV/c]
  m_tree->Branch("py1", &event.py1, "py1/D");
  m_tree->Branch("pz1", &event.pz1, "pz1/D");
  m_tree->Branch("p1",  &event.p1,  "p1/D");  // [MeV/c]
  m_tree->Branch("t1",  &event.t1,  "t1/D");  // [MeV]
#endif

  return;
  //  m_tree->Branch("t0",&event.t0,   "t0/D");
  // m_tree->Branch("Id",&event.Id, "Id/I");
  //   m_tree->Branch("nP",&event.nP,"nP/I");
  //   m_tree->Branch("nN",&event.nN,"nN/I");
  //   m_tree->Branch("nL",&event.nL,"nL/I");
  //   m_tree->Branch("nSm",&event.nSm,"nSm/I");
  //   m_tree->Branch("nSz",&event.nSz,"nSz/I");
  //   m_tree->Branch("nSp",&event.nSp,"nSp/I");
  //   m_tree->Branch("nXm",&event.nXm,"nXm/I");
  //   m_tree->Branch("nXz",&event.nXz,"nXz/I");
  //   m_tree->Branch("nXsm",&event.nXsm,"nXsm/I");
  //   m_tree->Branch("nXsz",&event.nXsz,"nXsz/I");
  //   m_tree->Branch("nPim",&event.nPim,"nPim/I");
  //   m_tree->Branch("nPiz",&event.nPiz,"nPiz/I");
  //   m_tree->Branch("nPip",&event.nPip,"nPip/I");
  //   m_tree->Branch("nKm",&event.nKm,"nKm/I");
  //   m_tree->Branch("nKp",&event.nKp,"nKp/I");

  m_tree->Branch("vdxp", event.SlituDeg, "vdxp[11]/D");
  m_tree->Branch("vdyp", event.SlitvDeg, "vdyp[11]/D");
  /*
    m_tree->Branch("SlitX", event.SlitX, "SlitX[8]/D");
    m_tree->Branch("SlitY", event.SlitY, "SlitY[8]/D");
    m_tree->Branch("Slitt", event.Slitt, "Slitt[8]/D");
    m_tree->Branch("SlitMom", event.SlitMom, "SlitMom[8]/D");
    m_tree->Branch("Slitp", event.Slitp, "Slitp[8]/D");
    m_tree->Branch("SlitNh", event.SlitNh, "SlitNh[8]/I");
    m_tree->Branch("SlitNP", event.SlitNP, "SlitNP[8]/I");
    m_tree->Branch("SlitNK", event.SlitNK, "SlitNK[8]/I");
    m_tree->Branch("SlitNPi", event.SlitNPi, "SlitNPi[8]/I");
    m_tree->Branch("SlitF", event.SlitF, "SlitF[8]/I");
  */

  // Order -->
  // vd8 | TOF | vd9 | AC | vd10 | WC | vd11 (vd[10]) | SDC1
  m_tree->Branch("vdx",   event.SlitX, Form("vdx[11][%d]/D", MaxHits));
  m_tree->Branch("vdy",   event.SlitY, "vdy[11]/D");
  m_tree->Branch("vdtime",event.Slitt, "vdtime[11]/D");
  m_tree->Branch("vdmom", event.SlitMom, "vdmom[11]/D");
  m_tree->Branch("vdpath",event.Slitp, "vdpath[11]/D");
  m_tree->Branch("vdNh",  event.SlitNh, "vdNh[11]/I");
  m_tree->Branch("vdNP",  event.SlitNP, "vdNP[11]/I");
  m_tree->Branch("vdNK",  event.SlitNK, "vdNK[11]/I");
  m_tree->Branch("vdNPi", event.SlitNPi, "vdNPi[11]/I");
  m_tree->Branch("vdPID", event.SlitF, "vdPID[11]/I");
  // m_tree->Branch("TOFAll", &event.TOFAll, "TOFAll/D");
  // m_tree->Branch("TOFt",event.TOFt,   "TOFt[17]/D");
  // m_tree->Branch("TOFtObs",event.TOFtObs,   "TOFtObs[17]/D");
  //m_tree->Branch("tofco",  &event.tofco,   "tofco[17]/I");

  m_tree->Branch("tofnhits",&event.TOFNhits, "tofnhits/I");
  m_tree->Branch("toftime", &event.toftime,  "toftime[18]/D");
  m_tree->Branch("tofdE",   &event.tofdE,    "tofdE[18]/D");

  m_tree->Branch("wcnhits",&event.WCNhits,  "wcnhits/I");
  //m_tree->Branch("wctime", &event.wctime,   "wctime[12]/D");
  //m_tree->Branch("wcdE",   &event.wcdE,     "wcdE[12]/D");
  //m_tree->Branch("wcnpe",  &event.wcnpe,    "wcnpe[12]/D");
  m_tree->Branch("wctime1", &event.wctime1,   "wctime1[6]/D");
  m_tree->Branch("wcdE1",   &event.wcdE1,     "wcdE1[6]/D");
  m_tree->Branch("wcnpe1",  &event.wcnpe1,    "wcnpe1[6]/D");
  m_tree->Branch("wctime2", &event.wctime2,   "wctime2[6]/D");
  m_tree->Branch("wcdE2",   &event.wcdE2,     "wcdE2[6]/D");
  m_tree->Branch("wcnpe2",  &event.wcnpe2,    "wcnpe2[6]/D");

  m_tree->Branch("Q1Trig", &event.Q1Trig,  "Q1Trig/B");
  m_tree->Branch("Q2Trig", &event.Q2Trig,  "Q2Trig/B");
  m_tree->Branch("VDTrig", &event.VDTrig,  "VDTrig/B");
  // m_tree->Branch("TOFTrig",&event.TOFTrig, "TOFTrig/B");

  m_tree->Branch("WCTrig", &event.WCTrig,  "WCTrig/B");
  // m_tree->Branch("TOFHit", &event.TOFHit, "TOFHit/I");
  // m_tree->Branch("ACX", &event.ACX, "ACX/D");
  // m_tree->Branch("ACY", &event.ACY, "ACY/D");
  // m_tree->Branch("ACt", &event.ACt, "ACt/D");
  // m_tree->Branch("ACp", &event.ACp, "ACp/D");
  // m_tree->Branch("ACNhits", &event.ACNhits, "ACNhits/I");
  // m_tree->Branch("ACHit", &event.ACHit, "ACHit/I");
  // m_tree->Branch("WCt", event.WCt,  "WCt[16]/D");
  // m_tree->Branch("WCp", event.WCp,  "WCp[16]/D");
  // m_tree->Branch("WCNhits", &event.WCNhits, "WCNhits/I");
  // m_tree->Branch("WCHit", &event.WCHit, "WCHit/I");
  //G4cout<<"tree and Branch is defined"<<G4endl;

}

//_____________________________________________________________________________
void
S2SAnaManager::SetPrimaryParticle(G4int id, G4int pdg,
                                  const G4LorentzVector& p,
                                  const G4LorentzVector& v,
                                  G4bool is_virtual_beam)
{
  G4int id1 = is_virtual_beam ? -1 : 1;
  G4int id2 = id;
  for (const auto& ptcl: event.hits.at("PRM")) {
    if (ptcl.GetMother(0) == id1 && ptcl.GetMother(1) == id2) {
      G4cerr << FUNC_NAME << " id1=" << id1 << ", id2=" << id2
             << " is already set" << G4endl;
    }
  }
  TParticle particle(pdg,
                     0, // fStatus
                     id1, // fMother[0]
                     id2, // fMother[1]
                     0, // fDaughter[0]
                     0, // fDaughter[1]
                     TLorentzVector(p.px(), p.py(), p.pz(), p.e()),
                     TLorentzVector(v.x(), v.y(), v.z(), v.t()));
  event.hits.at("PRM").push_back(particle);

  if(id == 0){
    hmap.at("PRMPThetaGen")->Fill(p.theta()/CLHEP::degree,
                                  p.v().mag()/CLHEP::GeV);
  }
}

//_____________________________________________________________________________
void S2SAnaManager::PrintHitsInformation(const G4Event* /* anEvent */,
                                         std::ostream &ost) const
{
  //   G4cout<<"PrintHits is called"<<G4endl;

  //  int GeomFlag = confMan->GeomFlag();

  // G4HCofThisEvent *HCE = anEvent->GetHCofThisEvent();
  // G4SDManager *SDMan = G4SDManager::GetSDMpointer();

  std::ios::fmtflags oldFlags = ost.flags();
  std::size_t preSiz = ost.precision();
  ost.setf( std::ios::fixed );

  ost.precision(5);
  ost << event.x0In << std::setw(10)
      << event.y0In << std::setw(10)
      <<" "<< -4320.62-event.z0In << std::setw(10);
  ost.precision(5);

  ost << event.p0     << std::setw(10)
      << event.theta0 << std::setw(10)
      << event.phi0   << std::setw(10);
  ost.precision(7);

  ost << 1.8 <<std::setw(10);
  ost << 0   <<std::setw(10);
  ost << 0   <<std::setw(10);

  ost.precision(5);

  //G4cout<<"Printed"<<G4endl;
  //G4cout<<"osf="<<ost<<G4endl;
  //G4cout<<"x0="<<event.x0In<<" y0="<<event.y0In<<G4endl;

  // const DCGeomMan & geomMan=DCGeomMan::GetInstance();

  // G4double nhTOF=0;
  // TOFHitsCollection    *TOFHC;
  // G4int colIdTOF = SDMan->GetCollectionID( "TOFSD"/*"DCCollection"*/ );
  // TOFHC =  dynamic_cast<TOFHitsCollection *>(  HCE->GetHC( colIdTOF ) );
  // if( TOFHC )     nhTOF     = TOFHC ->entries();

  // for( int i=0; i<nhTOF; ++i ){
  //   TOFHit *aHit = (*TOFHC)[i];
  //   double time = aHit->GetTime();
  //   double de = aHit->GetEdep();
  //   if(i==nhTOF-1){
  //     ost << std::setw(12) << 61
  //         << std::setw(12) << time;
  //     de=1.0;
  //     ost << std::setw(12) << 62 <<" "
  //         << std::setw(12) << de;
  //   }
  // }

  ost << std::setw(5) << -1 << std::endl;

  ost.flags( oldFlags );
  ost.precision( preSiz );

  ost << std::endl;

}
