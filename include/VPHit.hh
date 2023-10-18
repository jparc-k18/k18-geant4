// -*- C++ -*-

#ifndef VP_HIT_HH
#define VP_HIT_HH

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4ThreeVector.hh"
#include "G4Allocator.hh"

//_____________________________________________________________________________
class VPHit : public G4VHit
{
public:
  VPHit();
  ~VPHit(){}

private:
  VPHit(const VPHit&);
  VPHit& operator =(const VPHit&);

public:
  int operator ==(const VPHit&) const { return 0; }

  inline void * operator new(size_t size);
  inline void operator delete(void* aHit);

//   void Draw() const;
//   void Print() const;

private:
  G4int pass_;
  G4int layerID_;
  G4int wireID_;
  G4double path_;
  G4double time_;
  G4double edep_;
  G4ThreeVector pos_;
  G4ThreeVector mom_;
  G4int trackNo_;
  G4double xl_, yl_;
  std::vector <G4String> decayname_;

public:
  void SetPass() { pass_=1; }
  void SetPath(G4double path) { path_=path; }
  void SetLayerID(G4int id) { layerID_=id; }
  void SetWireID(G4int id) { wireID_=id; }
  void SetTime(G4double time) { time_=time; }
  void SetEdep(G4double edep) { edep_ = edep; }
  void AddEdep(G4double edep) { edep_ += edep; }
  void SetPos(const G4ThreeVector &pos) { pos_=pos; }
  void SetMom(const G4ThreeVector &mom) { mom_=mom; }
  void SetTrackNo(G4int no) { trackNo_=no; }
  void SetLocalPos(G4double x, G4double y) { xl_=x; yl_=y; }
  void SetDecayParticleName(G4String name) { decayname_.push_back(name); }

  G4int GetPass() const { return pass_; }
  G4double GetPath() const { return path_; }
  G4int GetLayerID() const { return layerID_; }
  G4int GetWireID() const { return wireID_; }
  G4double GetTime() const { return time_; }
  G4double GetEdep() const { return edep_; }
  const G4ThreeVector GetPos() const { return pos_; }
  const G4ThreeVector GetMom() const { return mom_; }
  G4int GetTrackNo() const { return trackNo_; }
  G4double GetXLocal() const { return xl_; }
  G4double GetYLocal() const { return yl_; }
  G4int NumOfHitParticles() { return decayname_.size(); }
  G4String GetDecayParticleName(int id) const { return decayname_[id]; }
};

//_____________________________________________________________________________
typedef G4THitsCollection<VPHit> VPHitsCollection;
extern G4Allocator<VPHit> VPHitAllocator;

//_____________________________________________________________________________
inline void *
VPHit::operator new(size_t)
{
  return static_cast<void *>(VPHitAllocator.MallocSingle());
}

//_____________________________________________________________________________
inline void
VPHit::operator delete(void* aHit)
{
  VPHitAllocator.FreeSingle(static_cast<VPHit*>(aHit));
}


#endif
