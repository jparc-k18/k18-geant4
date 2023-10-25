// -*- C++ -*-

#ifndef WC_HIT_HH
#define WC_HIT_HH 1

#include <G4VHit.hh>
#include <G4THitsCollection.hh>
#include <G4ThreeVector.hh>
#include <G4Allocator.hh>

#include "VHitInfo.hh"

class G4Step;

//_____________________________________________________________________________
class WCHit : public G4VHit, public VHitInfo
{
public:
  WCHit(const G4String& name, G4Step* step);
  virtual ~WCHit();

private:
  WCHit(const WCHit&);
  WCHit& operator =(const WCHit&);
public:
  void* operator new(size_t size);
  void operator delete(void* aHit);
  virtual void Draw();
  virtual void Print();

private:
  G4int pass_;
  G4int layerID_;
  G4int wireID_;
  G4double time_;
  G4double edep_;
  G4double npe_;
  G4ThreeVector pos_;
  G4ThreeVector mom_;
  G4int trackNo_;
  G4double xl_, yl_;
  std::vector <G4String> decayname_;

public:
  void SetPass( void ) { pass_=1; }
  void SetLayerID( G4int id ) { layerID_=id; }
  void SetWireID( G4int id ) { wireID_=id; }
  void SetTime( G4double time ) { time_=time; }
  void SetEdep( G4double edep ) { edep_ = edep; }
  void AddEdep( G4double edep ) { edep_ += edep; }
  void SetNPE( G4double npe )   { npe_  = npe; }
  void AddNPE( G4double npe )   { npe_ += npe; }
  void ClearNPE()   { npe_ = 0; }
  void SetPos( const G4ThreeVector &pos ) { pos_=pos; }
  void SetMom( const G4ThreeVector &mom ) { mom_=mom; }
  void SetTrackNo( G4int no ) { trackNo_=no; }
  void SetLocalPos( G4double x, G4double y ) { xl_=x; yl_=y; }
  void SetDecayParticleName( G4String name ) { decayname_.push_back( name ); }

  G4int GetPass() const { return pass_; }
  G4int GetLayerID( void ) const { return layerID_; }
  G4int GetWireID( void ) const { return wireID_; }
  G4double GetTime( void ) const { return time_; }
  G4double GetEdep() const { return edep_; }
  G4double GetNPE() const  { return npe_; }
  const G4ThreeVector GetPos( void ) const { return pos_; }
  const G4ThreeVector GetMom( void ) const { return mom_; }
  G4int GetTrackNo( void ) const { return trackNo_; }
  G4double GetXLocal( void ) const { return xl_; }
  G4double GetYLocal( void ) const { return yl_; }
  G4int NumOfHitParticles() { return decayname_.size(); }
  G4String GetDecayParticleName( int id ) const { return decayname_[id]; }
};

//_____________________________________________________________________________
using WCHitsCollection = G4THitsCollection<WCHit>;
extern G4Allocator<WCHit> WCHitAllocator;

//_____________________________________________________________________________
inline void*
WCHit::operator new(size_t)
{
  return static_cast<void*>(WCHitAllocator.MallocSingle());
}

//_____________________________________________________________________________
inline void
WCHit::operator delete(void* aHit)
{
  WCHitAllocator.FreeSingle(static_cast<WCHit*>(aHit));
}

#endif
