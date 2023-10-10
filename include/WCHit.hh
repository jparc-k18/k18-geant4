/*
  WCHit.hh
  2007/4  K.Shirotori
*/

#ifndef WCHit_h
#define WCHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4ThreeVector.hh"
#include "G4Allocator.hh"

class WCHit : public G4VHit
{
public:
  WCHit();
  ~WCHit() { }
private:
  WCHit( const WCHit & );
  WCHit & operator = ( const WCHit & );
public:
  int operator == ( const WCHit & ) const { return 0; }

  inline void * operator new ( size_t size );
  inline void operator delete( void * aHit );

//   void Draw() const;
//   void Print() const;

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

typedef G4THitsCollection<WCHit> WCHitsCollection;
extern G4Allocator<WCHit> WCHitAllocator;

inline void * WCHit::operator new( size_t )
{
  return static_cast<void *>( WCHitAllocator.MallocSingle() );
}

inline void WCHit::operator delete( void *aHit )
{
  WCHitAllocator.
    FreeSingle( static_cast<WCHit *>( aHit ) );
}
		       

#endif
