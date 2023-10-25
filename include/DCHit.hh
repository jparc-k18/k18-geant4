// -*- C++ -*-

#ifndef DC_HIT_HH
#define DC_HIT_HH 1

#include <G4VHit.hh>
#include <G4THitsCollection.hh>
#include <G4ThreeVector.hh>
#include <G4Allocator.hh>

#include "VHitInfo.hh"

class G4Step;

//_____________________________________________________________________________
class DCHit : public G4VHit, public VHitInfo
{
public:
  DCHit(const G4String& name, G4Step* step);
  virtual ~DCHit();
private:
  DCHit(const DCHit&);
  DCHit& operator =(const DCHit&);
public:
  void* operator new(size_t size);
  void operator delete(void* aHit);
  virtual void Draw();
  virtual void Print();
};

//_____________________________________________________________________________
using DCHitsCollection = G4THitsCollection<DCHit>;
extern G4Allocator<DCHit> DCHitAllocator;

//_____________________________________________________________________________
inline void*
DCHit::operator new(size_t)
{
  return static_cast<void*>(DCHitAllocator.MallocSingle());
}

//_____________________________________________________________________________
inline void
DCHit::operator delete(void* aHit)
{
  DCHitAllocator.FreeSingle(static_cast<DCHit*>(aHit));
}

#endif
