// -*- C++ -*-

#ifndef VP_HIT_HH
#define VP_HIT_HH

#include <G4VHit.hh>
#include <G4THitsCollection.hh>
#include <G4ThreeVector.hh>
#include <G4Allocator.hh>

#include "VHitInfo.hh"

//_____________________________________________________________________________
class VPHit : public G4VHit, public VHitInfo
{
public:
  VPHit(const G4String& name, G4Step* step);
  virtual ~VPHit();

private:
  VPHit(const VPHit&);
  VPHit& operator =(const VPHit&);

public:
  void* operator new(size_t size);
  void operator delete(void* aHit);
  virtual void Draw();
  virtual void Print();
};

//_____________________________________________________________________________
using VPHitsCollection = G4THitsCollection<VPHit>;
extern G4Allocator<VPHit> VPHitAllocator;

//_____________________________________________________________________________
inline void *
VPHit::operator new(size_t)
{
  return static_cast<void*>(VPHitAllocator.MallocSingle());
}

//_____________________________________________________________________________
inline void
VPHit::operator delete(void* aHit)
{
  VPHitAllocator.FreeSingle(static_cast<VPHit*>(aHit));
}


#endif
