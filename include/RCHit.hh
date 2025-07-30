// -*- C++ -*-

#ifndef RCHit_h
#define RCHit_h 1

#include <G4VHit.hh>
#include <G4THitsCollection.hh>
#include <G4ThreeVector.hh>
#include <G4Allocator.hh>

#include "VHitInfo.hh"

class G4Step;

//_____________________________________________________________________________
class RCHit : public G4VHit, public VHitInfo
{
public:
  RCHit(const G4String& name, G4Step* step);
  virtual ~RCHit();

private:
  RCHit(const RCHit &);
  RCHit& operator =(const RCHit&);

public:
  void* operator new(size_t size);
  void operator delete(void * aHit);
  virtual void Draw();
  virtual void Print();
};

//_____________________________________________________________________________
using RCHitsCollection = G4THitsCollection<RCHit>;
extern G4Allocator<RCHit> RCHitAllocator;

//_____________________________________________________________________________
inline void*
RCHit::operator new(size_t)
{
  return RCHitAllocator.MallocSingle();
}

//_____________________________________________________________________________
inline void
RCHit::operator delete(void *aHit)
{
  RCHitAllocator.FreeSingle(static_cast<RCHit*>(aHit));
}

#endif
