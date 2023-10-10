// -*- C++ -*-

#ifndef TOFHit_h
#define TOFHit_h 1

#include <G4VHit.hh>
#include <G4THitsCollection.hh>
#include <G4ThreeVector.hh>
#include <G4Allocator.hh>

#include "VHitInfo.hh"

class G4Step;

//_____________________________________________________________________________
class TOFHit : public G4VHit, public VHitInfo
{
public:
  TOFHit(const G4String& name, G4Step* step);
  virtual ~TOFHit();

private:
  TOFHit(const TOFHit &);
  TOFHit& operator =(const TOFHit&);

public:
  void* operator new(size_t size);
  void operator delete(void * aHit);
  virtual void Draw();
  virtual void Print();
};

using TOFHitsCollection = G4THitsCollection<TOFHit>;

extern G4Allocator<TOFHit> TOFHitAllocator;

//_____________________________________________________________________________
inline void*
TOFHit::operator new(size_t)
{
  return TOFHitAllocator.MallocSingle();
}

//_____________________________________________________________________________
inline void
TOFHit::operator delete(void *aHit)
{
  TOFHitAllocator.FreeSingle(static_cast<TOFHit*>(aHit));
}

#endif
