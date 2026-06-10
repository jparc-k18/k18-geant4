// -*- C++ -*-

#ifndef GE_HIT_HH
#define GE_HIT_HH 1

#include <G4Allocator.hh>
#include <G4THitsCollection.hh>
#include <G4VHit.hh>

#include "VHitInfo.hh"

class G4Step;

class GeHit : public G4VHit, public VHitInfo
{
public:
  GeHit(const G4String& name, G4Step* step);
  virtual ~GeHit();

private:
  GeHit(const GeHit&);
  GeHit& operator=(const GeHit&);

public:
  void* operator new(size_t size);
  void operator delete(void* hit);
  virtual void Draw();
  virtual void Print();
};

using GeHitsCollection = G4THitsCollection<GeHit>;
extern G4Allocator<GeHit> GeHitAllocator;

inline void*
GeHit::operator new(size_t)
{
  return GeHitAllocator.MallocSingle();
}

inline void
GeHit::operator delete(void* hit)
{
  GeHitAllocator.FreeSingle(static_cast<GeHit*>(hit));
}

#endif
